library(STbayes)
library(igraph)
library(dplyr)
library(NBDA)
library(ggplot2)

# Parameters
N <- 100
k <- 7
lambda_0 <- 0.001
A <- 1
s <- 3
t_max <- 1000

g <- sample_k_regular(N, k, directed = FALSE, multiple = FALSE)
V(g)$name <- 1:N
df <- data.frame(id=1:N, time=Inf)

# Optionally seed one demonstrator
# df[df$id == sample(1:N, 1), "time"] <- 0

time <- 0
while (time < t_max) {
    # determine informed / naive
    informed <- df[df$time < Inf, "id"]
    naive <- setdiff(1:N, informed)
    if (length(naive) == 0) break  # diffusion done

    # compute hazard for each naive individual
    hazards <- sapply(naive, function(x) {
        neighbors_x <- neighbors(g, x)
        C <- sum(neighbors_x$name %in% informed)
        lambda_0 * (A + s * C)
    })

    total_hazard <- sum(hazards)
    if (total_hazard == 0) break  # no possible learning

    # draw time until next event
    wait_time <- rexp(1, rate = total_hazard)
    time <- time + wait_time
    if (time > t_max) break

    # choose *who* learns based on hazard-weighted sampling
    learner_idx <- sample(1:length(naive), size = 1, prob = hazards)
    learner <- naive[learner_idx]
    df[df$id == learner, "time"] <- time
}

df$time[is.infinite(df$time)] <- t_max

event_data <- df %>%
    arrange(time) %>%
    group_by(time, .drop = T) %>%
    mutate(tie=ifelse(n()>1,1,0),
           seed=ifelse(time==0,1,0))

event_data <- event_data %>%
    mutate(t_end=max(event_data$time)+1)

hist(event_data$time)

# Define the adjacency matrix
adj_matrix <- as_adjacency_matrix(g, attr=NULL, sparse=FALSE)
dim(adj_matrix) = c(N,N,1)

tie_vec = event_data %>% arrange(time) %>% ungroup() %>% select(tie)
seed_vec = event_data %>% arrange(id) %>% ungroup() %>% select(seed)

#### Fit NBDA model ####
d = nbdaData(label="sim_data",
             assMatrix = adj_matrix,
             orderAcq = event_data$id,
             timeAcq = event_data$time,
             endTime = t_steps+1,
             ties = tie_vec$tie,
             demons = seed_vec$seed)
result = tadaFit(d)
data.frame(Variable=result@varNames,MLE=result@outputPar,SE=result@se)
est_rate = 1/result@outputPar[1]
est_rate

#import into STbayes
event_data <- event_data %>%
    mutate(trial=1) %>%
    select(-c(tie, seed))
edge_list <- as.data.frame(as_edgelist(g))
names(edge_list) = c("from","to")
edge_list$trial = 1
edge_list$assoc = 1 #assign named edgeweight since this is just an edge list

#save(event_data, file="../data/example_event_data.rda")
#save(edge_list, file="../data/example_edge_list.rda")

#generate STAN model from input data
data_list_user = import_user_STb(event_data, edge_list)
hist(data_list_user$D)

#generate STAN model from input data
model_obj = generate_STb_model(data_list_user, gq=T, est_acqTime = T)

# fit model
fit = fit_STb(data_list_user, "../inst/extdata/STAN_example_vanilla_ctada.stan", chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.9))

# check estimates
STb_summary(fit, digits=4)
#get data for estimated times
acqdata = extract_acqTime(fit, data_list_user)

#plot estimated times versus observed times w/ HPD
ggplot(acqdata, aes(x = observed_time, y = mean_time)) +
    geom_abline(slope = 1, intercept = 0,             # Line y = x
                color = "red", linetype = "dotted", linewidth = 1) +
    geom_pointrange(aes(ymin = lower_hpd, ymax = upper_hpd), size=.8) +
    facet_wrap(~trial, scales = "free_x") +
    labs(
        title = "Estimated acquisition time with 95% HPD Intervals",
        x = "Observed time",
        y = "Estimated time"
    ) +
    theme_minimal()

#plot estimated times versus observed times w/ residuals
ggplot(acqdata, aes(x = observed_time, y = mean_time)) +
    geom_segment(
        aes(x = observed_time, xend = observed_time, y = mean_time, yend = observed_time), # connect predicted to slope line
        color = "red",
        alpha = 0.2
    ) +
    geom_point(alpha = 0.6, size=2) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    facet_wrap(~trial, scales = "free_x") +
    labs(
        title = "Estimated acquisition time with residuals",
        x = "Observed time",
        y = "Estimated time"
    ) +
    theme_minimal()


# why not import from NBDA object? It works the same
data_list_nbda = import_NBDA_STb(d, network_names = c("assoc"))
model_obj = generate_STb_model(data_list_nbda, gq=T, est_acqTime = T)
fit = fit_STb(data_list_nbda, model_obj, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.99) )
STb_summary(fit)

