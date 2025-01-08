library(STbayes)
library(igraph)
library(dplyr)
library(NBDA)
library(ggplot2)
# Parameters
N <- 50  # Population size
k <- 4    # Degree of each node in the random regular graph
lambda_0=0.001 #baseline
A <- 1  # Individual learning rate
s <- 5  # Social learning rate per unit connection
t_steps <- 1000
max_time = t_steps+1

# create random regular graph
g <- sample_k_regular(N, k, directed = FALSE, multiple = FALSE)
V(g)$name <- 1:N

# initialize a dataframe to store the time of acquisition data
df <- data.frame(id=1:N, time=max_time, max_time = max_time)

# If you want to set a demonstrator, uncomment below
# seed <- sample(1:N, 1)
# df[df$id == seed, c("time", "informed_associates")] <- c(0, 0)

# simulate the diffusion
for (t in 1:t_steps) {
    # identify knowledgeable individuals
    informed <- df[df$time < max_time, "id"]

    # identify naive
    potential_learners <- c(1:N)
    potential_learners <- potential_learners[!(potential_learners %in% informed)]

    # break the loop if no one left to learn,
    if (length(potential_learners) == 0) break

    # calc the hazard
    learning_rates <- sapply(potential_learners, function(x) {
        neighbors <- neighbors(g, x)
        C <- sum(neighbors$name %in% informed)
        lambda <- lambda_0 * (A + s*C)
        return(lambda)
    })

    # convert hazard to probability
    learning_probs <- 1 - exp(-learning_rates)

    # for each potential learner, determine whether they learn the behavior
    learners_this_step <- rbinom(n=length(potential_learners), size=1, prob=learning_probs)

    # update their time of acquisition
    new_learners <- potential_learners[learners_this_step == 1]
    df[df$id %in% new_learners, "time"] <- t
}

diffusion_data <- df %>%
    arrange(time) %>%
    group_by(time, .drop = T) %>%
    mutate(tie=ifelse(n()>1,1,0),
           seed=ifelse(time==0,1,0))

hist(diffusion_data$time)

# Define the adjacency matrix
adj_matrix <- as_adjacency_matrix(g, attr=NULL, sparse=FALSE)
dim(adj_matrix) = c(N,N,1)

tie_vec = diffusion_data %>% arrange(time) %>% ungroup() %>% select(tie)
seed_vec = diffusion_data %>% arrange(id) %>% ungroup() %>% select(seed)

#### Fit NBDA model ####
d = nbdaData(label="sim_data",
             assMatrix = adj_matrix,
             orderAcq = diffusion_data$id,
             timeAcq = diffusion_data$time,
             endTime = max_time,
             ties = tie_vec$tie,
             demons = seed_vec$seed)
result = tadaFit(d)
data.frame(Variable=result@varNames,MLE=result@outputPar,SE=result@se)
est_rate = 1/result@outputPar[1]
est_rate

#import into STbayes
diffusion_data <- diffusion_data %>%
    mutate(trial=1) %>%
    select(-c(tie, seed))
edge_list <- as.data.frame(as_edgelist(g))
names(edge_list) = c("from","to")
edge_list$trial = 1
edge_list$assoc = 1 #assign named edgeweight since this is just an edge list

#save(diffusion_data, file="../data/example_diffusion_data.rda")
#save(edge_list, file="../data/example_edge_list.rda")

#generate STAN model from input data
data_list_user = import_user_STb(diffusion_data, edge_list)

#generate STAN model from input data
model_obj = generate_STb_model(data_list_user, gq=T, est_acqTime = T)

# Write to file for debugging? uncomment below why not
#write(model_obj, file = "../data/model_from_simulate_data.stan")

# fit model
fit = fit_STb(data_list_user, model_obj, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.99))

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
