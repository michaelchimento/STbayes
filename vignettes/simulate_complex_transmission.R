library(STbayes)
library(igraph)
library(dplyr)
library(NBDA)
library(ggplot2)

dini_func <- function (x, k){
    x = 2*x - 1
    y = ((x - k * x) / (k - 2 * k * abs(x) + 1) + 1) / 2
    return(y)
}

# Parameters
N <- 50  # Population size
n_trials <- 10
k <- 5   # Degree of each node
log_lambda_0_mean <- 7  # baseline
lambda_0 <- 1 / exp(log_lambda_0_mean)
s <- 4
s_prime <- lambda_0 * s
A <- 1   # Individual learning rate
k_shape <- -0.9
t_steps <- 3000
# run multiple trials
trial_results <- vector("list", n_trials)

# create random regular graph
g <- sample_k_regular(N, k, directed = FALSE, multiple = FALSE)
V(g)$name <- 1:N

for (trial in 1:n_trials) {


    # initialize df for this trial
    df <- data.frame(id = 1:N, time = t_steps + 1, t_end = t_steps)

    # uncomment below if you want a random seed informed individual each trial
    # seed <- sample(1:N, 1)
    # df[df$id == seed, "time"] <- 0

    # simulate diffusion
    for (t in 1:t_steps) {
        # get informed
        informed <- df[df$time <= t_steps, "id"]

        # get naive
        potential_learners <- setdiff(1:N, informed)

        # break if everyone knows
        if (length(potential_learners) == 0) break

        # calc hazard for each naive
        learning_rates <- sapply(potential_learners, function(x) {
            neighbors <- neighbors(g, x)
            C <- sum(neighbors$name %in% informed)
            prop_know <- C / k
            s_term <- dini_func(prop_know, k_shape)
            lambda <- lambda_0 * A + s_prime * s_term
            return(lambda)
        })

        # hazard to prob
        learning_probs <- 1 - exp(-learning_rates)

        # draw who learns
        learners_this_step <- rbinom(n = length(potential_learners), size = 1, prob = learning_probs)

        # update time of acquisition
        new_learners <- potential_learners[learners_this_step == 1]
        df[df$id %in% new_learners, "time"] <- t
    }

    # store trial
    trial_results[[trial]] <- df
}

# after the simulation loop
for (trial in 1:n_trials) {
    trial_results[[trial]]$trial <- trial
}

# combine all trials into one dataframe
combined_df <- bind_rows(trial_results)

# check the result
head(combined_df)

event_data <- combined_df %>%
    arrange(trial, time) %>%
    group_by(trial, time, .drop = T) %>%
    mutate(tie=ifelse(n()>1,1,0),
           seed=ifelse(time==0,1,0))

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
event_data <- event_data  %>%
    select(-c(tie, seed))
edge_list <- as.data.frame(as_edgelist(g))
names(edge_list) = c("from","to")
edge_list$assoc = 1 #assign named edge weight since this is just an edge list

# now replicate the edgelist for each trial and tag trial number
edge_list_trials <- bind_rows(lapply(1:n_trials, function(trial) {
    edge_list %>% mutate(trial = trial)
}))

#save(event_data, file="../data/example_event_data.rda")
#save(edge_list, file="../data/example_edge_list.rda")

#generate STAN model from input data
data_list_user = import_user_STb(event_data, edge_list_trials)

#generate STAN model from input data
model_obj = generate_STb_model(data_list_user, gq=T, est_acqTime = F, transmission_func = "standard")
model_obj1 = generate_STb_model(data_list_user, gq=T, est_acqTime = F, transmission_func = "freq-dep1")
model_obj2 = generate_STb_model(data_list_user, gq=T, est_acqTime = F, transmission_func = "freq-dep2")

# Write to file for debugging? uncomment below why not
write(model_obj1, file = "../inst/extdata/STAN_example_complex_f_transmission.stan")
write(model_obj2, file = "../inst/extdata/STAN_example_complex_k_transmission.stan")

# fit model
fit_simple = fit_STb(data_list_user, model_obj, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.95))
STb_summary(fit_simple)
fit_complex = fit_STb(data_list_user, model_obj2, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.95))
STb_summary(fit_complex)

fit_complex2 = fit_STb(data_list_user, model_obj1, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.95))
STb_summary(fit_complex2)

loo_output = STb_compare(fit_simple, fit_complex, fit_complex2)
# Extract comparison results for plotting
comparison_df <- as.data.frame(loo_output$comparison)
comparison_df$model <- rownames(comparison_df)
comparison_df$model <- as.factor(comparison_df$model)
levels(comparison_df$model) = c("complex 'k'","complex 'f'","simple")
# Plot comparison
# Best predictive model will be at top. If SE of difference crosses zero,
# you should be less certain that the best fit model actually provides a better fit
ggplot(comparison_df, aes(x = reorder(model, elpd_diff), y = elpd_diff)) +
    geom_point(size = 3) +  #mean elpd difference
    geom_errorbar(aes(ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff), width = 0.2) + #SE of elpd diff
    coord_flip() +
    theme_minimal() +
    labs(x = "Model", y = "ELPD LOO Difference") +
    theme(axis.text = element_text(size = 12))
ggsave("../docs/FigS_complex_model_comparison.png", width=7, height=6, scale=2, units="cm")
# Plot comparison
# Best predictive model will be at top. If SE of difference crosses zero,
# you should be less certain that the best fit model actually provides a better fit
ggplot(comparison_df, aes(x = reorder(model, looic), y = looic)) +
    geom_point(size = 3) +  #mean elpd difference
    geom_errorbar(aes(ymin = looic - se_looic, ymax = looic + se_looic), width = 0.2) + #SE of elpd diff
    coord_flip() +
    theme_minimal() +
    labs(x = "Model", y = "LOOIC", title = "Model Comparison") +
    theme(axis.text = element_text(size = 12))

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
