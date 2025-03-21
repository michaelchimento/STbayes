library(STbayes)
library(igraph)
library(dplyr)
library(ggplot2)

# Parameters
N <- 20  # Population size
k <- 4   # Degree of each node in the random regular graph
lambda_0_mean <- 8
lambda_0_sd <- .5
s_mean <- 2
s_sd <- .5
t_steps <- 1000
num_trials <- 25

# Assign individual parameters
individual_params <- data.frame(
    id = 1:N,
    lambda_0 = 1/exp(lambda_0_mean - rnorm(N, mean = 0, sd = lambda_0_sd)),
    s = exp(s_mean - rnorm(N, mean = 0, sd = s_sd))
)
individual_params %>% summarize(mean(lambda_0), mean(s))
plot(density(individual_params$lambda_0))
# storage for all trials
all_event_data <- data.frame(
    id = integer(),
    time = numeric(),
    t_end = numeric(),
    trial = integer()
)
all_edge_list <- data.frame(
    from = integer(),
    to = integer(),
    trial = integer(),
    assoc = numeric()
)

for (trial in 1:num_trials) {

    # Create random regular graph
    g <- sample_k_regular(N, k, directed = FALSE, multiple = FALSE)
    V(g)$name <- 1:N  # Name vertices

    # Initialize a dataframe to store the time of acquisition data
    df <- data.frame(id = 1:N, time = t_steps+1, t_end = t_steps)

    # Simulate the diffusion
    for (t in 1:t_steps) {
        # Identify knowledgeable individuals
        informed <- df[df$time <= t_steps, "id"]

        # identify naive individuals
        potential_learners <- c(1:N)
        potential_learners <- potential_learners[!(potential_learners %in% informed)]

        # break the loop if no one left to learn
        if (length(potential_learners) == 0) break

        # Calculate the hazard
        learning_rates <- sapply(potential_learners, function(x) {
            neighbors <- neighbors(g, x)
            C <- sum(neighbors$name %in% informed)
            lambda <- individual_params[individual_params$id == x, "lambda_0"] *
                (1 + individual_params[individual_params$id == x, "s"] * C)
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

    event_data <- df %>%
        arrange(time) %>%
        group_by(time, .drop = T) %>%
        mutate(tie=ifelse(n()>1,1,0),
               seed=ifelse(time==0,1,0))

    #### Create trial data and append to dataframe
    event_data <- event_data %>%
        mutate(trial=trial) %>%
        select(-c(tie, seed))
    edge_list <- as.data.frame(as_edgelist(g))
    names(edge_list) = c("from","to")
    edge_list$trial = trial
    edge_list$assoc = 1 #assign named edgeweight since this is just an edge list

    # add data
    all_event_data <- rbind(all_event_data, event_data)
    all_edge_list <- rbind(all_edge_list, edge_list)
    message("Trial", trial, "completed.\n")
}

# import data to STbayes and fit model
data_list_user <- import_user_STb(all_event_data, all_edge_list)
model_obj <- generate_STb_model(data_list_user, veff=c("lambda_0", "s"), priors = list(log_s="normal(2,3)"))
write(model_obj, file = "../data/STAN_example_veff.stan")
fit <- fit_STb(data_list_user, model_obj, chains = 5, cores = 5, iter = 5000, control = list(adapt_delta = 0.99))

#check to see if output matches simulation inputs
STb_summary(fit, digits=5, depth=1)

# extract posterior samples from the fit and calculate individual estimates
posterior_samples <- rstan::extract(fit)
lambda_0_estimates <- posterior_samples$lambda_0
s_estimates <- posterior_samples$s

#each col is for 1 ind
lambda_0_mean_estimates <- colMeans(lambda_0_estimates)
s_mean_estimates <- colMeans(s_estimates)

estimated <- data.frame(
    id = 1:N,
    est_lambda_0 = lambda_0_mean_estimates,
    est_s = s_mean_estimates
)

# join back with predefined parameter settings and reshape into long format for plotting
df_comparison = left_join(estimated, individual_params)
df_comparison <- reshape(
    df_comparison,
    varying = list(c("est_lambda_0", "est_s"), c("lambda_0", "s")),
    v.names = c("estimated", "real"),
    timevar = "parameter",
    times = c("lambda_0", "s"),
    direction = "long"
)
rownames(df_comparison) <- NULL

#plot real vs estimated
ggplot(df_comparison, aes(x=real, y=estimated))+
    facet_wrap(~parameter, scales = "free")+
    geom_point()+
    geom_abline(slope = 1, intercept=0, lty="dashed") +
    labs(title="Real versus estimated parameter values for individuals", x="Real value", y="Estimated value")

ggsave("../docs/ID_veff_real_vs_estimated.png", width=12, height=6, units="cm", scale=2)
