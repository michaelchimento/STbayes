library(STbayes)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Parameters
set.seed(123)
N <- 20 # Population size
k <- 3 # Degree of each node in the random regular graph
log_lambda_0_mean <- -5
lambda_0_sd <- .25
log_sprime_mean <- -3.5
s_sd <- .25
t_steps <- 1000
num_trials <- 100

# Assign individual parameters
individual_params <- data.frame(
  id = 1:N,
  lambda_0 = exp(log_lambda_0_mean + rnorm(N, mean = 0, sd = lambda_0_sd)),
  s_prime = exp(log_sprime_mean + rnorm(N, mean = 0, sd = s_sd))
) %>%
  mutate(s = s_prime / lambda_0)
individual_params %>% summarize(mean(lambda_0), mean(s_prime), mean(s))

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
  plot(g)
  V(g)$name <- 1:N # Name vertices

  # Initialize a dataframe to store the time of acquisition data
  df <- data.frame(id = 1:N, time = t_steps + 1, t_end = t_steps)

  # Simulate the diffusion
  for (t in 1:t_steps) {
    # Identify knowledgeable individuals
    informed <- df[df$time <= t_steps, "id"]

    # identify naive individuals
    potential_learners <- c(1:N)
    potential_learners <- potential_learners[!(potential_learners %in% informed)]

    # break the loop if no one left to learn
    if (length(potential_learners) == 0) {
      df$t_end <- t - 1
      break
    }

    # Calculate the hazard
    learning_rates <- sapply(potential_learners, function(x) {
      neighbors <- neighbors(g, x)
      C <- sum(neighbors$name %in% informed)
      lambda <- individual_params[individual_params$id == x, "lambda_0"] +
        (individual_params[individual_params$id == x, "s_prime"] * C)
      return(lambda)
    })

    # convert hazard to probability
    learning_probs <- 1 - exp(-learning_rates)

    # for each potential learner, determine whether they learn the behavior
    learners_this_step <- rbinom(n = length(potential_learners), size = 1, prob = learning_probs)

    # update their time of acquisition
    new_learners <- potential_learners[learners_this_step == 1]
    df[df$id %in% new_learners, "time"] <- t
  }


  event_data <- df %>%
    arrange(time) %>%
    group_by(time, .drop = T) %>%
    mutate(
      tie = ifelse(n() > 1, 1, 0),
      seed = ifelse(time == 0, 1, 0)
    )

  #### Create trial data and append to dataframe
  event_data <- event_data %>%
    mutate(trial = trial) %>%
    select(-c(tie, seed))
  edge_list <- as.data.frame(as_edgelist(g))
  names(edge_list) <- c("from", "to")
  edge_list$trial <- trial
  edge_list$assoc <- 1 # assign named edgeweight since this is just an edge list

  # add data
  all_event_data <- rbind(all_event_data, event_data)
  all_edge_list <- rbind(all_edge_list, edge_list)
  message("Trial", trial, "completed.\n")
}
hist(all_event_data$t_end)
# import data to STbayes and fit model
data_list_user <- import_user_STb(all_event_data, all_edge_list)
model_obj <- generate_STb_model(data_list_user, veff = c("lambda_0", "s"))

# write(model_obj, file = "../inst/extdata/STAN_example_veff.stan")
fit_veff <- fit_STb(data_list_user,
  "../inst/extdata/STAN_example_veff.stan",
  chains = 5,
  cores = 5,
  parallel_chains = 5,
  iter = 5000,
  control = list(adapt_delta = 0.99)
)

#### Plot mean estimates ####
summary_df <- STb_summary(fit_veff, digits = 5, depth = 1)

# True values to compare
true_vals <- c(
    log_lambda_0_mean = -5,
    log_s_mean = -3.5,
    `sigma_ID[1]` = 0.25,
    `sigma_ID[2]` = 0.25
)
param_order <- names(true_vals)

selected <- summary_df %>%
    filter(Parameter %in% param_order) %>%
    mutate(
        True = true_vals[Parameter],
        Label = factor(Parameter, levels = rev(param_order))
    )

label_exprs <- rev(c(
    expression(log(lambda[0])),
    expression(log(s * "'")),
    expression(sigma[id * ", " * lambda[0]]),
    expression(sigma[id * ", " * s * "'"])
))

g1 <- ggplot(selected, aes(y = Label, x = Median)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = HPDI_Lower, xmax = HPDI_Upper), height = 0.2) +
    geom_point(aes(x = True), shape = 4, size = 3, color = "steelblue") +
    scale_y_discrete(labels = label_exprs) +
    labs(
        x = "Value",
        y = NULL
    ) +
    theme_linedraw(base_size = 14)

#### Plot varying effects against true vals ####
draws <- fit_veff$draws(variables = c("lambda_0", "s_prime"), format = "draws_df")
# Identify number of individuals (assuming parameters are vectorized per individual)
lambda_names <- grep("^lambda_0\\[", colnames(draws), value = TRUE)
s_names <- grep("^s_prime\\[", colnames(draws), value = TRUE)

# Compute means
lambda_0_mean_estimates <- colMeans(draws[, lambda_names, drop = FALSE])
s_mean_estimates <- colMeans(draws[, s_names, drop = FALSE])

N <- length(lambda_0_mean_estimates)
estimated <- data.frame(
  id = seq_len(N),
  est_lambda_0 = lambda_0_mean_estimates,
  est_s_prime = s_mean_estimates
)

# Join with known/true parameter values
df_comparison <- dplyr::left_join(estimated, individual_params, by = "id")

p1 <- ggplot(df_comparison, aes(x = s_prime, y = est_s_prime)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = "dashed") +
  geom_segment(
    aes(x = s_prime, xend = s_prime, y = est_s_prime, yend = s_prime),
    color = "red",
    alpha = 0.4
  ) +
  labs(
    x = expression("Real " * s * "'"[i]),
    y = expression("Estimated " * s * "'"[i])
  ) +
  theme_linedraw()
p2 <- ggplot(df_comparison, aes(x = lambda_0, y = est_lambda_0)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = "dashed") +
  geom_segment(
    aes(x = lambda_0, xend = lambda_0, y = est_lambda_0, yend = lambda_0),
    color = "red",
    alpha = 0.4
  ) +
  labs(
    x = expression("Real " * lambda[0 * i]),
    y = expression("Estimated " * lambda[0 * i])
  ) +
  theme_linedraw()
g2 <- ggarrange(p1, p2, labels = c("B","C"))


ggarrange(g1, g2, nrow = 2, labels = "A")
ggsave("../docs/ID_veff.pdf", width = 10, height = 6, units = "cm", scale = 2)
