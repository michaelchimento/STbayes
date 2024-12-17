unloadNamespace("STbayes")
devtools::document()
devtools::install()
library(STbayes)
library(igraph)
library(dplyr)
library(NBDA)
library(ggplot2)
# Initialize final results dataframe
final_results <- data.frame(
    sim_id = integer(),
    model = character(),
    est_baseline = numeric(),
    est_s = numeric(),
    stringsAsFactors = FALSE
)

# Parameters
N <- 30  # Population size
k <- 4   # Degree of each node in the random regular graph
lambda_0 <- 0.0005
A <- 1    # Individual learning rate
s <- 5    # Social learning rate per unit connection
t_steps <- 1000
max_time <- t_steps + 1

# Run 100 simulations
num_simulations <- 100
for (sim in 21:100) {

    # Create random regular graph
    g <- sample_k_regular(N, k, directed = FALSE, multiple = FALSE)
    V(g)$name <- 1:N  # Name vertices

    # Initialize a dataframe to store the time of acquisition data
    df <- data.frame(id = 1:N, time = max_time, knowledgeable_associates = 0, max_time = max_time)

    # Simulate the diffusion
    for (t in 1:t_steps) {
        # Identify knowledgeable individuals
        learners <- df[df$time < max_time, "id"]

        # Identify naive individuals
        potential_learners <- c(1:N)
        potential_learners <- potential_learners[!(potential_learners %in% learners)]

        # Break the loop if no one left to learn
        if (length(potential_learners) == 0) break

        # Calculate the hazard
        learning_rates <- sapply(potential_learners, function(x) {
            neighbors <- neighbors(g, x)
            C <- sum(neighbors$name %in% learners)
            lambda <- lambda_0 * (A + s * C)
            return(lambda)
        })

        # Convert hazard to probability
        learning_probs <- 1 - exp(-learning_rates)

        # For each potential learner, determine whether they learn the behavior
        learners_this_step <- rbinom(n = length(potential_learners), size = 1, prob = learning_probs)

        # Update their time of acquisition
        new_learners <- potential_learners[learners_this_step == 1]
        df[df$id %in% new_learners, "time"] <- t
    }

    diffusion_data <- df %>%
        arrange(time) %>%
        group_by(time) %>%
        mutate(tie=ifelse(n()>1,1,0),
               seed=ifelse(time==0,1,0)) %>%
        ungroup() %>%
        select(-c(knowledgeable_associates))

    # Define the adjacency matrix
    adj_matrix <- as_adjacency_matrix(g, attr = NULL, sparse = FALSE)
    dim(adj_matrix) <- c(N, N, 1)

    # Define tie and seed vectors
    tie_vec = diffusion_data %>% arrange(time) %>% ungroup() %>% pull(tie)
    seed_vec = diffusion_data %>% arrange(id) %>% ungroup() %>% pull(seed)

    #### Fit NBDA model ####
    d <- nbdaData(
        label = "sim_data",
        assMatrix = adj_matrix,
        orderAcq = diffusion_data$id,
        timeAcq = diffusion_data$time,
        endTime = max_time,
        ties = tie_vec,
        demons = seed_vec
    )

    result <- tadaFit(d)
    result_NBDA <- data.frame(
        sim_id = sim,
        model = "NBDA",
        est_baseline = 1 / result@outputPar[1],
        est_s = result@outputPar[2]
    )

    # Add NBDA results to final results
    final_results <- rbind(final_results, result_NBDA)

    #### Fit STbayes model ####
    #import into STbayes
    diffusion_data <- diffusion_data %>%
        mutate(trial=1) %>%
        select(-c(tie, seed))
    edge_list <- as.data.frame(as_edgelist(g))
    names(edge_list) = c("from","to")
    edge_list$trial = 1
    edge_list$assoc = 1 #assign named edgeweight since this is just an edge list

    # Import data to STbayes
    data_list_user <- import_user_STb(diffusion_data, edge_list)
    model_obj <- generate_STb_model(data_list_user)

    # Fit model
    fit <- fit_STb(data_list_user, model_obj, chains = 4, cores = 4, iter = 2000, control = list(adapt_delta = 0.99))

    # Extract STbayes estimates
    STb_estimates <- STb_summary(fit, digits=5)
    transformed_s <- STb_estimates %>% filter(Parameter == "transformed_s") %>% pull(Mean)
    transformed_baseline <- STb_estimates %>% filter(Parameter == "transformed_baserate") %>% pull(Mean)

    result_STbayes <- data.frame(
        sim_id = sim,
        model = "STbayes",
        est_baseline = transformed_baseline,
        est_s = transformed_s
    )

    # Add STbayes results to final results
    final_results <- rbind(final_results, result_STbayes)

    # Print progress
    cat("Simulation", sim, "completed.\n")
}

# Save or view final results
print(final_results)

p1 = ggplot(final_results, aes(y=est_s, x=model))+
    geom_jitter()+
    stat_summary()+
    geom_hline(yintercept = 5, color="red") +
    labs(y="Est. S")

p2 = ggplot(final_results, aes(y=est_baseline, x=model))+
    geom_jitter()+
    stat_summary()+
    geom_hline(yintercept = .0005, color="red") +
    labs(y="Est. Baseline")

library(ggpubr)
ggarrange(p1,p2)
ggsave("../data/model_comparison.png", width=6, height=8, units="cm", scale=2)
