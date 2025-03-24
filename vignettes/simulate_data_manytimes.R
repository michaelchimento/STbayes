library(STbayes)
library(igraph)
library(dplyr)
library(NBDA)
library(ggplot2)
library(ggpubr)

# Initialize final results dataframe
final_results <- data.frame(
    sim_id = integer(),
    model = character(),
    mean_lambda_0 = numeric(),
    mean_s = numeric(),
    median_lambda_0 = numeric(),
    median_s = numeric(),
    stringsAsFactors = FALSE
)

# Parameters
N <- 30  # Population size
k <- 3   # Degree of each node in the random regular graph
lambda_0 <- 0.001
A <- 1    # Individual learning rate
s <- 5    # Social learning rate per unit connection
t_steps <- 1000

# Run 100 simulations
num_simulations <- 100
for (sim in 1:num_simulations) {
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

    event_data <- df %>%
        arrange(time) %>%
        group_by(time, .drop = T) %>%
        mutate(tie=ifelse(n()>1,1,0),
               seed=ifelse(time==0,1,0)) %>%
        ungroup() %>%
        mutate(tie_nbda = if_else(lag(time)==time,1,0, missing = 0))

    # Define the adjacency matrix
    adj_matrix <- as_adjacency_matrix(g, attr = NULL, sparse = FALSE)
    dim(adj_matrix) <- c(N, N, 1)

    # Define tie and seed vectors
    tie_vec = event_data %>% arrange(time) %>% ungroup() %>% pull(tie_nbda)
    seed_vec = event_data %>% arrange(id) %>% ungroup() %>% pull(seed)

    #### Fit NBDA model ####
    d <- nbdaData(
        label = "sim_data",
        assMatrix = adj_matrix,
        orderAcq = event_data$id,
        timeAcq = event_data$time,
        endTime = t_steps+1,
        ties = tie_vec
    )

    result <- tadaFit(d)
    result_NBDA <- data.frame(
        sim_id = sim,
        model = "NBDA",
        mean_lambda_0 = 1 / result@outputPar[1],
        mean_s = result@outputPar[2],
        median_lambda_0 = 1 / result@outputPar[1],
        median_s = result@outputPar[2]
    )
    final_results <- rbind(final_results, result_NBDA)

    #### Fit STbayes model ####
    event_data <- event_data %>%
        mutate(trial=1) %>%
        select(-c(tie, seed))
    edge_list <- as.data.frame(as_edgelist(g))
    names(edge_list) = c("from","to")
    edge_list$trial = 1
    edge_list$assoc = 1 #assign named edgeweight since this is just an edge list

    # import data to STbayes
    data_list_user <- import_user_STb(event_data, edge_list)
    model_obj <- generate_STb_model(data_list_user, gq=F)
    fit <- fit_STb(data_list_user, model_obj, chains = 5, cores = 5, iter = 5000, control = list(adapt_delta = 0.99))
    STb_estimates <- STb_summary(fit, digits=5)

    result_STbayes <- data.frame(
        sim_id = sim,
        model = "STbayes",
        median_lambda_0 = STb_estimates %>% filter(Parameter == "lambda_0") %>% pull(Median),
        mean_lambda_0 = STb_estimates %>% filter(Parameter == "lambda_0") %>% pull(Mean),
        median_s = STb_estimates %>% filter(Parameter == "s") %>% pull(Median),
        mean_s = STb_estimates %>% filter(Parameter == "s") %>% pull(Mean)
    )

    # add STbayes results to final results
    final_results <- rbind(final_results, result_STbayes)

    message("Simulation", sim, "completed.\n")
}

ggplot(final_results, aes(y=median_lambda_0, x=median_s))+
    facet_wrap(~model)+
    geom_point(alpha=0.4)+
    geom_density_2d()+
    geom_hline(yintercept = lambda_0, color="red", linetype="dashed") +
    geom_vline(xintercept = s, color="red", linetype="dashed") +
    labs(
        y = expression("Estimated " * lambda[0] ),
        x = expression("Estimated " * s )
    ) +
    theme_bw()
ggsave(file="../docs/Fig2_NBDA_STbayes_estimates.png", width=10, height=4, units="cm", scale=2)

final_results_wide <- final_results %>%
    select(sim_id, model, median_lambda_0, median_s) %>%
    tidyr::pivot_wider(
        names_from = model,
        values_from = c(median_lambda_0, median_s),
        names_glue = "{model}_{.value}",
        values_fn = mean
    )

p1 = ggplot(final_results_wide, aes(y=NBDA_median_lambda_0, x=STbayes_median_lambda_0))+
    geom_jitter()+
    labs(y = expression("NBDA estimated " * lambda[0] ),
         x = expression("STbayes estimated " * lambda[0] ))+
    theme_bw()

p2 = ggplot(final_results_wide, aes(y=NBDA_median_s, x=STbayes_median_s))+
    geom_jitter()+
    labs(y = expression("NBDA estimated " * s ),
         x = expression("STbayes estimated " * s ))+
    theme_bw()

ggarrange(p1,p2, labels = c("A","B"))
ggsave(file="../docs/FigS_NBDA_STbayes_correlations.png", width=10, height=4, units="cm", scale=2)

ggplot(final_results, aes(y=log(median_lambda_0), x=log(median_s)))+
    facet_wrap(~model)+
    geom_point(alpha=0.2)+
    geom_density_2d()+
    geom_hline(yintercept = log(lambda_0), color="red") +
    geom_vline(xintercept = log(s), color="red") +
    labs(y="Est. Baseline", x="Est_s")



