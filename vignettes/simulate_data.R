unloadNamespace("STbayes")
devtools::document()
devtools::install()
library(STbayes)
library(igraph)
library(dplyr)
library(NBDA)

# Parameters
N <- 30  # Population size
k <- 4    # Degree of each node in the random regular graph
lambda_0=0.001
A <- 1  # Individual learning rate
s <- 7  # Social learning rate per unit connection
t_steps <- 1000
max_time = t_steps+1

# create random regular graph
g <- sample_k_regular(N, k, directed = FALSE, multiple = FALSE)
V(g)$name <- 1:N  # Name vertices

# initialize a dataframe to store the time of acquisition data
df <- data.frame(id=1:N, time=max_time, knowledgeable_associates=0, max_time = max_time)

# one seed individual that already knows the behavior
#seed <- sample(1:N, 1)
#df[df$id == seed, c("time", "knowledgeable_associates")] <- c(0, 0)

# Simulate the diffusion
for (t in 1:t_steps) {
    # Identify knowledgeable individuals
    learners <- df[df$time < max_time, "id"]

    # Identify naive
    potential_learners <- c(1:N)
    potential_learners <- potential_learners[!(potential_learners %in% learners)]

    # break the loop if no one left to learn,
    if (length(potential_learners) == 0) break

    # calc the hazard
    learning_rates <- sapply(potential_learners, function(x) {
        neighbors <- neighbors(g, x)
        C <- sum(neighbors$name %in% learners)
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

    knowledgeable_associates <- sapply(1:N, function(x) {
        neighbors <- neighbors(g, x)
        C <- sum(neighbors$name %in% learners)
        return(C)
    })

    new_learner_associates <- knowledgeable_associates[match(new_learners, potential_learners)]
    df[df$id %in% new_learners, "knowledgeable_associates"] <- new_learner_associates
}

diffusion_data <- df %>%
    arrange(time) %>%
    group_by(time) %>%
    mutate(tie=ifelse(n()>1,1,0),
           seed=ifelse(time==0,1,0)) %>%
    ungroup() %>%
    select(-c(knowledgeable_associates))

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

#generate STAN model from input data
data_list_user = import_user_STb(diffusion_data, edge_list)
#generate STAN model from input data
model_obj = generate_STb_model(data_list_user)

# suggest writing to file for debugging
write(model_obj, file = "../inst/stan/model_from_simulate_data.stan")

# fit model
fit = fit_STb(data_list_user, "../inst/stan/model_from_simulate_data.stan", chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.99))

# check estimates
STb_summary(fit, digits=4)

# test import_NBDA_STb
data_list_nbda = import_NBDA_STb(d, network_names = c("assoc"))
model_obj = generate_STb_model(data_list_nbda)

#this should give the same results as fit above
fit = fit_STb(data_list_nbda, model_obj, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.99) )
STb_summary(fit)
