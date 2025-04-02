//stan
data {
    int<lower=0> K;                // Number of trials
    int<lower=0> Q;                // Number of individuals in each trial
    int<lower=1> P;                // Number of unique individuals
    int<lower=1> S;              // number of posterior samples for A edge weights
    array[K] int<lower=0> N;       // Number of individuals that learned during observation period
    array[K] int<lower=0> N_c;     // Number of right-censored individuals
    array[K, Q] int<lower=-1> ind_id; // IDs of individuals
    array[K] int<lower=1> T;       // Maximum time periods
    int<lower=1> T_max;            // Max timesteps reached
    array[K,P] int<lower=-1> t;     // Time of acquisition for each individual
    array[K, T_max] real<lower=0> D; // Scaled durations
    array[K, T_max, S] matrix[P, P] A_value; // Network matrices
    array[K] matrix[T_max, P] Z;   // Knowledge state slash cue matrix
    int<lower=0> N_veff;
}
parameters {
    real log_lambda_0_mean;  // Log baseline learning rate
    real log_s_mean; // Overall social transmission rate
}
transformed parameters {
   real<lower=0> lambda_0 = 1 / exp(log_lambda_0_mean);
real<lower=0> s_prime = 1 / exp(log_s_mean);
real<lower=0> s = s_prime/lambda_0;
}
model {
    log_lambda_0_mean ~ uniform(-10, 10);
    log_s_mean ~ uniform(-10, 10);
    vector[S] log_likelihoods;  // store likelihoods over samples
    for (trial in 1:K) {
        for (n in 1:N[trial]) {

            for (d in 1:S) {  // loop over posterior draws
                log_likelihoods[d] = 0;  // initialize log-likelihood for this draw
                int id = ind_id[trial, n];
                int learn_time = t[trial, id];
                if (learn_time > 0) {
                    for (time_step in 1:learn_time) {
                        real ind_term = 1.0;
                        real soc_term = s_prime * (sum(A_value[trial, time_step , d][id, ] .* Z[trial][time_step, ])) ;
                        real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                        log_likelihoods[d] += -lambda;
                        if (time_step == learn_time) {
                            log_likelihoods[d] += log( (lambda_0 * ind_term + soc_term));
                        }
                    }
                }
            } 
            target += log_sum_exp(log_likelihoods) - log(S);
        }
        if (N_c[trial] > 0) {
            for (c in 1:N_c[trial]) {
            for (d in 1:S) {  // loop over posterior draws
                log_likelihoods[d] = 0;  // initialize log-likelihood for this draw
                int id = ind_id[trial, N[trial] + c];
                for (time_step in 1:T[trial]) {
                    real ind_term = 1.0;
                    real soc_term = s_prime * (sum(A_value[trial, time_step , d][id, ] .* Z[trial][time_step, ])) ;
                    real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                    log_likelihoods[d] += -lambda;
                }
            }
            target += log_sum_exp(log_likelihoods) - log(S);
        }
       } 
    }
}

generated quantities {
    matrix[K, Q] log_lik_matrix = rep_matrix(0.0, K, Q);           // LL for each observation
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            vector[S] log_likelihoods;  // store likelihoods over samples
                for (d in 1:S) {  // loop over posterior draws
                log_likelihoods[d] = 0;  // initialize log-likelihood for this draw
            if (learn_time > 0){
                real cum_hazard = 0; //set val before adding
                for (time_step in 1:learn_time) {
                    real ind_term = 1.0;
                    real soc_term = s_prime * (sum(A_value[trial, time_step , d][id, ] .* Z[trial][time_step, ])) ;
                    real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                    cum_hazard += lambda; // accumulate hazard
                    //if it learn_time, record the ll
                    if (time_step == learn_time){
                        log_likelihoods[d] += log( (lambda_0 * ind_term + soc_term));
                    }
                }
            }
        }
           log_lik_matrix[trial, n] = log_sum_exp(log_likelihoods) - log(S);
        }
        // Contributions of censored individuals
        if (N_c[trial] > 0) {
            for (c in 1:N_c[trial]) {
                int id = ind_id[trial, N[trial] + c];
                int censor_time = T[trial]; // Censoring time (end of observation)
                vector[S] log_likelihoods;  // store likelihoods over samples
                    for (d in 1:S) {  // loop over posterior draws
                    log_likelihoods[d] = 0;  // initialize log-likelihood for this draw
                    // compute cumulative hazard up to the censoring time
                    real cum_hazard = 0;
                    for (time_step in 1:censor_time) {
                        real ind_term = 1.0;
                        real soc_term = s_prime * (sum(A_value[trial, time_step , d][id, ] .* Z[trial][time_step, ])) ;
                        real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                        cum_hazard += lambda; // accumulate hazard
                    }
                // Compute per-individual log likelihood
                log_likelihoods[d] += -cum_hazard;
           }
           log_lik_matrix[trial, N[trial] + c] = log_sum_exp(log_likelihoods) - log(S);
            }
        }
    }
    // Flatten log_lik_matrix into log_lik
    array[K * Q] real log_lik;
    int idx = 1;
    for (trial in 1:K) {
        for (n in 1:Q) {
            log_lik[idx] = log_lik_matrix[trial, n];
            idx += 1;
        }
    }
}
