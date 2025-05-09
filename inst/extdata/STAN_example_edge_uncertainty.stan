//stan
data {
    int<lower=0> K;                // Number of trials
    int<lower=0> Q;                // Number of individuals in each trial
    int<lower=1> P;                // Number of unique individuals
    array[K] int<lower=0> N;       // Number of individuals that learned during observation period
    array[K] int<lower=0> N_c;     // Number of right-censored individuals
    array[K, Q] int<lower=-1> ind_id; // IDs of individuals
    array[K] int<lower=1> T;       // Maximum time periods
    int<lower=1> T_max;            // Max timesteps reached
    array[K,P] int<lower=-1> t;     // Time of acquisition for each individual
    array[K, T_max] real<lower=0> D; // Scaled durations
    array[K] matrix[T_max, P] Z;   // Knowledge state slash cue matrix
    int<lower=0> N_veff;
    int N_dyad;  // number of dyads
    vector[N_dyad] logit_edge_mu_A_1;  // logit edge values for A_1
    matrix[N_dyad, N_dyad] logit_edge_cov_A_1;  // covariance matrix for A_1
}
parameters {
    vector[N_dyad] edge_logit_A_1; //edge weights
    real log_lambda_0_mean;  // Log baseline learning rate
    real log_s_mean; // Overall social transmission rate
}
transformed parameters {
    real<lower=0> lambda_0 = exp(log_lambda_0_mean);
    real<lower=0> s_prime = exp(log_s_mean);
    //convert to matrix for dot product l8r
    matrix[P, P] A_1;
    {
      int edge_idx = 1;
      for (i in 1:P) {
        for (j in 1:P) {
          if (i != j) {
            A_1[i, j] = inv_logit(edge_logit_A_1[edge_idx]);
            edge_idx += 1;
          } else {
            A_1[i, j] = 0;
          }
        }
      }
    }
}
model {
    //priors
    log_lambda_0_mean ~ normal(-4, 3);
    log_s_mean ~ normal(-4, 3);
    edge_logit_A_1 ~ multi_normal(logit_edge_mu_A_1, logit_edge_cov_A_1);

    //for each diffusion trial
    for (trial in 1:K) {
        //for each non-censored ind
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            //if not tutor/demo
            if (learn_time > 0) {
                //likelihood for each inter-event interval
                for (time_step in 1:learn_time) {
                    real ind_term = 1.0;
                    real soc_term = s_prime * (sum(A_1[id, ] .* Z[trial][time_step, ]));
                    real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                    target += -lambda;
                    if (time_step == learn_time) {
                        target += log( (lambda_0 * ind_term + soc_term));
                    }
                }
            }
        }
        //for each censored individual
        if (N_c[trial] > 0) {
            for (c in 1:N_c[trial]) {
                int id = ind_id[trial, N[trial] + c];
                    for (time_step in 1:T[trial]) {
                        real ind_term = 1.0;
                        real soc_term = s_prime * (sum(A_1[id, ] .* Z[trial][time_step, ])) ;
                        real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                        target += -lambda;
                    }
            }
        }
    }
}
generated quantities {
    real<lower=0> s = s_prime/lambda_0;
    matrix[K, Q] log_lik_matrix = rep_matrix(0.0, K, Q);           // LL for each observation
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            if (learn_time > 0){
                real cum_hazard = 0; //set val before adding
                for (time_step in 1:learn_time) {
                    real ind_term = 1.0;
                    real soc_term = s_prime * (sum(A_1[id, ] .* Z[trial][time_step, ])) ;
                    real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                    cum_hazard += lambda; // accumulate hazard
                    //if it learn_time, record the ll
                    if (time_step == learn_time){
                        log_lik_matrix[trial, n] = log( (lambda_0 * ind_term + soc_term)) - cum_hazard;
                    }
                }
            }
        }
        // Contributions of censored individuals
        if (N_c[trial] > 0) {
            for (c in 1:N_c[trial]) {
                int id = ind_id[trial, N[trial] + c];
                int censor_time = T[trial]; // Censoring time (end of observation)
                    // compute cumulative hazard up to the censoring time
                    real cum_hazard = 0;
                    for (time_step in 1:censor_time) {
                        real ind_term = 1.0;
                        real soc_term = s_prime * (sum(A_1[id, ] .* Z[trial][time_step, ])) ;
                        real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                        cum_hazard += lambda; // accumulate hazard
                    }
                // Compute per-individual log likelihood
                log_lik_matrix[trial, N[trial] + c] = -cum_hazard;
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
