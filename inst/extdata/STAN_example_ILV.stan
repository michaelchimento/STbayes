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
    array[K, T_max] matrix[P, P] A_kin, A_inverse_distance; // Network matrices
    array[K] matrix[T_max, P] Z;   // Knowledge state slash cue matrix
    array[P] real ILV_age;
array[K,Q,P] real ILV_dist_from_resource;
array[P] real ILV_sex;
    int<lower=0> N_veff;
}
parameters {
    real log_lambda_0_mean;  // Log baseline learning rate
    real log_s_mean; // Overall social transmission rate
    simplex[2] w; // Weights for networks
    real beta_ILVi_age;
real beta_ILVi_dist_from_resource;
    real beta_ILVs_sex;
}
transformed parameters {
   real<lower=0> lambda_0 = exp(log_lambda_0_mean);
real<lower=0> s_prime = exp(log_s_mean);
real<lower=0> s = s_prime/lambda_0;
}
model {
    log_lambda_0_mean ~ normal(-4, 3);
    log_s_mean ~ normal(-4, 3);
    w ~ dirichlet(rep_vector(0.5, 2));
    beta_ILVi_age ~ normal(0, 1);
beta_ILVi_dist_from_resource ~ normal(0, 1);
    beta_ILVs_sex ~ normal(0, 1);
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            if (learn_time > 0) {
                for (time_step in 1:learn_time) {
                    real ind_term = exp(beta_ILVi_age * ILV_age[id] + beta_ILVi_dist_from_resource * ILV_dist_from_resource[trial,time_step,id]);
                    real soc_term = s_prime * (w[1] * sum(A_kin[trial, time_step][id, ] .* Z[trial][time_step, ]) + w[2] * sum(A_inverse_distance[trial, time_step][id, ] .* Z[trial][time_step, ])) * exp(beta_ILVs_sex * ILV_sex[id]);
                    real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                    target += -lambda;
                    if (time_step == learn_time) {
                        target += log( (lambda_0 * ind_term + soc_term));
                    }
                }
            }
        }
        if (N_c[trial] > 0) {
            for (c in 1:N_c[trial]) {
                int id = ind_id[trial, N[trial] + c];
                    for (time_step in 1:T[trial]) {
                        real ind_term = exp(beta_ILVi_age * ILV_age[id] + beta_ILVi_dist_from_resource * ILV_dist_from_resource[trial,time_step,id]);
                        real soc_term = s_prime * (w[1] * sum(A_kin[trial, time_step][id, ] .* Z[trial][time_step, ]) + w[2] * sum(A_inverse_distance[trial, time_step][id, ] .* Z[trial][time_step, ])) * exp(beta_ILVs_sex * ILV_sex[id]);
                        real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                        target += -lambda;
                    }
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
            if (learn_time > 0){
                real cum_hazard = 0; //set val before adding
                for (time_step in 1:learn_time) {
                    real ind_term = exp(beta_ILVi_age * ILV_age[id] + beta_ILVi_dist_from_resource * ILV_dist_from_resource[trial,time_step,id]);
                    real soc_term = s_prime * (w[1] * sum(A_kin[trial, time_step][id, ] .* Z[trial][time_step, ]) + w[2] * sum(A_inverse_distance[trial, time_step][id, ] .* Z[trial][time_step, ])) * exp(beta_ILVs_sex * ILV_sex[id]);
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
                        real ind_term = exp(beta_ILVi_age * ILV_age[id] + beta_ILVi_dist_from_resource * ILV_dist_from_resource[trial,time_step,id]);
                        real soc_term = s_prime * (w[1] * sum(A_kin[trial, time_step][id, ] .* Z[trial][time_step, ]) + w[2] * sum(A_inverse_distance[trial, time_step][id, ] .* Z[trial][time_step, ])) * exp(beta_ILVs_sex * ILV_sex[id]);
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

