functions {
  real dini_func(real x, real k) {
    // transform x from [0,1] to [-1,1]
    real x_transformed = 2 * x - 1;
    real y = ((x_transformed - k * x_transformed) / (k - 2 * k * fabs(x_transformed) + 1) + 1) / 2;
    return y;
  }
}
data {
    int<lower=0> K;                // Number of trials
    int<lower=0> Q;                // Number of individuals in each trial
    int<lower=1> Z;                // Number of unique individuals
    array[K] int<lower=0> N;       // Number of individuals that learned during observation period
    array[K] int<lower=0> N_c;     // Number of right-censored individuals
    array[K, Q] int<lower=-1> ind_id; // IDs of individuals
    array[K] int<lower=1> T;       // Maximum time periods
    int<lower=1> T_max;            // Max timesteps reached
    array[K,Z] int<lower=-1> t;     // Time of acquisition for each individual
    array[K, T_max] real<lower=0> D; // Scaled durations
    array[K, T_max] matrix[Z, Z] A_assoc; // Network matrices
    array[K] matrix[T_max, Z] C;   // Knowledge state slash cue matrix
    int<lower=0> N_veff;
}
parameters {
    real log_lambda_0_mean;  // Log baseline learning rate
    real log_s_mean; // Overall social transmission rate
                                   real k_raw;
}
transformed parameters {
   real<lower=0> lambda_0 = 1 / exp(log_lambda_0_mean);
real<lower=0> s_prime = 1 / exp(log_s_mean);
real<lower=0> s = s_prime/lambda_0;
real<lower=-1, upper=1> k_shape = 2 / (1 + exp(-k_raw)) - 1;
}
model {
    log_lambda_0_mean ~ uniform(-10, 10);
    log_s_mean ~ uniform(-10, 10);
                              k_raw ~ normal(0,1);
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            if (learn_time > 0) {
                for (time_step in 1:learn_time) {
                    real ind_term = 1.0;
                    real sum_ass = sum(A_assoc[trial, time_step][id, ]);
                                       real prop_know = 0.0;
                                       if (sum_ass>0) prop_know = sum(A_assoc[trial, time_step][id, ] .* C[trial][time_step, ])/sum_ass;
                                       real dini_transformed = dini_func(prop_know,k_shape);
real soc_term = s_prime * (dini_transformed) ;
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
                    real ind_term = 1.0;
                    real sum_ass = sum(A_assoc[trial, time_step][id, ]);
                                       real prop_know = 0.0;
                                       if (sum_ass>0) prop_know = sum(A_assoc[trial, time_step][id, ] .* C[trial][time_step, ])/sum_ass;
                                       real dini_transformed = dini_func(prop_know,k_shape);
real soc_term = s_prime * (dini_transformed) ;
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
                for (time_step in 1:T[trial]) {
                    real ind_term = 1.0;
                    real sum_ass = sum(A_assoc[trial, time_step][id, ]);
                                       real prop_know = 0.0;
                                       if (sum_ass>0) prop_know = sum(A_assoc[trial, time_step][id, ] .* C[trial][time_step, ])/sum_ass;
                                       real dini_transformed = dini_func(prop_know,k_shape);
real soc_term = s_prime * (dini_transformed) ;
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
                    real sum_ass = sum(A_assoc[trial, time_step][id, ]);
                                       real prop_know = 0.0;
                                       if (sum_ass>0) prop_know = sum(A_assoc[trial, time_step][id, ] .* C[trial][time_step, ])/sum_ass;
                                       real dini_transformed = dini_func(prop_know,k_shape);
real soc_term = s_prime * (dini_transformed) ;
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
