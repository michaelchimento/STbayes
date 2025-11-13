//stan
functions {
    real dini_func(real x, real k) {
        // transform x from [0,1] to [-1,1]
        real x_transformed = 2 * x - 1;
        real y = ((x_transformed - k * x_transformed) / (k - 2 * k * abs(x_transformed) + 1) + 1) / 2;
        return y;
    }
}
data {
    int<lower=0> K;                // Number of trials
    int<lower=0> Q;                // Number of individuals in each trial
    int<lower=1> P;                // Number of unique individuals
    array[K] int<lower=0> N;       // Number of individuals that learned during observation period
    array[K] int<lower=0> N_c;     // Number of right-censored individuals
    array[K, Q] int<lower=-1> ind_id; // IDs of individuals
    array[K] int<lower=1> T;       // Maximum time periods
    int<lower=1> T_max;            // Max timesteps reached
    array[K,P] int t;     // Time of acquisition for each individual
    array[K, T_max] real<lower=0> D; // Scaled durations
    int<lower=1> N_networks;
    array[N_networks, K, T_max] matrix[P, P] A;  // network matrices
    array[K] matrix[T_max, P] Z;   // Knowledge state * cue matrix
    array[K] matrix[T_max, P] Zn;   // Knowledge state
    int<lower=0> N_veff;
}
parameters {
    real log_lambda_0_mean;  // Log baseline learning rate
    real log_s_prime_mean;
    real k_raw;
}
transformed parameters {
    real s_prime;
    real<lower=0> lambda_0;
    real<lower=-1, upper=1> k_shape;
    s_prime = exp(log_s_prime_mean);
    lambda_0 = exp(log_lambda_0_mean);
    k_shape = 2 / (1 + exp(-k_raw)) - 1;
}
model {
    log_lambda_0_mean ~ normal(-4, 2);
    log_s_prime_mean ~ normal(-4, 2);
    k_raw ~ normal(0,3);
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            if (learn_time > 0) {
                for (time_step in 1:learn_time) {
                    real ind_term = 1.0;
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                        real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
                        real prop = denom > 0 ? numer / denom : 0.0;
                        real dini_transformed = dini_func(prop, k_shape);
                        net_effect += s_prime * dini_transformed;
                    }
                    real soc_term = net_effect;
                    real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step] ;
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
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                        real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
                        real prop = denom > 0 ? numer / denom : 0.0;
                        real dini_transformed = dini_func(prop, k_shape);
                        net_effect += s_prime * dini_transformed;
                    }
                    real soc_term = net_effect;
                    real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step] ;
                    target += -lambda;
                }
            }
        }
    }
}
generated quantities {
    real<lower=0> s = s_prime/lambda_0;
    matrix[K, Q] log_lik_matrix = rep_matrix(0.0, K, Q);           // LL for each observation
    //for %ST
    int count_ST = 0;
    vector[N_networks] psocn_sum = rep_vector(0.0, N_networks);
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            if (learn_time > 0){
                real cum_hazard = 0; //set val before adding
                for (time_step in 1:learn_time) {
                    real ind_term = 1.0;
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                        real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
                        real prop = denom > 0 ? numer / denom : 0.0;
                        real dini_transformed = dini_func(prop, k_shape);
                        net_effect += s_prime * dini_transformed;
                    }
                    real soc_term = net_effect;
                    real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step] ;
                    cum_hazard += lambda; // accumulate hazard
                    //if it learn_time, record the ll
                    if (time_step == learn_time){
                        log_lik_matrix[trial, n] = log( (lambda_0 * ind_term + soc_term)) - cum_hazard;
                        for (network in 1:N_networks) {
                            real numer = dot_product(A[network, trial, time_step][id, ], Z[trial][time_step, ]);
                            real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
                            real prop = denom > 0 ? numer / denom : 0.0;
                            real dini = dini_func(prop, k_shape);
                            psocn_sum[network] += (s_prime * D[trial, time_step]   * dini) / lambda;
                        }
                        count_ST += 1;
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
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                        real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
                        real prop = denom > 0 ? numer / denom : 0.0;
                        real dini_transformed = dini_func(prop, k_shape);
                        net_effect += s_prime * dini_transformed;
                    }
                    real soc_term = net_effect;
                    real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step] ;
                    cum_hazard += lambda; // accumulate hazard
                }
                // Compute per-individual log likelihood
                log_lik_matrix[trial, N[trial] + c] = -cum_hazard;
            }
        }
    }
    vector[N_networks] percent_ST = psocn_sum / count_ST;
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
