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
    array[K,P] int t;     // Time of acquisition for each individual
    array[K, T_max] real<lower=0> D; // Scaled durations
    int<lower=1> N_networks;
    array[N_networks, K, T_max] matrix[P, P] A;  // network matrices
    array[K] matrix[T_max, P] Z;   // Knowledge state * cue matrix
    array[K] matrix[T_max, P] Zn;   // Knowledge state
    matrix[P,1] ILV_bool_ILV;
    vector[P] ILV_cont_ILV;
    matrix[P,3] ILV_cat_ILV;
    int<lower=0> N_veff;
}
parameters {
    real log_lambda_0_mean;  // Log baseline learning rate
    real log_s_prime_mean;
    vector[1] beta_ILVi_bool_ILV; //boolean ILV (additive, intrinsic)
    real beta_ILVs_cont_ILV; //continuous ILV (additive, social)
    vector[3] beta_ILVm_cat_ILV; //categorical ILV (multiplicative)
}
transformed parameters {
    real s_prime;
    real<lower=0> lambda_0;
    vector[P] bool_ILV_i;
    vector[P] cont_ILV_s;
    vector[P] cat_ILV_m;
    s_prime = exp(log_s_prime_mean);
    lambda_0 = exp(log_lambda_0_mean);
    bool_ILV_i = ILV_bool_ILV * beta_ILVi_bool_ILV;
    cont_ILV_s = ILV_cont_ILV * beta_ILVs_cont_ILV;
    cat_ILV_m = ILV_cat_ILV * beta_ILVm_cat_ILV;
}
model {
    log_lambda_0_mean ~ normal(-4, 2);
    log_s_prime_mean ~ normal(-4, 2);
    beta_ILVi_bool_ILV ~ normal(0,1);
    beta_ILVs_cont_ILV ~ normal(0,1);
    beta_ILVm_cat_ILV ~ normal(0,1);
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            if (learn_time > 0) {
                for (time_step in 1:learn_time) {
                    real ind_term = exp(bool_ILV_i[id]);
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        net_effect += s_prime * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                    }
                    real soc_term = net_effect* exp(cont_ILV_s[id]);
                    real lambda = exp(cat_ILV_m[id]) * (lambda_0 * ind_term + soc_term) * D[trial, time_step] ;
                    target += -lambda;
                    if (time_step == learn_time) {
                        target += log(exp(cat_ILV_m[id]) * (lambda_0 * ind_term + soc_term));
                    }
                }
            }
        }
        if (N_c[trial] > 0) {
            for (c in 1:N_c[trial]) {
                int id = ind_id[trial, N[trial] + c];
                for (time_step in 1:T[trial]) {
                    real ind_term = exp(bool_ILV_i[id]);
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        net_effect += s_prime * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                    }
                    real soc_term = net_effect* exp(cont_ILV_s[id]);
                    real lambda = exp(cat_ILV_m[id]) * (lambda_0 * ind_term + soc_term) * D[trial, time_step] ;
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
                    real ind_term = exp(bool_ILV_i[id]);
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        net_effect += s_prime * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                    }
                    real soc_term = net_effect* exp(cont_ILV_s[id]);
                    real lambda = exp(cat_ILV_m[id]) * (lambda_0 * ind_term + soc_term) * D[trial, time_step] ;
                    cum_hazard += lambda; // accumulate hazard
                    //if it learn_time, record the ll
                    if (time_step == learn_time){
                        log_lik_matrix[trial, n] = log(exp(cat_ILV_m[id]) * (lambda_0 * ind_term + soc_term)) - cum_hazard;
                        for (network in 1:N_networks) {
                            real Tn = dot_product(A[network, trial, time_step][id, ], Z[trial][time_step, ]);
                            psocn_sum[network] += (s_prime * D[trial, time_step] * exp(cont_ILV_s[id])  * Tn) / lambda;
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
                    real ind_term = exp(bool_ILV_i[id]);
                    real net_effect = 0;
                    for (network in 1:N_networks) {
                        net_effect += s_prime * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
                    }
                    real soc_term = net_effect* exp(cont_ILV_s[id]);
                    real lambda = exp(cat_ILV_m[id]) * (lambda_0 * ind_term + soc_term) * D[trial, time_step] ;
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
