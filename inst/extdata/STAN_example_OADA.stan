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
    int<lower=0> N_veff;
}
parameters {
    real log_s_mean; // Overall social transmission rate
}
transformed parameters {
    real<lower=0> s_prime = exp(log_s_mean);
}
model {
    log_s_mean ~ normal(1,3);
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            int time_step = learn_time;
            if (learn_time > 0) {
                real i_ind = 1.0;
                real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += 1.0 * sum(A[network, trial, time_step][id, ] .* Z[trial][time_step, ]);
}
                real i_soc = s_prime * (net_effect);
                real i_lambda = 1.0 * (i_ind + i_soc);
                vector[Q] j_rates = rep_vector(0.0, Q);
                for (j in 1:Q) {
                    real j_ind = 1.0;
                    real net_effect_j = 0;
for (network in 1:N_networks) {
  net_effect_j += 1.0 * sum(A[network, trial, time_step][j, ] .* Z[trial][time_step, ]);
}
                    real j_soc = s_prime * (net_effect_j);
                    real j_lambda = 1.0 * (j_ind + j_soc);
                    j_rates[j] += j_lambda * (1-Z[trial][learn_time, j]); //only include those who haven't learned in denom
                }
                target += log(i_lambda) - log(sum(j_rates));
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
                int time_step = learn_time;
                if (learn_time > 0) {
                    real i_ind = 1.0;
                    real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += 1.0 * sum(A[network, trial, time_step][id, ] .* Z[trial][time_step, ]);
}
                    real i_soc = s_prime * (net_effect);
                    real i_lambda = 1.0 * (i_ind + i_soc);
                    vector[Q] j_rates = rep_vector(0.0, Q);
                    for (j in 1:Q) {
                        real j_ind = 1.0;
                        real net_effect_j = 0;
for (network in 1:N_networks) {
  net_effect_j += 1.0 * sum(A[network, trial, time_step][j, ] .* Z[trial][time_step, ]);
}
                        real j_soc = s_prime * (net_effect_j);
                        real j_lambda = 1.0 * (j_ind + j_soc);
                        j_rates[j] += j_lambda * (1-Z[trial][learn_time, j]);
                    }
                    log_lik_matrix[trial, n] = log(i_lambda) - log(sum(j_rates));
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

