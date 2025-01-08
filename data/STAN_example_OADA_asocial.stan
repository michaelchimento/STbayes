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
    array[K, T_max] matrix[Z, Z] A_default; // Network matrices
    array[K] matrix[T_max, Z] C;   // Knowledge state slash cue matrix
    int<lower=0> N_veff;
}
parameters {
    real log_s_mean;         // Overall social transmission rate
}
transformed parameters {
   real<lower=0> s = 0.0;
}
model {
    log_s_mean ~ normal(0,2);

    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            int time_step = learn_time;
            if (learn_time > 0) {
                real i_ind = 1.0;
                real i_soc = s * (sum(A_default[trial, time_step][id, ] .* C[trial][time_step, ])) ;
                real i_lambda = 1.0 * (i_ind + i_soc);

                vector[Q] j_rates = rep_vector(0.0, Q);

                for (j in 1:Q) {
                    real j_ind = 1.0;
                    real j_soc = s * (sum(A_default[trial, time_step][j, ] .* C[trial][time_step, ])) ;
                    real j_lambda = 1.0 * (j_ind + j_soc);
                    j_rates[j] += j_lambda * (1-C[trial][learn_time, j]); //only include those who haven't learned in denom
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
                    real i_soc = s * (sum(A_default[trial, time_step][id, ] .* C[trial][time_step, ])) ;
                    real i_lambda = 1.0 * (i_ind + i_soc);

                    vector[Q] j_rates = rep_vector(0.0, Q);

                    for (j in 1:Q) {
                        real j_ind = 1.0;
                        real j_soc = s * (sum(A_default[trial, time_step][j, ] .* C[trial][time_step, ])) ;
                        real j_lambda = 1.0 * (j_ind + j_soc);
                        j_rates[j] += j_lambda * (1-C[trial][learn_time, j]);
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
