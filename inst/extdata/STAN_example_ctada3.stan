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
    array[K] int<lower=0> time_max; //Duration of obs period for each trial
    array[K, T_max] int<lower=0> D_int; // integer durations
}
parameters {
    real<lower=0, upper=1> log_lambda_0_mean;  // Log baseline learning rate
    real<lower=0, upper=1> scaled_s;  // Scaled social term
}
transformed parameters {
   real<lower=0> lambda_0 = log_lambda_0_mean;
   real<lower=0> s_prime = scaled_s;
   real<lower=0> s = s_prime/lambda_0;
}
model {
    log_lambda_0_mean ~ uniform(0, 1);
    scaled_s ~ uniform(0, 1);
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            if (learn_time > 0) {
                for (time_step in 1:learn_time) {
                    real ind_term = 1.0;
                    real soc_term = s_prime * (sum(A_assoc[trial, time_step][id, ] .* C[trial][time_step, ])) ;
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
                    real soc_term = s_prime * (sum(A_assoc[trial, time_step][id, ] .* C[trial][time_step, ])) ;
                    real lambda =  (lambda_0 * ind_term + soc_term) * D[trial, time_step];
                    target += -lambda;
                }
            }
        }
    }
}
