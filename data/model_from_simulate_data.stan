data {
    int<lower=0> K;                // Number of trials
    int<lower=0> Q;                // Number of individuals in each trial
    int<lower=1> Z;                // Number of unique individuals
    array[K] int<lower=1> N;       // Number of individuals that learned during observation period
    array[K] int<lower=0> N_c;     // Number of right-censored individuals
    array[K, Q] int<lower=1> ind_id; // IDs of individuals
    array[K] int<lower=1> T;                // Maximum time periods
    int<lower=1> T_max;             //max timesteps reached
    array[K,Z] int<lower=0> t;              // Time of acquisition for each individual
    array[K, T_max] real<lower=0> D;   // Scaled durations
    array[K, T_max] matrix[Z, Z] A_assoc; // Network matrices
    array[K] matrix[T_max, Z] C;        // knowledge state slash cue matrix
    
}

parameters {
    real log_lambda_0_mean;  // Log baseline learning rate
    real log_s_mean; // Overall social transmission rate
    
    
    
    

}

transformed parameters {
    real<lower=0> lambda_0 = 1 / exp(log_lambda_0_mean);
    real<lower=0> s = exp(log_s_mean);
}

model {
    log_lambda_0_mean ~ normal(6, 2);
    log_s_mean ~ normal(1, 2);
    
    
    
    
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];

            if (learn_time>0){
                for (time_step in 1:learn_time) {
                    real ind_term = 1;
                    real soc_term = s * (sum(A_assoc[trial, time_step][id, ] .* C[trial][time_step, ])) ;
                    real lambda =  lambda_0 * (ind_term + soc_term) * D[trial, time_step];
                    target += -lambda; //cumulative hazard (accumulating)
                    if (time_step == learn_time) {
                        target += log( lambda_0 * (ind_term + soc_term)); //inst. hazard
                    }
                }
            }
        }

        // Contributions of censored individuals
        if (N_c[trial] > 0) {
            for (c in 1:N_c[trial]) {
                int id = ind_id[trial, N[trial] + c];

                for (time_step in 1:T[trial]) {
                    real ind_term = 1;
                    real soc_term = s * (sum(A_assoc[trial, time_step][id, ] .* C[trial][time_step, ])) ;
                    real lambda =  lambda_0 * (ind_term + soc_term) * D[trial, time_step];
                    target += -lambda;
                }
            }
        }
    }
}
