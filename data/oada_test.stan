//stan
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
   real<lower=0> s = exp(log_s_mean);
}

model {
    log_s_mean ~ normal(0,2);
    
    for (trial in 1:K) {
        for (n in 1:N[trial]) {
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            real ind = 1;
            real soc = s * (sum(A_default[trial, learn_time][id, ] .* C[trial][learn_time, ]));
            real lambda = 1.0 * (ind + soc);
            target += log(lambda);
            
            vector[Q] rates = rep_vector(0.0, Q);
            
            for (j in 1:Q) {
                real temp_ind = 1;
                real temp_soc = s * (sum(A_default[trial, learn_time][j, ] .* C[trial][learn_time, ]));
                real temp_lambda = 1.0 * (temp_ind + temp_soc);
                rates[j] += temp_lambda * (1-C[trial][learn_time, j]);
            }
                        //print(rates[id]/sum(rates));
            target += -log(sum(rates));
        }
    }       
}

