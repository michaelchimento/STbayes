#' Dynamically generate STAN model based on input data
#'
#' @param STb_data a list of formatted data returned from the STbayes_data() function
#' @param N_veff Number of varying effects. Defaults to 2 (one for strength of SL, one for asocial rate)
#' @param gq Boolean to indicate whether the generated quantities block is added (incl. ll for WAIC)
#' @param est_acqTime Boolean to indicate whether gq block includes estimates for acquisition time. At the moment this uses 'one weird trick' to accomplish this and does not support estimates for non-integer learning times.
#'
#' @return A STAN model (character) that is customized to the input data.
#' @export
#'
#' @examples
#' #very mock data
#' diffusion_data <- data.frame(
#'   id = c("A", "B", "C", "D", "E", "F"),
#'   trial = c(1, 1, 1, 2, 2, 2),
#'   time = c(0, 1, 2, 0, 1, 4),
#'   max_time = c(3, 3, 3, 4, 4, 4)
#' )
#' networks <- data.frame(
#'   trial = c(1, 1, 1, 2, 2, 2),
#'   from = c("A", "A", "B", "D", "D", "E"),
#'   to = c("B", "C", "C", "E", "F", "F"),
#'   kin = c(1, 0, 1, 0, 1, 1),
#'   inverse_distance = c(0, 1, .5, .25, .1, 0)
#' )
#' ILV_metadata <- data.frame(
#'   id = c("A", "B", "C", "D", "E", "F"),
#'   age = c(2, 3, 4, 2, 5, 6),
#'   sex = c(0, 1, 1, 0, 1, 0) # Factor ILVs must be input as numeric
#' )
#' data_list <- import_user_STb(
#'   diffusion_data = diffusion_data,
#'   networks = networks,
#'   ILV_metadata = ILV_metadata,
#'   ILVi = c("age"), # Use only 'age' for asocial learning
#'   ILVs = c("sex") # Use only 'sex' for social learning
#' )
#'
#' model = generate_STb_model(data_list)
#' model = generate_STb_model(data_list, N_veff=2) # estimate varying effects for s and lambda_0.
#' print(model)
generate_STb_model <- function(STb_data, N_veff = 0, gq = TRUE, est_acqTime = FALSE) {
    network_names = STb_data$network_names
    ILVi_vars = STb_data$ILVi_names[!STb_data$ILVi_names %in% "ILVabsent"]
    ILVs_vars = STb_data$ILVs_names[!STb_data$ILVs_names %in% "ILVabsent"]
    ILVm_vars = STb_data$ILVm_names[!STb_data$ILVm_names %in% "ILVabsent"]

    combined_ILV_vars = unique(c(ILVi_vars, ILVs_vars, ILVm_vars))

    ILV_declaration = ""
    if (length(combined_ILV_vars)>0){
        ILV_declaration = paste0('array[Z] real ', combined_ILV_vars, ';', collapse = '\\n')
    }

    # Placeholders for dynamic components
    num_networks <- length(STb_data$network_names)
    num_ILVi = length(ILVi_vars)
    num_ILVs = length(ILVs_vars)
    num_ILVm = length(ILVm_vars)

    if (num_networks == 1) {
        network_term <- paste0("sum(A_", network_names[1], "[trial, time_step][id, ] .* C[trial][time_step, ])")
        w_param <- ""
        w_prior <- ""
    } else {
        network_term <- paste0("w[", 1:num_networks, "] * sum(A_", network_names, "[trial, time_step][id, ] .* C[trial][time_step, ])", collapse = " + ")
        w_param <- paste0("simplex[", num_networks, "] w; // Weights for networks")
        w_prior <- paste0("w ~ dirichlet(rep_vector(0.5, ", num_networks, "));")
    }

    # Handle individual-level information (ILVi)
    if (num_ILVi < 1) {
        ILVi_param <- ""
        ILVi_prior <- ""
        ILVi_variable_effects <- "1"
    } else {
        ILVi_param <- paste0("real beta_ILVi_", ILVi_vars, ";", sep = "\n")
        ILVi_prior <- paste0("beta_ILVi_", ILVi_vars, " ~ normal(0, 1);", sep = "\n")
        ILVi_variable_effects <- paste0("exp(", paste0("beta_ILVi_", ILVi_vars, " * ", ILVi_vars, "[id]", collapse = " + "),")")
    }

    # Handle social-level information (ILVs)
    if (num_ILVs < 1) {
        ILVs_param <- ""
        ILVs_prior <- ""
        ILVs_variable_effects <- ""
    } else {
        ILVs_param <- paste0("real beta_ILVs_", ILVs_vars, ";", sep = "\n")
        ILVs_prior <- paste0("beta_ILVs_", ILVs_vars, " ~ normal(0, 1);", sep = "\n")
        ILVs_variable_effects <- paste0("* exp(", paste0("beta_ILVs_", ILVs_vars, " * ", ILVs_vars, "[id]", collapse = " + "),")")
    }

    # Handle social-level information (ILVs)
    if (num_ILVm < 1) {
        ILVm_param <- ""
        ILVm_prior <- ""
        ILVm_variable_effects <- ""
    } else {
        ILVm_param <- paste0("real beta_ILVm_", ILVm_vars, ";", sep = "\n")
        ILVm_prior <- paste0("beta_ILVm_", ILVm_vars, " ~ normal(0, 1);", sep = "\n")
        ILVm_variable_effects <- paste0("exp(", paste0("beta_ILVm_", ILVm_vars, " * ", ILVm_vars, "[id] *", collapse = " + "),")")
    }

    data_block <- glue::glue("
data {{
    int<lower=0> K;                // Number of trials
    int<lower=0> Q;                // Number of individuals in each trial
    int<lower=1> Z;                // Number of unique individuals
    array[K] int<lower=1> N;       // Number of individuals that learned during observation period
    array[K] int<lower=0> N_c;     // Number of right-censored individuals
    array[K, Q] int<lower=1> ind_id; // IDs of individuals
    array[K] int<lower=1> T;       // Maximum time periods
    int<lower=1> T_max;            // Max timesteps reached
    {if (est_acqTime) 'array[K] int<lower=0> time_max; //Duration of obs period for each trial' else ''}
    array[K,Z] int<lower=0> t;     // Time of acquisition for each individual
    array[K, T_max] real<lower=0> D; // Scaled durations
    array[K, T_max] matrix[Z, Z] {paste0('A_', network_names, collapse = ', ')}; // Network matrices
    array[K] matrix[T_max, Z] C;   // Knowledge state slash cue matrix
    {ILV_declaration}
    int<lower=0> N_veff;
    array[K, T_max] int<lower=0> D_int; // integer durations
}}
")

    # Parameters block
    parameters_block <- glue::glue("
parameters {{
    real log_lambda_0_mean;  // Log baseline learning rate
    real log_s_mean;         // Overall social transmission rate
    {w_param}
    {ILVi_param}
    {ILVs_param}
    {ILVm_param}
    {if (N_veff == 2) 'matrix[N_veff,Z] z_ID; vector<lower=0>[N_veff] sigma_ID; cholesky_factor_corr[N_veff] Rho_ID;' else ''}
}}
")

    # Transformed parameters block
    transformed_parameters_block <- glue::glue("
transformed parameters {{
    {if (N_veff == 2) '
        matrix[Z,N_veff] v_ID;
        v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)\\';
        vector<lower=0>[Z] lambda_0 = 1 / exp(log_lambda_0_mean + v_ID[,1]);
        vector<lower=0>[Z] s = exp(log_s_mean + v_ID[,2]);
    ' else '
        real<lower=0> lambda_0 = 1 / exp(log_lambda_0_mean);
        real<lower=0> s = exp(log_s_mean);
    '}
}}
")

    # Model block
    model_block <- glue::glue("
model {{
    log_lambda_0_mean ~ normal(6, 2);
    log_s_mean ~ normal(1, 2);
    {w_prior}
    {ILVi_prior}
    {ILVs_prior}
    {ILVm_prior}

    for (trial in 1:K) {{
        for (n in 1:N[trial]) {{
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];

            if (learn_time > 0) {{
                for (time_step in 1:learn_time) {{
                    real ind_term = {ILVi_variable_effects};
                    real soc_term = {if (N_veff == 2) 's[id]' else 's'} * ({network_term}) {ILVs_variable_effects};
                    real lambda = {ILVm_variable_effects} {if (N_veff == 2) 'lambda_0[id]' else 'lambda_0'} * (ind_term + soc_term) * D[trial, time_step];
                    target += -lambda;
                    if (time_step == learn_time) {{
                        target += log({ILVm_variable_effects} {if (N_veff == 2) 'lambda_0[id]' else 'lambda_0'} * (ind_term + soc_term));
                    }}
                }}
            }}
        }}

        if (N_c[trial] > 0) {{
            for (c in 1:N_c[trial]) {{
                int id = ind_id[trial, N[trial] + c];

                for (time_step in 1:T[trial]) {{
                    real ind_term = {ILVi_variable_effects};
                    real soc_term = {if (N_veff == 2) 's[id]' else 's'} * ({network_term}) {ILVs_variable_effects};
                    real lambda = {ILVm_variable_effects} {if (N_veff == 2) 'lambda_0[id]' else 'lambda_0'} * (ind_term + soc_term) * D[trial, time_step];
                    target += -lambda;
                }}
            }}
        }}
    }}
}}
")

est_acqTime_code <- if (est_acqTime==TRUE) glue::glue("
    matrix[K, Q] acquisition_time;         // simulated acquisition times
    for (trial in 1:K) {{
        for (n in 1:Q) {{ //have to loop through bc stan
            acquisition_time[trial, n] = time_max[trial];
        }}
        for (n in 1:Q) {{
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            if (learn_time > 0){{
                real cum_hazard = 0; //set val before adding
                int global_time = 1;
                for (time_step in 1:T[trial]) {{
                    for (micro_time in 1:D_int[trial, time_step]){{
                        real ind_term = {ILVi_variable_effects};
                        real soc_term = {if (N_veff == 2) 's[id]' else 's'} * ({network_term}) {ILVs_variable_effects};
                        real lambda = {ILVm_variable_effects} {if (N_veff == 2) 'lambda_0[id]' else 'lambda_0'} * (ind_term + soc_term);
                        real prob = 1-exp(-lambda);
                        if (bernoulli_rng(prob) && acquisition_time[trial, n]>=time_max[trial]) {{
                            acquisition_time[trial, n] = global_time;
                        }}
                        global_time += 1;
                    }}
                }}
            }}
        }}
    }}"
) else ""

    generated_quantities_block <- glue::glue("
generated quantities {{
    matrix[K, Q] log_lik_matrix = rep_matrix(0.0, K, Q);           // LL for each observation

    for (trial in 1:K) {{
        for (n in 1:N[trial]) {{
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];

            if (learn_time > 0){{
                real cum_hazard = 0; //set val before adding
                for (time_step in 1:T[trial]) {{
                    real ind_term = {ILVi_variable_effects};
                    real soc_term = {if (N_veff == 2) 's[id]' else 's'} * ({network_term}) {ILVs_variable_effects};
                    real lambda = {ILVm_variable_effects} {if (N_veff == 2) 'lambda_0[id]' else 'lambda_0'} * (ind_term + soc_term) * D[trial, time_step];
                    cum_hazard += lambda; // accumulate hazard
                    //if it learn_time, record the ll
                    if (time_step == learn_time){{
                        log_lik_matrix[trial, n] = log({ILVm_variable_effects} {if (N_veff == 2) 'lambda_0[id]' else 'lambda_0'} * (ind_term + soc_term)) - cum_hazard;
                    }}
                }}
            }}
        }}

        // Contributions of censored individuals
        if (N_c[trial] > 0) {{
            for (c in 1:N_c[trial]) {{
                int id = ind_id[trial, N[trial] + c];
                int censor_time = T[trial]; // Censoring time (end of observation)

                // compute cumulative hazard up to the censoring time
                real cum_hazard = 0;
                for (time_step in 1:censor_time) {{
                    real ind_term = {ILVi_variable_effects};
                    real soc_term = {if (N_veff == 2) 's[id]' else 's'} * ({network_term}) {ILVs_variable_effects};
                    real lambda = {ILVm_variable_effects} {if (N_veff == 2) 'lambda_0[id]' else 'lambda_0'} * (ind_term + soc_term) * D[trial, time_step];
                    cum_hazard += lambda; // accumulate hazard
                }}
                // Compute per-individual log likelihood
                log_lik_matrix[trial, N[trial] + c] = -cum_hazard;
            }}
        }}
    }}

    {est_acqTime_code}

    // Flatten log_lik_matrix into log_lik
    array[K * Q] real log_lik;
    int idx = 1;
    for (trial in 1:K) {{
        for (n in 1:Q) {{
            log_lik[idx] = log_lik_matrix[trial, n];
            idx += 1;
        }}
    }}
}}
")

    # combine all blocks
    stan_model <- glue::glue("{data_block}
                             {parameters_block}
                             {transformed_parameters_block}
                             {model_block}
                             {if (gq==T) {generated_quantities_block} else ''}")



    return(stan_model)
}
