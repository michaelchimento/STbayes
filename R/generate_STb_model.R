#' Dynamically generate STAN model based on input data
#'
#' @param STb_data a list of formatted data returned from the STbayes_data() function
#' @param N_veff Number of varying effects. Defaults to 2 (one for strength of SL, one for asocial rate)
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
generate_STb_model <- function(STb_data, N_veff = 0) {
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
        ILVi_variable_effects <- paste0("exp(", paste0("beta_ILVs_", ILVi_vars, " * ", ILVi_vars, "[id])", collapse = " + "),")")
    }

    # Handle social-level information (ILVs)
    if (num_ILVs < 1) {
        ILVs_param <- ""
        ILVs_prior <- ""
        ILVs_variable_effects <- ""
    } else {
        ILVs_param <- paste0("real beta_ILVs_", ILVs_vars, ";", sep = "\n")
        ILVs_prior <- paste0("beta_ILVs_", ILVs_vars, " ~ normal(0, 1);", sep = "\n")
        ILVs_variable_effects <- paste0("* exp(", paste0("beta_ILVs_", ILVs_vars, " * ", ILVs_vars, "[id])", collapse = " + "),")")
    }

    # Handle social-level information (ILVs)
    if (num_ILVm < 1) {
        ILVm_param <- ""
        ILVm_prior <- ""
        ILVm_variable_effects <- ""
    } else {
        ILVm_param <- paste0("real beta_ILVm_", ILVm_vars, ";", sep = "\n")
        ILVm_prior <- paste0("beta_ILVm_", ILVm_vars, " ~ normal(0, 1);", sep = "\n")
        ILVm_variable_effects <- paste0("exp(", paste0("beta_ILVm_", ILVm_vars, " * ", ILVm_vars, "[id]) *", collapse = " + "),")")
    }

    stan_model = ""

    if (N_veff==2){
        # Build the Stan model
        stan_model <- glue::glue("
        data {{
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
            array[K, T_max] matrix[Z, Z] {paste0('A_', network_names, collapse = ', ')}; // Network matrices
            array[K] matrix[T_max, Z] C;        // knowledge state slash cue matrix
            {ILV_declaration}
            int<lower=0> N_veff;
        }}

        parameters {{
            real log_lambda_0_mean;  // Log baseline learning rate
            real log_s_mean; // Overall social transmission rate
            {w_param}
            {ILVi_param}
            {ILVs_param}
            {ILVm_param}

            // ID effects
            matrix[N_veff,Z] z_ID;               // Matrix of uncorrelated z-values
            vector<lower=0>[N_veff] sigma_ID;    // SD of parameters among individuals
            cholesky_factor_corr[N_veff] Rho_ID; // Cholesky factor for correlation matrix

        }}

        transformed parameters {{
            matrix[Z,N_veff] v_ID; // Matrix of varying effects for each individual
            v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)';

            vector<lower=0>[Z] lambda_0 = 1 / exp(log_lambda_0_mean + v_ID[,1]);
            vector<lower=0>[Z] s = exp(log_s_mean + v_ID[,2]);
        }}

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

                    if (learn_time>0){{
                        for (time_step in 1:learn_time) {{
                            real ind_term = {ILVi_variable_effects};
                            real soc_term = s[id] * ({network_term}) {ILVs_variable_effects};
                            real lambda = {ILVm_variable_effects} lambda_0[id] * (ind_term + soc_term) * D[trial, time_step];
                            target += -lambda; //cumulative hazard (accumulating)
                            if (time_step == learn_time) {{
                                target += log({ILVm_variable_effects} lambda_0[id] * (ind_term + soc_term)); //inst. hazard
                            }}
                        }}
                    }}

                }}

                // Contributions of censored individuals
                if (N_c[trial] > 0) {{
                    for (c in 1:N_c[trial]) {{
                        int id = ind_id[trial, N[trial] + c];

                        for (time_step in 1:T[trial]) {{
                            real ind_term = {ILVi_variable_effects};
                            real soc_term = s[id] * ({network_term}) {ILVs_variable_effects};
                            real lambda = {ILVm_variable_effects} lambda_0[id] * (ind_term + soc_term) * D[trial, time_step];
                            target += -lambda;
                        }}
                    }}
                }}
            }}
        }}
        ")
    }

    else{
        # Build the Stan model
        stan_model <- glue::glue("
        data {{
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
            array[K, T_max] matrix[Z, Z] {paste0('A_', network_names, collapse = ', ')}; // Network matrices
            array[K] matrix[T_max, Z] C;        // knowledge state slash cue matrix
            {ILV_declaration}
        }}

        parameters {{
            real log_lambda_0_mean;  // Log baseline learning rate
            real log_s_mean; // Overall social transmission rate
            {w_param}
            {ILVi_param}
            {ILVs_param}
            {ILVm_param}

        }}

        transformed parameters {{
            real<lower=0> lambda_0 = 1 / exp(log_lambda_0_mean);
            real<lower=0> s = exp(log_s_mean);
        }}

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

                    if (learn_time>0){{
                        for (time_step in 1:learn_time) {{
                            real ind_term = {ILVi_variable_effects};
                            real soc_term = s * ({network_term}) {ILVs_variable_effects};
                            real lambda = {ILVm_variable_effects} lambda_0 * (ind_term + soc_term) * D[trial, time_step];
                            target += -lambda; //cumulative hazard (accumulating)
                            if (time_step == learn_time) {{
                                target += log({ILVm_variable_effects} lambda_0 * (ind_term + soc_term)); //inst. hazard
                            }}
                        }}
                    }}
                }}

                // Contributions of censored individuals
                if (N_c[trial] > 0) {{
                    for (c in 1:N_c[trial]) {{
                        int id = ind_id[trial, N[trial] + c];

                        for (time_step in 1:T[trial]) {{
                            real ind_term = {ILVi_variable_effects};
                            real soc_term = s * ({network_term}) {ILVs_variable_effects};
                            real lambda = {ILVm_variable_effects} lambda_0 * (ind_term + soc_term) * D[trial, time_step];
                            target += -lambda;
                        }}
                    }}
                }}
            }}
        }}
        ")
    }

    return(stan_model)
}
