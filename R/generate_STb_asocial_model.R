#' Dynamically generate asocial model based on input data
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
#' model = generate_STb_asocial_model(data_list)
#' model = generate_STb_asocial_model(data_list, N_veff=1) # estimate varying effects lambda_0.
#' print(model)
generate_STb_asocial_model <- function(STb_data, N_veff = 0) {
    ILVi_vars = STb_data$ILVi_names[!STb_data$ILVi_names %in% "ILVabsent"]

    combined_ILV_vars = unique(c(ILVi_vars))

    ILV_declaration = ""
    if (length(combined_ILV_vars)>0){
        ILV_declaration = paste0('array[Z] real ', combined_ILV_vars, ';', collapse = '\\n')
    }

    # Placeholders for dynamic components
    num_ILVi = length(ILVi_vars)

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

    stan_model = ""

    if (N_veff==1){
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
            array[K] matrix[T_max, Z] C;        // knowledge state slash cue matrix
            {ILV_declaration}
            int<lower=0> N_veff;
        }}

        parameters {{
            real log_lambda_0_mean;  // Log baseline learning rate
            {ILVi_param}

            // ID effects
            matrix[N_veff,Z] z_ID;               // Matrix of uncorrelated z-values
            vector<lower=0>[N_veff] sigma_ID;    // SD of parameters among individuals
            cholesky_factor_corr[N_veff] Rho_ID; // Cholesky factor for correlation matrix

        }}

        transformed parameters {{
            matrix[Z,N_veff] v_ID; // Matrix of varying effects for each individual
            v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)';

            vector<lower=0>[Z] lambda_0 = 1 / exp(log_lambda_0_mean + v_ID[,1]);
        }}

        model {{
            log_lambda_0_mean ~ normal(6, 2);
            {ILVi_prior}

            for (trial in 1:K) {{
                for (n in 1:N[trial]) {{
                    int id = ind_id[trial, n];
                    int learn_time = t[trial, id];

                    if (learn_time>0){{
                        for (time_step in 1:learn_time) {{
                            real ind_term = {ILVi_variable_effects};
                            real lambda = lambda_0[id] * ind_term * D[trial, time_step];
                            target += -lambda; //cumulative hazard (accumulating)
                            if (time_step == learn_time) {{
                                target += log(lambda_0[id] * (ind_term)); //inst. hazard
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
                            real lambda = lambda_0[id] * ind_term * D[trial, time_step];
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
            array[K] matrix[T_max, Z] C;        // knowledge state slash cue matrix
            {ILV_declaration}
            int<lower=0> N_veff;
        }}

        parameters {{
            real log_lambda_0_mean;  // Log baseline learning rate
            {ILVi_param}
        }}

        transformed parameters {{
            real lambda_0 = 1 / exp(log_lambda_0_mean);
        }}

        model {{
            log_lambda_0_mean ~ normal(6, 2);
            {ILVi_prior}

            for (trial in 1:K) {{
                for (n in 1:N[trial]) {{
                    int id = ind_id[trial, n];
                    int learn_time = t[trial, id];

                    if (learn_time>0){{
                        for (time_step in 1:learn_time) {{
                            real ind_term = {ILVi_variable_effects};
                            real lambda = lambda_0 * ind_term * D[trial, time_step];
                            target += -lambda; //cumulative hazard (accumulating)
                            if (time_step == learn_time) {{
                                target += log(lambda_0 * ind_term); //inst. hazard
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
                            real lambda = lambda_0 * ind_term * D[trial, time_step];
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
