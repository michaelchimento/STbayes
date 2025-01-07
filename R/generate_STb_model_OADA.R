#' Dynamically generate STAN model based on input data
#'
#' @param STb_data a list of formatted data returned from the STbayes_data() function
#' @param veff_ID Parameters for which to estimate varying effects by individuals. Default is no varying effects.
#' @param gq Boolean to indicate whether the generated quantities block is added (incl. ll for WAIC)
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
#'   sex = c(0, 1, 1, 0, 1, 0), # Factor ILVs must be input as numeric
#'   weight = c(0.5, .25, .3, 0, -.2, -.4)
#' )
#' data_list <- import_user_STb(
#'   diffusion_data = diffusion_data,
#'   networks = networks,
#'   ILV_metadata = ILV_metadata,
#'   ILVi = c("age"), # Use only 'age' for asocial learning
#'   ILVs = c("sex"), # Use only 'sex' for social learning
#'   ILVm = c("weight") # Use weight for multiplicative effect on asocial and social learning
#' )
#'
#' model = generate_STb_model_OADA(data_list) # no varying effects
#' model = generate_STb_model_OADA(data_list, veff_ID = c("s")) # estimate varying effects by ID for strength of social transmission. Baseline is not estimated for OADA type models.
#' print(model)
generate_STb_model_OADA <- function(STb_data, veff_ID = c(), gq = TRUE, est_acqTime = FALSE) {

    # declare network variables and weight parameter if multi-network
    network_names = STb_data$network_names
    num_networks <- length(STb_data$network_names)
    if (num_networks == 1) {
        network_term <- paste0("sum(A_", network_names[1], "[trial, time_step][id, ] .* C[trial][time_step, ])")
        network_term_j <- paste0("sum(A_", network_names[1], "[trial, time_step][j, ] .* C[trial][time_step, ])")
        w_param <- ""
        w_prior <- ""
    } else {
        network_term <- paste0("w[", 1:num_networks, "] * sum(A_", network_names, "[trial, time_step][id, ] .* C[trial][time_step, ])", collapse = " + ")
        network_term_j <- paste0("w[", 1:num_networks, "] * sum(A_", network_names, "[trial, time_step][J, ] .* C[trial][time_step, ])", collapse = " + ")
        w_param <- paste0("simplex[", num_networks, "] w; // Weights for networks")
        w_prior <- paste0("w ~ dirichlet(rep_vector(0.5, ", num_networks, "));")
    }



    # Declare variables that will be used for ILVs
    ILVi_vars = STb_data$ILVi_names[!STb_data$ILVi_names %in% "ILVabsent"]
    ILVs_vars = STb_data$ILVs_names[!STb_data$ILVs_names %in% "ILVabsent"]
    ILVm_vars = STb_data$ILVm_names[!STb_data$ILVm_names %in% "ILVabsent"]

    ILVi_vars_clean = ILVi_vars
    ILVs_vars_clean = ILVs_vars
    ILVm_vars_clean= ILVm_vars

    # placeholders for dynamic components
    num_ILVi = length(ILVi_vars)
    num_ILVs = length(ILVs_vars)
    num_ILVm = length(ILVm_vars)

    combined_ILV_vars = unique(c(ILVi_vars, ILVs_vars, ILVm_vars))

    ILV_declaration = ""
    if (length(combined_ILV_vars)>0){
        ILV_declaration = paste0('array[Z] real ', combined_ILV_vars, ';', collapse = '\n')
    }

    # deal with varying effects in transformed parameters

    N_veff = length(veff_ID)
    #start w empty list
    transformed_params = c()

    #if user didn't specify veff_ID for s
    if (!is.element('s', veff_ID)){
        transformed_params = append(transformed_params,"real<lower=0> s = exp(log_s_mean);")
    }
    count = 1

    if (N_veff > 0){
        for (parameter in veff_ID) {
            if (parameter=="s"){
                transformed_params = append(transformed_params, paste0("vector<lower=0>[Z] s = exp(log_s_mean + v_ID[,",count,"]);"))
                count = count + 1
            }
        }
    }

    # Handle individual-level information (ILVi)
    ILVi_variable_effects = c()
    if (num_ILVi < 1) {
        ILVi_param <- ""
        ILVi_prior <- ""
        ILVi_variable_effects <- "1.0"
    } else {
        #for each ilv
        for (ilv in ILVi_vars){
            #if user specified this should be include a varying effect for id
            if (is.element(ilv, veff_ID)){
                #add declaration in transformed parameters
                transformed_params = append(transformed_params, paste0("vector<lower=0>[Z] ", ilv, " = beta_ILVi_", parameter," + v_ID[,",count,"]);"))
                count = count + 1
                #rename with [id] so it can be indexed in the main model loop
                ILVi_vars[ILVi_vars == ilv] <- paste0(ilv, "[id]")
            } else {
                #otherwise append prefix of beta_ILVx_ so it can be accessed directly
                ILVi_vars[ILVi_vars == ilv] <- paste0("beta_ILVi_", ilv)
            }
        }
        #creates the line to insert into likelihood calculation
        ILVi_param <- paste0("real beta_ILVi_", ILVi_vars_clean, ";", sep = "\n")
        ILVi_prior <- paste0("beta_ILVi_", ILVi_vars_clean, " ~ normal(0, 1);", sep = "\n")
        ILVi_variable_effects <- paste0("exp(", paste0(ILVi_vars, " * ", ILVi_vars_clean, "[id]", collapse = " + "),")")
    }

    # Handle social-level information (ILVs)
    if (num_ILVs < 1) {
        ILVs_param <- ""
        ILVs_prior <- ""
        ILVs_variable_effects <- ""
    } else {

        #for each ilv
        for (ilv in ILVs_vars){
            #if user specified this should be include a varying effect for id
            if (is.element(ilv, veff_ID)){
                #add declaration in transformed parameters
                transformed_params = append(transformed_params, paste0("vector<lower=0>[Z] ", ilv, " = beta_ILVs_", parameter," + v_ID[,",count,"]);"))
                count = count + 1
                #rename with [id] so it can be indexed in the main model loop
                ILVs_vars[ILVs_vars == ilv] <- paste0(ilv, "[id]")
            } else {
                #otherwise append prefix of beta_ILVx_ so it can be accessed directly
                ILVs_vars[ILVs_vars == ilv] <- paste0("beta_ILVs_", ilv)
            }
        }
        #creates the line to insert into likelihood calculation
        ILVs_param <- paste0("real beta_ILVs_", ILVs_vars_clean, ";", sep = "\n")
        ILVs_prior <- paste0("beta_ILVs_", ILVs_vars_clean, " ~ normal(0, 1);", sep = "\n")
        ILVs_variable_effects <- paste0("* exp(", paste0(ILVs_vars, " * ", ILVs_vars_clean, "[id]", collapse = " + "),")")
    }

    # Handle social-level information (ILVs)
    if (num_ILVm < 1) {
        ILVm_param <- ""
        ILVm_prior <- ""
        ILVm_variable_effects <- "1.0"
    } else {
        #for each ilv
        for (ilv in ILVm_vars){
            #if user specified this should be include a varying effect for id
            if (is.element(ilv, veff_ID)){
                #add declaration in transformed parameters
                transformed_params = append(transformed_params, paste0("vector<lower=0>[Z] ", ilv, " = beta_ILVm_", parameter," + v_ID[,",count,"]);"))
                count = count + 1
                #rename with [id] so it can be indexed in the main model loop
                ILVm_vars[ILVm_vars == ilv] <- paste0(ilv, "[id]")
            } else {
                #otherwise append prefix of beta_ILVx_ so it can be accessed directly
                ILVm_vars[ILVm_vars == ilv] <- paste0("beta_ILVm_", ilv)
            }
        }
        ILVm_param <- paste0("real beta_ILVm_", ILVm_vars_clean, ";", sep = "\n")
        ILVm_prior <- paste0("beta_ILVm_", ILVm_vars_clean, " ~ normal(0, 1);", sep = "\n")
        ILVm_variable_effects <- paste0("exp(", paste0(ILVm_vars, " * ", ILVm_vars_clean, "[id]", collapse = " + "),") *")
    }

    transformed_params_declaration = paste0(transformed_params, collapse = "\n")

    data_block <- glue::glue("
data {{
    int<lower=0> K;                // Number of trials
    int<lower=0> Q;                // Number of individuals in each trial
    int<lower=1> Z;                // Number of unique individuals
    array[K] int<lower=0> N;       // Number of individuals that learned during observation period
    array[K] int<lower=0> N_c;     // Number of right-censored individuals
    array[K, Q] int<lower=-1> ind_id; // IDs of individuals
    array[K] int<lower=1> T;       // Maximum time periods
    int<lower=1> T_max;            // Max timesteps reached
    {if (est_acqTime) 'array[K] int<lower=0> time_max; //Duration of obs period for each trial' else ''}
    array[K,Z] int<lower=-1> t;     // Time of acquisition for each individual
    array[K, T_max] matrix[Z, Z] {paste0('A_', network_names, collapse = ', ')}; // Network matrices
    array[K] matrix[T_max, Z] C;   // Knowledge state slash cue matrix
    int<lower=0> N_veff;
    {ILV_declaration}
}}
")

    # Parameters block
    parameters_block <- glue::glue("
parameters {{
    real log_s_mean;         // Overall social transmission rate
    {w_param}
    {ILVi_param}
    {ILVs_param}
    {ILVm_param}
    {if (N_veff > 0) '
    matrix[N_veff,Z] z_ID;
    vector<lower=0>[N_veff] sigma_ID;
    cholesky_factor_corr[N_veff] Rho_ID;
    ' else ''}
}}
")

    # Transformed parameters block
    transformed_parameters_block <- glue::glue("
transformed parameters {{
   {if (N_veff > 0) '
    matrix[Z,N_veff] v_ID;
    v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)\\';
   ' else ''}
   {transformed_params_declaration}
}}
")

    # Model block
    model_block <- glue::glue("
model {{
    log_s_mean ~ normal(0,2);
    {w_prior}
    {ILVi_prior}
    {ILVs_prior}
    {ILVm_prior}

    {if (N_veff > 0) '
    to_vector(z_ID) ~ normal(0,1);
    sigma_ID ~ exponential(1);
    Rho_ID ~ lkj_corr_cholesky(1);
    ' else ''}

    for (trial in 1:K) {{
        for (n in 1:N[trial]) {{
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            int time_step = learn_time;
            if (learn_time > 0) {{
                real i_ind = {ILVi_variable_effects};
                real i_soc = {if (is.element('s', veff_ID)) 's[id]' else 's'} * ({network_term}) {ILVs_variable_effects};
                real i_lambda = {ILVm_variable_effects} * (i_ind + i_soc);

                vector[Q] j_rates = rep_vector(0.0, Q);

                for (j in 1:Q) {{
                    real j_ind = {ILVi_variable_effects};
                    real j_soc = {if (is.element('s', veff_ID)) 's[j]' else 's'} * ({network_term_j}) {ILVs_variable_effects};
                    real j_lambda = {ILVm_variable_effects} * (j_ind + j_soc);
                    j_rates[j] += j_lambda * (1-C[trial][learn_time, j]); //only include those who haven't learned in denom
                }}
                target += log(i_lambda) - log(sum(j_rates));
            }}
        }}
    }}
}}
")

#I'll have to come back to this, if it's even possible. as it is, the rng can pick the same ind twice, and C won't match the simulated order
#if this ever gets uncommented, be sure to add back {est_acqTime_code} into the generated quantities glue.
# est_acqTime_code <- if (est_acqTime==TRUE) glue::glue("
#     matrix[K, Q] acquisition_time;         // simulated acquisition times
#     for (trial in 1:K) {{
#         for (n in 1:Q) {{ //have to loop through bc stan
#             acquisition_time[trial, n] = time_max[trial];
#         }}
#
#         int global_time = 1;
#         for (time_step in 1:T[trial]) {{
#             vector[Q] j_rates = rep_vector(0.0, Q);
#             for (j in 1:Q) {{
#                 real j_ind = {ILVi_variable_effects};
#                 real j_soc = {if (is.element('s', veff_ID)) 's[j]' else 's'} * ({network_term_j}) {ILVs_variable_effects};
#                 real j_lambda = {ILVm_variable_effects} * (j_ind + j_soc);
#                 j_rates[j] += j_lambda * (1-C[trial][time_step, j]);
#             }}
#
#             vector[Q] probs = j_rates / sum(j_rates); // normalize to probabilities
#             int learner = categorical_rng(probs); // sample index proportional to probabilities
#             acquisition_time[trial, learner] = time_step;
#             global_time += 1;
#         }}
#     }}"
# ) else ""

    generated_quantities_block <- glue::glue("
generated quantities {{
    matrix[K, Q] log_lik_matrix = rep_matrix(0.0, K, Q);           // LL for each observation

    for (trial in 1:K) {{
        for (n in 1:N[trial]) {{
                int id = ind_id[trial, n];
                int learn_time = t[trial, id];
                int time_step = learn_time;
                if (learn_time > 0) {{
                    real i_ind = {ILVi_variable_effects};
                    real i_soc = {if (is.element('s', veff_ID)) 's[id]' else 's'} * ({network_term}) {ILVs_variable_effects};
                    real i_lambda = {ILVm_variable_effects} * (i_ind + i_soc);

                    vector[Q] j_rates = rep_vector(0.0, Q);

                    for (j in 1:Q) {{
                        real j_ind = {ILVi_variable_effects};
                        real j_soc = {if (is.element('s', veff_ID)) 's[j]' else 's'} * ({network_term_j}) {ILVs_variable_effects};
                        real j_lambda = {ILVm_variable_effects} * (j_ind + j_soc);
                        j_rates[j] += j_lambda * (1-C[trial][learn_time, j]);
                    }}
                    log_lik_matrix[trial, n] = log(i_lambda) - log(sum(j_rates));
                }}
            }}
    }}

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
