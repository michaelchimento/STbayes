#' Dynamically generate STAN model based on input data
#'
#' @param STb_data a list of formatted data returned from the STbayes_data() function
#' @param model_type string specifying the model type: "asocial" or "full"
#' @param transmission_func string specifying transmission function: "standard", "freqdep_f" or "freqdep_k" for frequency dependent complex contagion. Defaults to "standard".
#' @param veff_ID Parameters for which to estimate varying effects by individuals. Default is no varying effects.
#' @param gq Boolean to indicate whether the generated quantities block is added (incl. ll for WAIC)
#' @param priors named list with strings containing the prior for log s or f. defaults to list(log_s = "uniform(-10, 10)", log_f = "normal(0,1)")

#' @return A STAN model (character) that is customized to the input data.
#'
#' @examples
#' #very mock data
#' event_data <- data.frame(
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
#'   event_data = event_data,
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
#'
generate_STb_model_OADA <- function(STb_data,
                                    model_type="full",
                                    transmission_func="standard",
                                    multi_network_s = c("shared", "separate"),
                                    veff_ID = c(),
                                    gq = TRUE,
                                    priors = list()) {
    #maybe one day this can be an argument :(
    est_acqTime = FALSE

    if (!model_type %in% c("asocial", "full")) {
        stop("Invalid model_type. Choose 'asocial' or 'full'.")
    }

    prior_baserate <- priors[["log_lambda_0"]]
    prior_s <- priors[["log_s"]]
    prior_f <- priors[["log_f"]]
    prior_k <- priors[["k_raw"]]
    prior_z_ID <- priors[["z_ID"]]
    prior_sigma_ID <- priors[["sigma_ID"]]
    prior_rho_ID <- priors[["rho_ID"]]
    prior_beta <- priors[["beta_ILV"]]

    # check if edgeweights are sampled from posterior distribution (import func should have created S variable)
    if ("N_dyad" %in% names(STb_data)) is_distribution <- TRUE else is_distribution <- FALSE
    if ("N_dyad" %in% names(STb_data)) est_acqTime <- FALSE # don't want to deal with that rn

    network_names <- STb_data$network_names
    num_networks <- length(network_names)
    multi_network_s <- match.arg(multi_network_s)
    separate_s <- (multi_network_s == "separate" & num_networks > 1)

    if (multi_network_s == "separate" & num_networks == 1) {
      stop("multi_network_s = 'separate' requires more than one network.")
    }

    # make custom declarations for distributions:
    if (is_distribution & model_type=="full") {
      # data declaration
      distribution_data_declaration <- glue::glue("    matrix[N_networks, N_dyad] logit_edge_mu;  // logit edge values
                      array[N_networks] matrix[N_dyad, N_dyad] logit_edge_cov;  // covariance matrix")
      distribution_data_declaration <- glue::glue("int N_dyad;  // number of dyads\n\n{distribution_data_declaration}")

      # param declaration
      distribution_param_declaration <-
        glue::glue_collapse(glue::glue("matrix[N_networks, N_dyad] edge_logit;"), sep = "\n")

      # transformed param declaration
      distribution_transformed_declaration <- {
        # Declare matrices
        matrix_decls <- glue::glue_collapse(glue::glue("array[N_networks] matrix[P, P] A;"), sep = "\n")

        glue::glue("{matrix_decls}
                      {{
                        for (network in 1:N_networks){{
                        int edge_idx = 1;
                            for (i in 1:P) {{
                              for (j in 1:P) {{
                                if (i != j) {{
                                    A[network, i, j] = inv_logit(edge_logit[network, edge_idx]);
                                  edge_idx += 1;
                                }} else {{
                                    A[network, i, j]  = 0;
                                }}
                              }}
                            }}
                        }}

                      }}")
      }
      # model declaration
      distribution_model_block <- "
        for (n in 1:N_networks) {
            edge_logit[n] ~ multi_normal(logit_edge_mu[n], logit_edge_cov[n]);
        }"

    } else {
      distribution_data_declaration <- ""
      distribution_param_declaration <- ""
      distribution_transformed_declaration <- ""
      distribution_model_block <- ""
    }

    if (model_type == "full") {
      if (separate_s) {
        if ("s" %in% veff_ID) {
          # Varying s per network and individual
          s_param <- paste0("matrix[N_networks, P] log_s_mean;")
        } else {
          # Static s per network
          s_param <- "vector[N_networks] log_s_mean;"
        }
      } else {
        if ("s" %in% veff_ID) {
          # Varying s per individual, shared network weights
          s_param <- "real log_s_mean;"
        } else {
          # Static global s
          s_param <- "real log_s_mean;"
        }
      }
    } else {
      s_param <- ""
    }


    # check if user wants to fit f parameter w complex contagion
    if (transmission_func == "freqdep_f") {
      f_param <- "real log_f_mean;"
      f_prior <- paste0("log_f_mean ~ ", prior_f, ";")
      if (is.element("f", veff_ID)) f_statement <- "f[id]" else f_statement <- "f"
    } else {
      f_param <- ""
      f_prior <- ""
    }

    # check if user wants to fit k parameter w complex contagion
    if (transmission_func == "freqdep_k") {
      k_param <- "real k_raw;"
      k_prior <- paste0("k_raw ~ ", prior_k, ";")
      if (is.element("k", veff_ID)) k_statement <- "k_shape[id]" else k_statement <- "k_shape"
    } else {
      k_param <- ""
      k_prior <- ""
    }

    # Declare network variables and weight parameter if multi-network (only for social model)
    if (model_type == "full") {
      if (!is_distribution) {
        network_declaration <- "array[N_networks, K, T_max] matrix[P, P] A;  // network matrices"
      } else {
        network_declaration <- ""
      }

      network_term <- get_network_term(
        transmission_func = transmission_func,
        is_distribution = is_distribution,
        num_networks = num_networks,
        separate_s = separate_s,
        veff_ID = veff_ID,
        id_var = "id"
      )

      network_term_j <- get_network_term(
        transmission_func = transmission_func,
        is_distribution = is_distribution,
        num_networks = num_networks,
        separate_s = separate_s,
        veff_ID = veff_ID,
        id_var = "j"
      )

      # If shared s (i.e., use w[] weights), declare w
      if (!separate_s & num_networks>1) {
        w_param <- paste0("simplex[", num_networks, "] w; // Weights for networks")
        w_prior <- paste0("w ~ dirichlet(rep_vector(0.5, ", num_networks, "));")
      } else {
        w_param <- ""
        w_prior <- ""
      }
    } else{
      network_term=""
      network_term_j =""
    }

    # Declare variables that will be used for ILVs
    ILVi_vars = STb_data$ILVi_names[!STb_data$ILVi_names %in% "ILVabsent"]
    ILVs_vars = if (model_type == "full") STb_data$ILVs_names[!STb_data$ILVs_names %in% "ILVabsent"] else character(0)
    ILVm_vars = if (model_type == "full") STb_data$ILVm_names[!STb_data$ILVm_names %in% "ILVabsent"] else character(0)

    ILVi_vars_clean = ILVi_vars
    ILVs_vars_clean = ILVs_vars
    ILVm_vars_clean= ILVm_vars

    # placeholders for dynamic components
    num_ILVi = length(ILVi_vars)
    num_ILVs = length(ILVs_vars)
    num_ILVm = length(ILVm_vars)

    combined_ILV_vars = unique(c(ILVi_vars, ILVs_vars, ILVm_vars))

    ILV_declaration = ""
    if (length(combined_ILV_vars) > 0) {
        ILV_declaration = paste0(
            sapply(combined_ILV_vars, function(var) {
                if (!is.null(dim(STb_data[[paste0("ILV_", var)]]))) {
                    paste0('array[K,T_max,P] real ILV_', var, ';')
                } else {
                    paste0('array[P] real ILV_', var, ';')
                }
            }),
            collapse = '\n'
        )
    }

    # deal with varying effects in transformed parameters

    N_veff = length(veff_ID)
    #start w empty list
    transformed_params <- c()
    gq_transformed_params <- c()

    #if user didn't specify veff_ID for s
    if (!is.element('s', veff_ID)  & model_type=="full"){
      if (separate_s) {
        # static s[n]
        transformed_params <- append(transformed_params, "vector<lower=0>[N_networks] s_prime = exp(log_s_mean);")
        #gq_transformed_params <- append(gq_transformed_params, "vector<lower=0> s = s_prime ./ lambda_0;")
      } else {
        # static scalar s
        transformed_params <- append(transformed_params, "real<lower=0> s_prime = exp(log_s_mean);")
        #gq_transformed_params <- append(gq_transformed_params, "real<lower=0> s = s_prime/lambda_0;")
      }
    }
    #if user didn't specify veff_ID for f
    if (!is.element('f', veff_ID) & model_type=="full" & transmission_func=="freqdep_f"){
        transformed_params = append(transformed_params,"real<lower=0> f = exp(log_f_mean);")
    }
    #if user didn't specify veff_ID for k
    if (!is.element('k', veff_ID) & model_type=="full" & transmission_func=="freqdep_k"){
        transformed_params = append(transformed_params,"real<lower=-1, upper=1> k_shape = 2 / (1 + exp(-k_raw)) - 1;")
    }

    count = 1

    if (N_veff > 0){
        for (parameter in veff_ID) {
          if (parameter == "s" & model_type == "full") {
            #this is still a bit weird, adding single varying effect for all networks..
            if (separate_s) {
              # s[n, id] as exp(log_s_mean[n, id] + v_ID[,i])
              transformed_params <- append(transformed_params, glue::glue(
                "matrix<lower=0>[N_networks, P] s_prime;
for (n in 1:N_networks) {{
  for (id in 1:P) {{
    s_prime[n, id] = exp(log_s_mean[n] + v_ID[id, {count}]);
  }}
}}"
              ))
              gq_transformed_params <- append(gq_transformed_params, "real s_mean = mean(to_vector(s_prime));")
            } else {
              # s[id] = exp(log_s_mean + v_ID[,i])
              transformed_params <- append(transformed_params, paste0("vector<lower=0>[P] s_prime = exp(log_s_mean + v_ID[,", count, "]);"))
              #gq_transformed_params <- append(gq_transformed_params, paste0("vector<lower=0>[P] s = s_prime ./ lambda_0;"))
              gq_transformed_params <- append(gq_transformed_params, paste0("real sprime_mean = exp(log_s_mean);"))
              #gq_transformed_params <- append(gq_transformed_params, paste0("real<lower=0> s_mean = (exp(log_s_mean)) / (exp(log_lambda_0_mean));"))
            }
            count <- count + 1
          }
            else if (parameter=="f" & model_type=="full"){
                transformed_params = append(transformed_params, paste0("vector<lower=0>[P] f = exp(log_f_mean + v_ID[,",count,"]);"))
                transformed_params = append(transformed_params, paste0("real<lower=0> f_mean = exp(log_f_mean);"))
            }
            else if (parameter=="k" & model_type=="full"){
                transformed_params = append(transformed_params, paste0("vector<lower=0>[P] k_shape = 2 / (1 + exp(-k_raw+ v_ID[,",count,"])) - 1;"))
                transformed_params = append(transformed_params, paste0("real<lower=-1, upper=1> k_shape_mean = 2 / (1 + exp(-k_raw)) - 1;"))
            }
            count = count + 1
        }
    }

    # Handle asocial ILV (ILVi)
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
                transformed_params = append(transformed_params, paste0("vector[P] ", ilv, " = beta_ILVi_", ilv," + v_ID[,",count,"]);"))
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
        ILVi_prior <- paste0("beta_ILVi_", ILVi_vars_clean, " ~ ",prior_beta,";", sep = "\n")
        ILVi_variable_effects <- paste0(
            "exp(",
            paste0(
                ILVi_vars,
                " * ",
                sapply(ILVi_vars_clean, function(var_clean) {
                    if (!is.null(dim(STb_data[[paste0("ILV_", var_clean)]]))) {
                        paste0("ILV_", var_clean,"[trial,time_step,id]") #if has dimensions, assume its time-varying
                    } else {
                        paste0("ILV_", var_clean,"[id]") #else time constant
                    }
                }),
                collapse = " + "
            ),")")
    }

    # Handle social ILV (ILVs)
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
                transformed_params = append(transformed_params, paste0("vector[P] ", ilv, " = beta_ILVs_", ilv," + v_ID[,",count,"]);"))
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
        ILVs_prior <- paste0("beta_ILVs_", ILVs_vars_clean, " ~ ",prior_beta,";", sep = "\n")
        ILVs_variable_effects <- paste0(
            "* exp(",
            paste0(
                ILVs_vars,
                " * ",
                sapply(ILVs_vars_clean, function(var_clean) {
                    if (!is.null(dim(STb_data[[paste0("ILV_", var_clean)]]))) {
                        paste0("ILV_", var_clean,"[trial,time_step,id]")
                    } else {
                        paste0("ILV_", var_clean,"[id]")
                    }
                }),
                collapse = " + "
            ),")")
    }

    # Handle multiplicative ILV (ILVm)
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
                transformed_params = append(transformed_params, paste0("vector[P] ", ilv, " = beta_ILVm_", ilv," + v_ID[,",count,"]);"))
                count = count + 1
                #rename with [id] so it can be indexed in the main model loop
                ILVm_vars[ILVm_vars == ilv] <- paste0(ilv, "[id]")
            } else {
                #otherwise append prefix of beta_ILVx_ so it can be accessed directly
                ILVm_vars[ILVm_vars == ilv] <- paste0("beta_ILVm_", ilv)
            }
        }
        ILVm_param <- paste0("real beta_ILVm_", ILVm_vars_clean, ";", sep = "\n")
        ILVm_prior <- paste0("beta_ILVm_", ILVm_vars_clean, " ~ ",prior_beta,";", sep = "\n")
        ILVm_variable_effects <- paste0(
            "exp(",
            paste0(
                ILVm_vars,
                " * ",
                sapply(ILVm_vars_clean, function(var_clean) {
                    # Check if dim() is NULL for the variable
                    if (!is.null(dim(STb_data[[paste0("ILV_", var_clean)]]))) {
                        paste0("ILV_", var_clean,"[trial,time_step,id]")
                    } else {
                        paste0("ILV_", var_clean,"[id]")
                    }
                }),
                collapse = " + "
            ),") *")
    }

    #for oada we'll have to index j for part of the likelihood
    #don't gsub "id" in case id is part of a variable name
    ILVi_variable_effects_j = gsub("[id]", "[j]", ILVi_variable_effects, fixed = TRUE)
    ILVi_variable_effects_j = gsub("[trial,time_step,id]", "[trial,time_step,j]", ILVi_variable_effects_j, fixed = TRUE)

    ILVs_variable_effects_j = gsub("[id]", "[j]", ILVs_variable_effects, fixed = TRUE)
    ILVs_variable_effects_j = gsub("[trial,time_step,id]", "[trial,time_step,j]", ILVs_variable_effects_j, fixed = TRUE)

    ILVm_variable_effects_j = gsub("[id]", "[j]", ILVm_variable_effects, fixed = TRUE)
    ILVm_variable_effects_j = gsub("[trial,time_step,id]", "[trial,time_step,j]", ILVm_variable_effects_j, fixed = TRUE)

    #collapse lists into multiline statements
    ILVi_param = paste0(ILVi_param, collapse = "\n")
    ILVi_prior = paste0(ILVi_prior, collapse = "\n")
    ILVs_param = paste0(ILVs_param, collapse = "\n")
    ILVs_prior = paste0(ILVs_prior, collapse = "\n")
    ILVm_param = paste0(ILVm_param, collapse = "\n")
    ILVm_prior = paste0(ILVm_prior, collapse = "\n")
    transformed_params_declaration = paste0(transformed_params, collapse = "\n")
    gq_transformed_params_declaration <- paste0(gq_transformed_params, collapse = "\n")

    N_veff_priors = if (N_veff > 0) glue::glue("
    to_vector(z_ID) ~ {prior_z_ID};
    sigma_ID ~ {prior_sigma_ID};
    Rho_ID ~ {prior_rho_ID};") else ''

    functions_block <- if (transmission_func=="freqdep_k") glue::glue("
functions {{
  real dini_func(real x, real k) {{
    // transform x from [0,1] to [-1,1]
    real x_transformed = 2 * x - 1;
    real y = ((x_transformed - k * x_transformed) / (k - 2 * k * fabs(x_transformed) + 1) + 1) / 2;
    return y;
  }}
}}") else ""

    data_block <- glue::glue("
data {{
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
    {if (model_type=='full') {network_declaration} else ''}
    array[K] matrix[T_max, P] Z;   // Knowledge state * cue matrix
    array[K] matrix[T_max, P] Zn;   // Knowledge state
    {ILV_declaration}
    int<lower=0> N_veff;
    {if (est_acqTime) 'array[K] int<lower=0> time_max; //Duration of obs period for each trial' else ''}
    {if (est_acqTime) 'array[K, T_max] int<lower=0> D_int; // integer durations' else ''}
    {distribution_data_declaration}

}}
")

    # Parameters block
    parameters_block <- glue::glue("
parameters {{
    {distribution_param_declaration}
    {if (model_type=='full') 'real log_s_mean; // Overall social transmission rate' else ''}
    {if (model_type=='full') {w_param} else ''}
    {ILVi_param}
    {if (model_type=='full') {ILVs_param} else ''}
    {if (model_type=='full') {ILVm_param} else ''}
    {if (N_veff > 0) '
    matrix[N_veff,P] z_ID;
    vector<lower=0, upper=2>[N_veff] sigma_ID;
    cholesky_factor_corr[N_veff] Rho_ID;
    ' else ''}
}}
")

    # Transformed parameters block
    transformed_parameters_block <- glue::glue("
transformed parameters {{
   {if (N_veff > 0) '
    matrix[P,N_veff] v_ID;
    v_ID = (diag_pre_multiply(sigma_ID, Rho_ID) * z_ID)\\';
   ' else ''}
    {transformed_params_declaration}
    {distribution_transformed_declaration}
}}
")

    #create string inputs cuz recursion don't work 2 levels down in glue
    if (model_type=='full'){
        i_social_info_statement = glue::glue("real i_soc = {if (is.element('s', veff_ID) & !separate_s) 's_prime[id]' else if (!is.element('s', veff_ID) & !separate_s) 's_prime' else '1.0'} * (net_effect{ILVs_variable_effects});")
        i_lambda_statement = glue::glue("real i_lambda = {ILVm_variable_effects} * (i_ind + i_soc);")
        j_social_info_statement = glue::glue("real j_soc = {if (is.element('s', veff_ID) & !separate_s) 's_prime[j]' else if (!is.element('s', veff_ID) & !separate_s) 's_prime' else '1.0'} * (net_effect_j{ILVs_variable_effects_j});")
        j_lambda_statement = glue::glue("real j_lambda = {ILVm_variable_effects_j} * (j_ind + j_soc);")

        target_increment_statement = glue::glue("target += log(i_lambda) - log(sum(j_rates));")
        log_lik_statement = glue::glue("log_lik_matrix[trial, n] = log(i_lambda) - log(sum(j_rates));")

    } else if (model_type=='asocial'){
        i_social_info_statement = ""
        i_lambda_statement = glue::glue("real i_lambda = {ILVm_variable_effects} * i_ind;")
        j_social_info_statement = ""
        j_lambda_statement = glue::glue("real j_lambda = {ILVm_variable_effects_j} * j_ind;")
        target_increment_statement = glue::glue("target += log(i_lambda) - log(sum(j_rates));")
        log_lik_statement = glue::glue("log_lik_matrix[trial, n] = log(i_lambda) - log(sum(j_rates));")
    }

    # Model block
    model_block <- glue::glue("
model {{
    {if (model_type=='full') paste0('log_s_mean ~ ',prior_s,';') else ''}
    {if (model_type=='full') {w_prior} else ''}
    {if (model_type=='full') {f_prior} else ''}
    {if (model_type=='full') {k_prior} else ''}
    {ILVi_prior}
    {if (model_type=='full') {ILVs_prior} else ''}
    {if (model_type=='full') {ILVm_prior} else ''}
    {distribution_model_block}

    {N_veff_priors}

    for (trial in 1:K) {{

        for (n in 1:N[trial]) {{
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            int time_step = learn_time;
            if (learn_time > 0) {{
                real i_ind = {ILVi_variable_effects};
                {network_term}
                {i_social_info_statement}
                {i_lambda_statement}

                vector[Q] j_rates = rep_vector(0.0, Q);

                for (j in 1:Q) {{
                    real j_ind = {ILVi_variable_effects_j};
                    {network_term_j}
                    {j_social_info_statement}
                    {j_lambda_statement}
                    j_rates[j] += j_lambda * (1-Z[trial][learn_time, j]); //only include those who haven't learned in denom
                }}
                {target_increment_statement}
            }}
        }}
    }}
}}
")

#I'll have to come back to this, if it's even possible. as it is, the rng can pick the same ind twice, and Z won't match the simulated order
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
#                 real j_ind = {ILVi_variable_effects_j};
#                 {j_social_info_statement}
#                 {j_lambda_statement}
#                 j_rates[j] += j_lambda * (1-Z[trial][time_step, j]);
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
    {gq_transformed_params_declaration}
    matrix[K, Q] log_lik_matrix = rep_matrix(0.0, K, Q);           // LL for each observation

    for (trial in 1:K) {{
        for (n in 1:N[trial]) {{
                int id = ind_id[trial, n];
                int learn_time = t[trial, id];
                int time_step = learn_time;
                if (learn_time > 0) {{
                    real i_ind = {ILVi_variable_effects};
                    {network_term}
                    {i_social_info_statement}
                    {i_lambda_statement}

                    vector[Q] j_rates = rep_vector(0.0, Q);

                    for (j in 1:Q) {{
                        real j_ind = {ILVi_variable_effects_j};
                        {network_term_j}
                        {j_social_info_statement}
                        {j_lambda_statement}
                        j_rates[j] += j_lambda * (1-Z[trial][learn_time, j]);
                    }}
                    {log_lik_statement}
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
    stan_model <- glue::glue("//stan
                              {functions_block}
                              {data_block}
                             {parameters_block}
                             {transformed_parameters_block}
                             {model_block}
                             {if (gq==T) {generated_quantities_block} else ''}")
    stan_model = gsub("(?m)^[ \\t]*\\n", "", stan_model, perl = TRUE)
    stan_model <- paste0(stan_model, "\n")
    return(stan_model)
}
