#' Dynamically generate STAN model based on input data
#' @param STb_data a list of formatted data returned from the STbayes_data() function
#' @param model_type string specifying the model type: "asocial" or "full"
#' @param transmission_func string specifying transmission function: "standard", "freqdep_f" or "freqdep_k" for complex contagion. Defaults to "standard". Complex contagion with multi-network model is not supported.
#' @param dTADA boolean indicating whether dTADA should be used.
#' @param veff_ID Parameters for which to estimate varying effects by individuals. Default is no varying effects.
#' @param gq Boolean to indicate whether the generated quantities block is added (incl. ll for WAIC)
#' @param est_acqTime Boolean to indicate whether gq block includes estimates for acquisition time. At the moment this uses 'one weird trick' to accomplish this and does not support estimates for non-integer learning times.
#' @param priors named list with strings containing the prior for log baserate, s, f, k.
#' @return A STAN model (character) that is customized to the input data.
#'
#' @examples
#' # very mock data
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
#' ILV_c <- data.frame(
#'   id = LETTERS[1:6],
#'   age = c(-1, -2, 0, 1, 2), # continuous variables should be normalized
#'   sex = c(0, 1, 1, 0, 1, 0), # Factor ILVs must be input as numeric
#'   weight = c(0.5, .25, .3, 0, -.2, -.4)
#' )
#' data_list <- import_user_STb(
#'   event_data = event_data,
#'   networks = networks,
#'   ILV_c = ILV_c,
#'   ILVi = c("age"), # Use only 'age' for asocial learning
#'   ILVs = c("sex"), # Use only 'sex' for social learning
#'   ILVm = c("weight") # Use weight for multiplicative effect on asocial and social learning
#' )
#'
#' model <- generate_STb_model(data_list) # no varying effects
#' model <- generate_STb_model(data_list, veff_ID = c("lambda_0", "s")) # estimate varying effects by ID for baseline learning rate and strength of social transmission.
#' print(model)
generate_STb_model_TADA <- function(STb_data,
                                    model_type = "full",
                                    intrinsic_rate = "constant",
                                    transmission_func = "standard",
                                    dTADA = F,
                                    veff_ID = c(),
                                    gq = TRUE,
                                    est_acqTime = FALSE,
                                    priors = list()) {
  if (!model_type %in% c("asocial", "full")) {
    stop("Invalid model_type. Choose 'asocial' or 'full'.")
  }

  if (est_acqTime == TRUE & min(check_integer(STb_data$time)) == 0) {
    message("WARNING: You have input float times, and unfortunately estimating acquisition times in the GQ block is only possible with integer times at the moment.\nThe model will be created with est_acqTime=F.")
    est_acqTime <- FALSE
  }

  high_res_k = if (STb_data$high_res==T & transmission_func=="freqdep_k") T else F
  high_res_f = if (STb_data$high_res==T & transmission_func=="freqdep_f") T else F

  prior_lambda0 <- priors[["log_lambda0"]]
  prior_s <- priors[["log_sprime"]]
  prior_f <- priors[["log_f"]]
  prior_k <- priors[["k_raw"]]
  prior_z_ID <- priors[["z_ID"]]
  prior_sigma_ID <- priors[["sigma_ID"]]
  prior_rho_ID <- priors[["rho_ID"]]
  prior_beta <- priors[["beta_ILV"]]
  prior_gamma <- priors[["gamma"]]

  # check if edgeweights are sampled from posterior distribution (import func should have created S variable)
  if ("N_dyad" %in% names(STb_data)) is_distribution <- TRUE else is_distribution <- FALSE
  if ("N_dyad" %in% names(STb_data)) est_acqTime <- FALSE # don't want to deal with that rn

  network_names <- STb_data$network_names
  num_networks <- length(network_names)

  separate_s <- (STb_data$multinetwork_s == "separate" & num_networks > 1)

  # make custom declarations for distributions:
  if (is_distribution & model_type == "full") {
    # data declaration
    distribution_data_declaration <- glue::glue("    matrix[N_networks, N_dyad] logit_edge_mu;  // logit edge values
                      array[N_networks] matrix[N_dyad, N_dyad] logit_edge_cov;  // covariance matrix
                                                    array[N_dyad] int<lower=1> from_ID;
                                                    array[N_dyad] int<lower=1> to_ID;")
    distribution_data_declaration <- glue::glue("int N_dyad;  // number of dyads\n\n{distribution_data_declaration}")

    # param declaration
    distribution_param_declaration <-
      glue::glue_collapse(glue::glue("matrix[N_networks, N_dyad] edge_logit;"), sep = "\n")

    # transformed param declaration
    distribution_transformed_declaration <- {
      # Declare matrices
      matrix_decls <- glue::glue_collapse(glue::glue("array[N_networks] matrix[P, P] A;
                                                           for (network in 1:N_networks) {{
                                                              A[network] = rep_matrix(0, P, P);
                                                            }}"), sep = "\n")

      if (STb_data$directed == T) {
        glue::glue("{matrix_decls}
                        for (network in 1:N_networks){{
                          for (edge_idx in 1:N_dyad) {{
                                real w = inv_logit(edge_logit[network, edge_idx]);
                                A[network, from_ID[edge_idx], to_ID[edge_idx]] = w;
                          }}
                        }}")
      } else {
        glue::glue("{matrix_decls}
                        for (network in 1:N_networks){{
                       for (edge_idx in 1:N_dyad) {{
                                real w = inv_logit(edge_logit[network, edge_idx]);
                                A[network, from_ID[edge_idx], to_ID[edge_idx]] = w;
                                A[network, to_ID[edge_idx], from_ID[edge_idx]] = w;
                       }}
                        }}")
      }
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
      s_param <- "vector[N_networks] log_s_prime_mean;"
    } else {
      s_param <- "real log_s_prime_mean;"
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

  # check if user wants to fit f parameter w complex contagion
  if (intrinsic_rate == "weibull") {
    gamma_param <- "real log_gamma;"
    gamma_prior <- paste0("log_gamma ~ ", prior_gamma, ";")
    if (is.element("gamma", veff_ID)) gamma_term <- "gamma[id]" else gamma_term <- "gamma"
  } else {
    gamma_param <- ""
    gamma_prior <- ""
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
      high_res = high_res_k
    )

    # If shared s (i.e., use w[] weights), declare w
    # if (!separate_s & num_networks>1) {
    #     w_param <- paste0("simplex[", num_networks, "] w; // Weights for networks")
    #     w_prior <- paste0("w ~ dirichlet(rep_vector(0.5, ", num_networks, "));")
    # } else {
    #     w_param <- ""
    #     w_prior <- ""
    # }
  } else {
    network_term <- ""
  }

  # Process ILVs, differing by model type
  ILVi_vars <- STb_data$ILVi_names[!STb_data$ILVi_names %in% "ILVabsent"]
  ILVs_vars <- if (model_type == "full") STb_data$ILVs_names[!STb_data$ILVs_names %in% "ILVabsent"] else character(0)
  ILVm_vars <- if (model_type == "full") STb_data$ILVm_names[!STb_data$ILVm_names %in% "ILVabsent"] else character(0)

  ILVi_vars_clean <- ILVi_vars
  ILVs_vars_clean <- ILVs_vars
  ILVm_vars_clean <- ILVm_vars

  # placeholders for dynamic components
  num_ILVi <- length(ILVi_vars)
  num_ILVs <- length(ILVs_vars)
  num_ILVm <- length(ILVm_vars)

  combined_ILV_vars <- unique(c(ILVi_vars, ILVs_vars, ILVm_vars))

  # create declaration for data block, check dimensions of ILVs for timevarying or constant
  ILV_declaration <- ""
  if (length(combined_ILV_vars) > 0) {
    ILV_declaration <- paste0(
      sapply(combined_ILV_vars, function(var) {
        if (!is.null(dim(STb_data[[paste0("ILV_", var)]]))) {
          paste0("array[K,T_max,P] real ILV_", var, ";")
        } else {
          paste0("array[P] real ILV_", var, ";")
        }
      }),
      collapse = "\n"
    )
  }
  # deal with varying effects in transformed parameters
  N_veff <- length(veff_ID)
  # start w empty list
  transformed_params <- c()
  gq_transformed_params <- c()

  # if user didn't specify veff_ID for baseline
  if (!is.element("lambda_0", veff_ID)) {
    transformed_params <- append(transformed_params, "real<lower=0> lambda_0 = exp(log_lambda_0_mean);")
  }
  # if user didn't specify veff_ID for s
  if (!is.element("s", veff_ID) & model_type == "full") {
    if (separate_s) {
      transformed_params <- append(transformed_params, "vector<lower=0>[N_networks] s_prime = exp(log_s_prime_mean);")
      gq_transformed_params <- append(gq_transformed_params, "vector<lower=0>[N_networks] s = s_prime ./ lambda_0;")
    } else {
      # static scalar s
      transformed_params <- append(transformed_params, "real<lower=0> s_prime = exp(log_s_prime_mean);")
      if (!is.element("lambda_0", veff_ID)) {
        gq_transformed_params <- append(gq_transformed_params, "real<lower=0> s = s_prime/lambda_0;")
      } else {
        gq_transformed_params <- append(gq_transformed_params, "vector<lower=0>[P] s = s_prime/lambda_0;")
      }
    }
  }
  # if user didn't specify veff_ID for f
  if (!is.element("f", veff_ID) & model_type == "full" & transmission_func == "freqdep_f") {
    transformed_params <- append(transformed_params, "real<lower=0> f = exp(log_f_mean);")
  }
  # if user didn't specify veff_ID for k
  if (!is.element("k", veff_ID) & model_type == "full" & transmission_func == "freqdep_k") {
    transformed_params <- append(transformed_params, "real<lower=-1, upper=1> k_shape = 2 / (1 + exp(-k_raw)) - 1;")
  }

  if (!is.element("gamma", veff_ID) & intrinsic_rate == "weibull") {
    transformed_params <- append(transformed_params, "real<lower=0> gamma = exp(log_gamma);")
  }

  count <- 1

  if (N_veff > 0) {
    if ("s" %in% veff_ID & model_type == "full"){
        # this is still a bit weird, adding single varying effect for all networks..
        if (separate_s) {
          # s[n, id] as exp(log_s_prime_mean[n, id] + v_ID[,i])
          transformed_params <- append(transformed_params, glue::glue(
            "matrix<lower=0>[N_networks, P] s_prime;
for (n in 1:N_networks) {{
  for (id in 1:P) {{
    s_prime[n, id] = exp(log_s_prime_mean[n] + v_ID[id, n]);
  }}
}}"
          ))
          gq_transformed_params <- append(gq_transformed_params, glue::glue(
            "matrix<lower=0>[N_networks, P] s_id;
for (n in 1:N_networks) {{
  for (id in 1:P) {{
    s_id[n, id] = s_prime[n, id] / lambda_0[id];
  }}
}}"
          ))
          gq_transformed_params <- append(gq_transformed_params, "vector[N_networks] s_mean = exp(log_s_prime_mean) / exp(log_lambda_0_mean);")
          count <- count + num_networks
        } else {
          # s[id] = exp(log_s_prime_mean + v_ID[,i])
          transformed_params <- append(transformed_params, paste0("vector<lower=0>[P] s_prime = exp(log_s_prime_mean + v_ID[,", count, "]);"))
          gq_transformed_params <- append(gq_transformed_params, paste0("vector<lower=0>[P] s_id = s_prime ./ lambda_0;"))
          gq_transformed_params <- append(gq_transformed_params, paste0("real sprime_mean = exp(log_s_prime_mean);"))
          gq_transformed_params <- append(gq_transformed_params, paste0("real<lower=0> s_mean = (exp(log_s_prime_mean)) / (exp(log_lambda_0_mean));"))
          count <- count + 1
        }
    }
    for (parameter in veff_ID) {
      if (parameter =="s"){
        next
      } else if (parameter == "lambda_0") {
        transformed_params <- append(transformed_params, paste0("vector<lower=0>[P] lambda_0 = exp(log_lambda_0_mean + v_ID[,", count, "]);"))
        gq_transformed_params <- append(gq_transformed_params, paste0("real lambda_0_mean = exp(log_lambda_0_mean);"))
        count <- count + 1
      } else if (parameter == "f" & model_type == "full") {
        transformed_params <- append(transformed_params, paste0("vector<lower=0>[P] f = exp(log_f_mean + v_ID[,", count, "]);"))
        transformed_params <- append(transformed_params, paste0("real<lower=0> f_mean = exp(log_f_mean);"))
        count <- count + 1
      } else if (parameter == "k" & model_type == "full") {
        transformed_params <- append(transformed_params, paste0("vector<lower=0>[P] k_shape = 2 / (1 + exp(-k_raw+ v_ID[,", count, "])) - 1;"))
        transformed_params <- append(transformed_params, paste0("real<lower=-1, upper=1> k_shape_mean = 2 / (1 + exp(-k_raw)) - 1;"))
        count <- count + 1
      } else if (parameter == "gamma") {
        transformed_params <- append(transformed_params, paste0("vector<lower=0>[P] gamma = exp(log_gamma + v_ID[,", count, "]);"))
        transformed_params <- append(transformed_params, paste0("real<lower=0> mean_gamma = exp(log_gamma);"))
        count <- count + 1
      }
    }
  }

  # Handle asocial ILV (ILVi)
  ilvi_result <- process_ILVs(ILVi_vars, ILVi_vars_clean, veff_ID, "i", STb_data, count, prior_beta)
  ILVi_param <- ilvi_result$param
  ILVi_prior <- ilvi_result$prior
  ILVi_variable_effects <- ilvi_result$term
  transformed_params <- append(transformed_params, ilvi_result$transformed)
  count <- ilvi_result$count


  # Handle social ILV (ILVs)
  ilvs_result <- process_ILVs(ILVs_vars, ILVs_vars_clean, veff_ID, "s", STb_data, count, prior_beta)
  ILVs_param <- ilvs_result$param
  ILVs_prior <- ilvs_result$prior
  ILVs_variable_effects <- ilvs_result$term
  transformed_params <- append(transformed_params, ilvs_result$transformed)
  count <- ilvs_result$count


  # Handle multiplicative ILV
  ilvm_result <- process_ILVs(ILVm_vars, ILVm_vars_clean, veff_ID, "m", STb_data, count, prior_beta)
  ILVm_param <- ilvm_result$param
  ILVm_prior <- ilvm_result$prior
  ILVm_variable_effects <- ilvm_result$term
  transformed_params <- append(transformed_params, ilvm_result$transformed)
  count <- ilvm_result$count

  # collapse lists into multiline statements
  transformed_params_declaration <- paste0(transformed_params, collapse = "\n")
  gq_transformed_params_declaration <- paste0(gq_transformed_params, collapse = "\n")

  N_veff_priors <- if (N_veff > 0) glue::glue("
    to_vector(z_ID) ~ {prior_z_ID};
    sigma_ID ~ {prior_sigma_ID};
    Rho_ID ~ {prior_rho_ID};") else ""

  functions_block <- if (transmission_func == "freqdep_k") glue::glue("
functions {{
  real dini_func(real x, real k) {{
    // transform x from [0,1] to [-1,1]
    real x_transformed = 2 * x - 1;
    real y = ((x_transformed - k * x_transformed) / (k - 2 * k * abs(x_transformed) + 1) + 1) / 2;
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
    {if (high_res_k) 'array[N_networks, K, T_max, P] real<lower=0, upper=1> prop_k;' else ''}
}}
")
  # Parameters block
  parameters_block <- glue::glue("
parameters {{
    {distribution_param_declaration}
    real log_lambda_0_mean;  // Log baseline learning rate
    {s_param}
    {if (model_type=='full') {f_param} else ''}
    {if (model_type=='full') {k_param} else ''}
    {gamma_param}
    {ILVi_param}
    {if (model_type=='full') {ILVs_param} else ''}
    {if (model_type=='full') {ILVm_param} else ''}
    {if (N_veff > 0) '
    matrix[N_veff,P] z_ID;
    vector<lower=0, upper=3>[N_veff] sigma_ID;
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
  # create string inputs cuz recursion don't work 2 levels down in glue
  gamma_statement <- if (intrinsic_rate == "weibull") glue::glue("* pow(elapsed_time, {gamma_term} - 1)") else ""
  eA_gamma_statement <- if (intrinsic_rate == "weibull") glue::glue("* pow(global_time, {gamma_term} - 1)") else ""


  if (model_type == "full") {
    psoc_code <- if (model_type=="full") {
      get_ST_prob_term(
        transmission_func = transmission_func,
        is_distribution = is_distribution,
        separate_s = separate_s,
        veff_ID = veff_ID,
        num_networks = num_networks,
        ILVs_variable_effects = ILVs_variable_effects,
        weibull_term = gamma_statement,
        high_res=high_res_k)
    }

    social_info_statement <- glue::glue(
      "real soc_term = net_effect{ILVs_variable_effects};"
    )

    lambda_statement <- glue::glue(
      "real lambda = {ILVm_variable_effects} ({if (is.element('lambda_0', veff_ID)) 'lambda_0[id]' else 'lambda_0'} * ind_term + soc_term) * D[trial, time_step] {gamma_statement};"
    )

    lambda_statement_estAcq <- glue::glue(
      "real lambda = {ILVm_variable_effects} ({if (is.element('lambda_0', veff_ID)) 'lambda_0[id]' else 'lambda_0'} * ind_term + soc_term) {eA_gamma_statement};"
    )

    if (dTADA) {
      target_increment_statement <- glue::glue(
        "target += lambda + log1m_exp(-lambda);"
      )

      log_lik_statement <- glue::glue(
        "// dTADA: probability of learning within interval\n" %+%
          "log_lik_matrix[trial, n] = lambda + log1m_exp(-lambda) - cum_hazard;"
      )
    } else {
      target_increment_statement <- glue::glue(
        "target += log({ILVm_variable_effects} ({if (is.element('lambda_0', veff_ID)) 'lambda_0[id]' else 'lambda_0'} * ind_term + soc_term){gamma_statement});"
      )

      log_lik_statement <- glue::glue(
        "log_lik_matrix[trial, n] = log({ILVm_variable_effects} ({if (is.element('lambda_0', veff_ID)) 'lambda_0[id]' else 'lambda_0'} * ind_term + soc_term){gamma_statement}) - cum_hazard;"
      )
    }
  } else if (model_type == "asocial") {
    psoc_code <- ""
    social_info_statement <- ""
    lambda_statement <- glue::glue("real lambda =  {if (is.element('lambda_0', veff_ID)) 'lambda_0[id]' else 'lambda_0'} * ind_term * D[trial, time_step]{gamma_statement};")
    lambda_statement_estAcq <- glue::glue("real lambda =  {if (is.element('lambda_0', veff_ID)) 'lambda_0[id]' else 'lambda_0'} * ind_term{eA_gamma_statement};")
    target_increment_statement <- glue::glue("target += log( {if (is.element('lambda_0', veff_ID)) 'lambda_0[id]' else 'lambda_0'} * ind_term{gamma_statement});")
    log_lik_statement <- glue::glue("log_lik_matrix[trial, n] = log({if (is.element('lambda_0', veff_ID)) 'lambda_0[id]' else 'lambda_0'} * ind_term{gamma_statement}) - cum_hazard;")
  }

  # Model block
  model_block <- glue::glue("
model {{
    log_lambda_0_mean ~ {prior_lambda0};
    {if (model_type=='full') paste0('log_s_prime_mean ~ ',prior_s,';') else ''}
    {if (model_type=='full') {f_prior} else ''}
    {if (model_type=='full') {k_prior} else ''}
    {gamma_prior}
    {ILVi_prior}
    {if (model_type=='full') {ILVs_prior} else ''}
    {if (model_type=='full') {ILVm_prior} else ''}
    {distribution_model_block}

    {N_veff_priors}

    for (trial in 1:K) {{


        for (n in 1:N[trial]) {{
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];

            if (learn_time > 0) {{
                for (time_step in 1:learn_time) {{
                    {if (intrinsic_rate=='weibull') 'real elapsed_time = sum(D[trial, 1:time_step]);' else ''}
                    real ind_term = {ILVi_variable_effects};
                    {network_term}
                    {social_info_statement}
                    {lambda_statement}
                    target += -lambda;
                    if (time_step == learn_time) {{
                        {target_increment_statement}
                    }}
                }}
            }}
        }}

        if (N_c[trial] > 0) {{
            for (c in 1:N_c[trial]) {{
                int id = ind_id[trial, N[trial] + c];
                    for (time_step in 1:T[trial]) {{
                        {if (intrinsic_rate=='weibull') 'real elapsed_time = sum(D[trial, 1:time_step]);' else ''}
                        real ind_term = {ILVi_variable_effects};
                        {network_term}
                        {social_info_statement}
                        {lambda_statement}
                        target += -lambda;
                    }}
            }}
        }}
    }}
}}
")
  est_acqTime_code <- if (est_acqTime == TRUE) glue::glue("
    matrix[K, Q] acquisition_time;         // simulated acquisition times
    for (trial in 1:K) {{
        for (n in 1:Q) {{ //have to loop through bc stan
            acquisition_time[trial, n] = time_max[trial];
        }}
        for (n in 1:Q) {{
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            // if demonstrator, skip simulation
            if (learn_time < 0) {{
                acquisition_time[trial, n] = 0;
                continue;
            }}

            real cum_hazard = 0; //set val before adding
            int global_time = 1;
            for (time_step in 1:T[trial]) {{
                for (micro_time in 1:D_int[trial, time_step]){{
                    real ind_term = {ILVi_variable_effects};
                    {network_term}
                    {social_info_statement}
                    {lambda_statement_estAcq}
                    real prob = 1-exp(-lambda);
                    if (bernoulli_rng(prob) && acquisition_time[trial, n]>=time_max[trial]) {{
                        acquisition_time[trial, n] = global_time;
                    }}
                    global_time += 1;
                }}
             }}
        }}
    }}") else ""

  generated_quantities_block <- glue::glue("
generated quantities {{
                                             {gq_transformed_params_declaration}
    {if (N_veff > 0) 'corr_matrix[N_veff] Rho;
    Rho = multiply_lower_tri_self_transpose(Rho_ID);' else ''}
    matrix[K, Q] log_lik_matrix = rep_matrix(0.0, K, Q);           // LL for each observation

    //for %ST
    int count_ST = 0;
    vector[N_networks] psocn_sum = rep_vector(0.0, N_networks);

    for (trial in 1:K) {{
        for (n in 1:N[trial]) {{
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];

            if (learn_time > 0){{
                real cum_hazard = 0; //set val before adding
                for (time_step in 1:learn_time) {{
                    {if (intrinsic_rate=='weibull') 'real elapsed_time = sum(D[trial, 1:time_step]);' else ''}
                    real ind_term = {ILVi_variable_effects};
                    {network_term}
                    {social_info_statement}
                    {lambda_statement}
                    cum_hazard += lambda; // accumulate hazard
                    //if it learn_time, record the ll
                    if (time_step == learn_time){{
                                             {log_lik_statement}
                                             {psoc_code}
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
                        {if (intrinsic_rate=='weibull') 'real elapsed_time = sum(D[trial, 1:time_step]);' else ''}
                        real ind_term = {ILVi_variable_effects};
                        {network_term}
                        {social_info_statement}
                        {lambda_statement}
                        cum_hazard += lambda; // accumulate hazard
                    }}
                // Compute per-individual log likelihood
                log_lik_matrix[trial, N[trial] + c] = -cum_hazard;
            }}
        }}

    }}

    vector[N_networks] percent_ST = psocn_sum / count_ST;

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
  stan_model <- glue::glue("//stan
                             {functions_block}
                             {data_block}
                             {parameters_block}
                             {transformed_parameters_block}
                             {model_block}
                             {if (gq==T) {generated_quantities_block} else ''}")
  stan_model <- gsub("(?m)^[ \\t]*\\n", "", stan_model, perl = TRUE)
  stan_model <- paste0(stan_model, "\n")
  return(stan_model)
}
