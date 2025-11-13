#' generate_STb_model_OADA()
#'
#' Helper function that creates stan code for OADA type models.
#'
#' @param STb_data a list of formatted data returned from the STbayes_data() function
#' @param model_type string specifying the model type: "asocial" or "full"
#' @param transmission_func string specifying transmission function: "standard", "freqdep_f" or "freqdep_k" for frequency dependent complex contagion. Defaults to "standard".
#' @param veff_params Vector of parameter names (string) for which to estimate varying effects. For the moment, this is mutually exclusive for varying effects on trial or ID, you must choose one or the other using argument `veff_type`. Default is no varying effects.
#' @param veff_type string specifying whether varying effects should be applied to "ID" or "trial". Default is ID.
#' @param gq Boolean to indicate whether the generated quantities block is added (incl. ll for WAIC)
#' @param priors named list with strings containing the prior for log s or f. defaults to list(log_s = "uniform(-10, 10)", log_f = "normal(0,1)")

#' @return A STAN model (character) that is customized to the input data.
generate_STb_model_OADA <- function(STb_data,
                                    model_type = "full",
                                    transmission_func = "standard",
                                    veff_params = c(),
                                    veff_type = "id",
                                    gq = TRUE,
                                    priors = list()) {
    # maybe one day this can be an argument :(
    est_acqTime <- FALSE

    if (!model_type %in% c("asocial", "full")) {
        stop("Invalid model_type. Choose 'asocial' or 'full'.")
    }

    # set var for index string to be used for veff
    veff_idx <- paste0(veff_type, collapse = ",")

    # set var for declaring vectors holding estimates, either pop size or num trials
    veff_decl <- list(id = "P", trial = "K")

    # set vars holding priors
    prior_baserate <- priors[["log_lambda0"]]
    prior_s <- priors[["log_sprime"]]
    prior_f <- priors[["log_f"]]
    prior_k <- priors[["k_raw"]]
    prior_z_veff <- priors[["z_veff"]]
    prior_sigma_veff <- priors[["sigma_veff"]]
    prior_rho_veff <- priors[["rho_veff"]]
    prior_beta <- priors[["beta_ILV"]]

    # check if edgeweights are sampled from posterior distribution (import func should have created S variable)
    if ("N_dyad" %in% names(STb_data)) is_distribution <- TRUE else is_distribution <- FALSE
    if ("N_dyad" %in% names(STb_data)) est_acqTime <- FALSE # don't want to deal with that rn

    network_names <- STb_data$network_names
    num_networks <- length(network_names)
    separate_s <- num_networks > 1
    network_key <- if (separate_s) "multi_network" else "single_network"

    # make custom declarations for distributions:
    if (is_distribution && model_type == "full") {
        # data declaration
        distribution_data_declaration <- glue::glue("    matrix[N_networks, N_dyad] logit_edge_mu;  // logit edge values
                      array[N_networks] matrix[N_dyad, N_dyad] logit_edge_cov;  // covariance matrix
                                                    array[N_dyad] int<lower=1> focal_ID;
                                                    array[N_dyad] int<lower=1> other_ID;")
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
                          int edge_idx = 1;
                          for (edge_idx in 1:N_dyad) {{
                                real w = inv_logit(edge_logit[network, edge_idx]);
                                A[network, focal_ID[edge_idx], other_ID[edge_idx]] = w;
                          }}
                        }}")
            } else {
                glue::glue("{matrix_decls}
                        for (network in 1:N_networks){{
                       for (edge_idx in 1:N_dyad) {{
                                real w = inv_logit(edge_logit[network, edge_idx]);
                                A[network, focal_ID[edge_idx], other_ID[edge_idx]] = w;
                                A[network, other_ID[edge_idx], focal_ID[edge_idx]] = w;
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
        f_param <- if (separate_s) "vector[N_networks] log_f_mean;" else "real log_f_mean;"
        f_prior <- paste0("log_f_mean ~ ", prior_f, ";")
    } else {
        f_param <- ""
        f_prior <- ""
    }

    # check if user wants to fit k parameter w complex contagion
    if (transmission_func == "freqdep_k") {
        k_param <- if (separate_s) "vector[N_networks] k_raw;" else "real k_raw;"
        k_prior <- paste0("k_raw ~ ", prior_k, ";")
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
            veff_params = veff_params,
            veff_type = veff_type,
            id_var = "id",
            veff_idx = veff_idx,
            high_res = F
        )

        network_term_j <- get_network_term(
            transmission_func = transmission_func,
            is_distribution = is_distribution,
            num_networks = num_networks,
            veff_params = veff_params,
            veff_type = veff_type,
            id_var = "j",
            veff_idx = veff_idx,
            high_res = F
        )
    } else {
        network_term <- ""
        network_term_j <- ""
    }

    veff_is_id_trial <- setequal(veff_type, c("id", "trial"))

    veff_templates <- list(
        default = list(
            id    = "v_id[,{count}]",
            trial = "v_trial[,{count}]"
        ),
        indexed = list(
            id    = "v_id[id,{count}]",
            trial = "v_trial[trial,{count}]"
        ),
        multinetwork_default = list(
            id    = "v_id[,n]",
            trial = "v_trial[,n]"
        ),
        multinetwork_indexed = list(
            id    = "v_id[id,n]",
            trial = "v_trial[trial,n]"
        ),
        multinetwork_fk = list(
            id    = "v_id[,n+{count}-1]",
            trial = "v_trial[trial,n+{count}-1]"
        )
    )

    # do we have both id and veff
    veff_is_id_trial <- identical(veff_type, c("trial", "id"))

    if (veff_is_id_trial) {
        v_term <- paste(veff_templates$default$id, veff_templates$indexed$trial, sep = " + ")
        v_term_multinetwork <- paste(veff_templates$multinetwork_default$id, veff_templates$multinetwork_indexed$trial, sep = " + ")
        v_term_multinetwork_fk <- paste(veff_templates$multinetwork_fk$id, veff_templates$multinetwork_fk$trial, sep = " + ")
    } else {
        v_term <- paste(veff_templates$default[[veff_type]])
        v_term_multinetwork <- paste(veff_templates$multinetwork_default[[veff_type]])
        v_term_multinetwork_fk <- paste(veff_templates$multinetwork_fk[[veff_type]])
    }

    # deal with varying effects in transformed parameters
    N_veff <- length(veff_params)

    count <- 1

    # SOCIAL TRANSMISSION RATE (S)
    parameter_template_s <- list(
        # default is if there are no veff, no multinetwork
        veff_0 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "real s_prime;",
                    calc = "s_prime = exp(log_s_prime_mean);"
                ),
                gq_block = ""
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "vector<lower=0>[N_networks] s_prime;",
                    calc = "s_prime = exp(log_s_prime_mean);"
                ),
                gq_block = ""
            )
        ),
        veff_1 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "vector<lower=0>[{veff_decl[[veff_type]]}] s_prime;",
                    calc = "s_prime = exp(log_s_prime_mean + {v_term});"
                ),
                gq_block = "real s_prime_mean = exp(log_s_prime_mean);"
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "array[N_networks] vector<lower=0>[{veff_decl[[veff_type]]}] s_prime;",
                    calc = "for (n in 1:N_networks) s_prime[n] = exp(log_s_prime_mean[n] + {v_term_multinetwork});"
                ),
                gq_block = "vector<lower=0>[N_networks] s_prime_mean = exp(log_s_prime_mean);"
            )
        ),
        veff_2 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "array[K] vector<lower=0>[P] s_prime;",
                    calc = "for (trial in 1:K) s_prime[trial] = exp(log_s_prime_mean + {v_term});"
                ),
                gq_block = "real s_prime_mean = exp(log_s_prime_mean);"
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "array[N_networks, K] vector<lower=0>[P] s_prime;",
                    calc = "for (n in 1:N_networks) for (trial in 1:K) s_prime[n, trial] = exp(log_s_prime_mean[n] + {v_term_multinetwork});"
                ),
                gq_block = "vector<lower=0>[N_networks] s_prime_mean = exp(log_s_prime_mean);"
            )
        )
    )

    veff_case <- function(parameter, veff_params, veff_type) {
        is_veff <- parameter %in% veff_params
        both_it <- setequal(veff_type, c("id", "trial"))
        if (!is_veff) {
            return("veff_0")
        }
        if (both_it) {
            return("veff_2")
        }
        return("veff_1")
    }

    if (model_type == "full") {
        veff_key <- veff_case("s", veff_params, veff_type)
        strings <- parameter_template_s[[veff_key]][[network_key]]
        s_decl <- glue::glue(strings$transformed_params_block$decl)
        tmp <- glue::glue(strings$transformed_params_block$calc, .envir = environment())
        s_calc <- glue::glue(tmp, .envir = environment())
        s_gq <- glue::glue(strings$gq_block)
        if (veff_key != "veff_0" && !separate_s) count <- count + 1 else if (veff_key != "veff_0" && separate_s) count <- count + num_networks
    } else {
        s_decl <- s_calc <- s_gq <- ""
    }

    parameter_template_f <- list(
        # default is if there are no veff, no multinetwork
        veff_0 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "real<lower=0> f;",
                    calc = "f = exp(log_f_mean);"
                ),
                gq_block = ""
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "vector<lower=0>[N_networks] f;",
                    calc = "f = exp(log_f_mean);"
                ),
                gq_block = ""
            )
        ),
        veff_1 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "vector<lower=0>[{veff_decl[[veff_type]]}] f;",
                    calc = "f = exp(log_f_mean + {v_term});"
                ),
                gq_block = "real<lower=0> f_mean = exp(log_f_mean);"
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "array[N_networks] vector<lower=0>[{veff_decl[[veff_type]]}] f;",
                    calc = "for (n in 1:N_networks) f[n] = exp(log_f_mean + {v_term_multinetwork_fk});"
                ),
                gq_block = "vector<lower=0>[N_networks] f_mean = exp(log_f_mean);"
            )
        ),
        veff_2 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "array[K] vector<lower=0>[P] f;",
                    calc = "for (trial in 1:K) f[trial] = exp(log_f_mean + {v_term});"
                ),
                gq_block = "real<lower=0> f_mean = exp(log_f_mean);"
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "array[N_networks, K] vector<lower=0>[P] f;",
                    calc = "for (n in 1:N_networks) for (trial in 1:K) f[n, trial] = exp(log_f_mean + {v_term_multinetwork_fk});"
                ),
                gq_block = "vector<lower=0>[N_networks] f_mean = exp(log_f_mean);"
            )
        )
    )

    if (model_type == "full" && transmission_func == "freqdep_f") {
        veff_key <- veff_case("f", veff_params, veff_type)
        strings <- parameter_template_f[[veff_key]][[network_key]]
        f_decl <- glue::glue(strings$transformed_params_block$decl)
        tmp <- glue::glue(strings$transformed_params_block$calc, .envir = environment())
        f_calc <- glue::glue(tmp, .envir = environment())
        f_gq <- glue::glue(strings$gq_block)
        if (veff_key != "veff_0" && !separate_s) count <- count + 1 else if (veff_key != "veff_0" && separate_s) count <- count + num_networks
    } else {
        f_decl <- f_calc <- f_gq <- ""
    }

    parameter_template_k <- list(
        # default is if there are no veff, no multinetwork
        veff_0 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "real<lower=-1, upper=1> k_shape;",
                    calc = "k_shape = 2 / (1 + exp(-k_raw)) - 1;"
                ),
                gq_block = ""
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "vector<lower=0>[N_networks] k_shape;",
                    calc = "k_shape = 2 / (1 + exp(-k_raw)) - 1;"
                ),
                gq_block = ""
            )
        ),
        veff_1 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "vector<lower=0>[{veff_decl[[veff_type]]}] k_shape;",
                    calc = "k_shape = 2 / (1 + exp(-k_raw + {v_term})) - 1;"
                ),
                gq_block = "real<lower=-1, upper=1> k_shape_mean = 2 / (1 + exp(-k_raw)) - 1;"
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "array[N_networks] vector<lower=0>[{veff_decl[[veff_type]]}] k_shape;",
                    calc = "for (n in 1:N_networks) k_shape[n] = 2 / (1 + exp(-k_raw + {v_term_multinetwork_fk})) - 1;"
                ),
                gq_block = "vector<lower=0>[N_networks] k_shape_mean = 2 / (1 + exp(-k_raw)) - 1;"
            )
        ),
        veff_2 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "array[K] vector<lower=0>[P] k_shape;",
                    calc = "for (trial in 1:K) k_shape[trial] = 2 / (1 + exp(-k_raw + {v_term})) - 1;"
                ),
                gq_block = "real<lower=0> k_shape_mean = 2 / (1 + exp(-k_raw)) - 1;"
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "array[N_networks, K] vector<lower=0>[P] k_shape;",
                    calc = "for (n in 1:N_networks) for (trial in 1:K) k_shape[n, trial] = 2 / (1 + exp(-k_raw + {v_term_multinetwork_fk})) - 1;"
                ),
                gq_block = "vector<lower=0>[N_networks] k_shape_mean = 2 / (1 + exp(-k_raw)) - 1;"
            )
        )
    )

    if (model_type == "full" && transmission_func == "freqdep_k") {
        veff_key <- veff_case("k", veff_params, veff_type)
        strings <- parameter_template_k[[veff_key]][[network_key]]
        k_decl <- glue::glue(strings$transformed_params_block$decl)
        tmp <- glue::glue(strings$transformed_params_block$calc, .envir = environment())
        k_calc <- glue::glue(tmp, .envir = environment())
        k_gq <- glue::glue(strings$gq_block)
        if (veff_key != "veff_0" && !separate_s) count <- count + 1 else if (veff_key != "veff_0" && separate_s) count <- count + num_networks
    } else {
        k_decl <- k_calc <- k_gq <- ""
    }



    #### Process ILVs ####
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

    ilv_datatypes <- STb_data$ILV_datatypes
    ilv_n_levels <- STb_data$ILV_n_levels
    ilv_timevarying <- STb_data$ILV_timevarying

    # create declaration for data block, check dimensions of ILVs for timevarying or constant
    ILV_declaration <- ""
    if (length(combined_ILV_vars) > 0) {
        ILV_declaration <- paste0(
            sapply(combined_ILV_vars, function(var) {
                if (!ilv_timevarying[[paste0("ILV_", var)]]) {
                    if (ilv_datatypes[[paste0("ILV_", var)]] == "continuous") {
                        paste0("vector[P] ILV_", var, ";")
                    } else {
                        paste0("matrix[P,", ilv_n_levels[[paste0("ILV_", var)]] - 1, "] ILV_", var, ";")
                    }
                } else {
                    if (ilv_datatypes[[paste0("ILV_", var)]] == "continuous") {
                        paste0("array[K, T_max] vector[P] ILV_", var, ";")
                    } else {
                        paste0("array[K,T_max] matrix[P, ", ilv_n_levels[[paste0("ILV_", var)]] - 1, "] ILV_", var, ";")
                    }
                }
            }),
            collapse = "\n"
        )
    }

    # Handle asocial ILV (ILVi)
    ilvi_result <- process_ILVs(ILVi_vars, ILVi_vars_clean, ilv_datatypes, ilv_n_levels, ilv_timevarying, veff_params, veff_type, "i", STb_data, count, prior_beta)
    ILVi_param <- ilvi_result$param
    ILVi_prior <- ilvi_result$prior
    ILVi_variable_effects <- ilvi_result$term
    count <- ilvi_result$count


    # Handle social ILV (ILVs)
    ilvs_result <- process_ILVs(ILVs_vars, ILVs_vars_clean, ilv_datatypes, ilv_n_levels, ilv_timevarying, veff_params, veff_type, "s", STb_data, count, prior_beta)
    ILVs_param <- ilvs_result$param
    ILVs_prior <- ilvs_result$prior
    ILVs_variable_effects <- ilvs_result$term
    count <- ilvs_result$count

    # Handle multiplicative ILV
    ilvm_result <- process_ILVs(ILVm_vars, ILVm_vars_clean, ilv_datatypes, ilv_n_levels, ilv_timevarying, veff_params, veff_type, "m", STb_data, count, prior_beta)
    ILVm_param <- ilvm_result$param
    ILVm_prior <- ilvm_result$prior
    ILVm_variable_effects <- ilvm_result$term
    count <- ilvm_result$count

    # PUT TOGETHER TRANSFORMED PARAMS and GQ block
    transformed_params <- c(
        s_decl, f_decl, k_decl,
        ilvi_result$transformed_decl,
        ilvs_result$transformed_decl, ilvm_result$transformed_decl
    )

    # append calculations
    transformed_params <- c(
        transformed_params, s_calc, f_calc,
        k_calc, ilvi_result$transformed_calc,
        ilvs_result$transformed_calc, ilvm_result$transformed_calc
    )

    # CALCULATIONS IN GQ
    gq_transformed_params <- c(s_gq, f_gq, k_gq)

    # collapse lists into multiline statements
    transformed_params_declaration <- optimize_transformed_params(transformed_params)
    gq_transformed_params_declaration <- paste0(gq_transformed_params, collapse = "\n")

    # for oada we'll have to index j for part of the likelihood
    # don't gsub "id" in case id is part of a variable name
    ILVi_variable_effects_j <- gsub("[id]", "[j]", ILVi_variable_effects, fixed = TRUE)
    ILVi_variable_effects_j <- gsub("[trial,id]", "[trial,j]", ILVi_variable_effects, fixed = TRUE)
    ILVi_variable_effects_j <- gsub("[trial,time_step,id]", "[trial,time_step,j]", ILVi_variable_effects_j, fixed = TRUE)

    ILVs_variable_effects_j <- gsub("[id]", "[j]", ILVs_variable_effects, fixed = TRUE)
    ILVs_variable_effects_j <- gsub("[trial,id]", "[trial,j]", ILVs_variable_effects, fixed = TRUE)
    ILVs_variable_effects_j <- gsub("[trial,time_step,id]", "[trial,time_step,j]", ILVs_variable_effects_j, fixed = TRUE)

    ILVm_variable_effects_j <- gsub("[id]", "[j]", ILVm_variable_effects, fixed = TRUE)
    ILVm_variable_effects_j <- gsub("[trial,id]", "[trial,j]", ILVm_variable_effects, fixed = TRUE)
    ILVm_variable_effects_j <- gsub("[trial,time_step,id]", "[trial,time_step,j]", ILVm_variable_effects_j, fixed = TRUE)

    # collapse lists into multiline statements
    ILVi_param <- paste0(ILVi_param, collapse = "\n")
    ILVi_prior <- paste0(ILVi_prior, collapse = "\n")
    ILVs_param <- paste0(ILVs_param, collapse = "\n")
    ILVs_prior <- paste0(ILVs_prior, collapse = "\n")
    ILVm_param <- paste0(ILVm_param, collapse = "\n")
    ILVm_prior <- paste0(ILVm_prior, collapse = "\n")


    id_veff_priors <- if (N_veff > 0 && is.element("id", veff_type)) glue::glue("
    to_vector(z_id) ~ {prior_z_veff};
    sigma_id ~ {prior_sigma_veff};
    rho_id ~ {prior_rho_veff};") else ""

    trial_veff_priors <- if (N_veff > 0 && is.element("trial", veff_type)) glue::glue("
    to_vector(z_trial) ~ {prior_z_veff};
    sigma_trial ~ {prior_sigma_veff};
    rho_trial ~ {prior_rho_veff};") else ""

    veff_priors <- glue::glue_collapse(
        c(id_veff_priors, trial_veff_priors),
        sep = "\n"
    )

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

}}
")
    # Parameters block
    id_veff_block <- if (N_veff > 0 && is.element("id", veff_type)) {
        glue::glue("
    matrix[N_veff, {veff_decl$id}] z_id;
    vector<lower=0, upper=3>[N_veff] sigma_id;
    cholesky_factor_corr[N_veff] rho_id;
  ")
    } else {
        ""
    }

    trial_veff_block <- if (N_veff > 0 && is.element("trial", veff_type)) {
        glue::glue("
    matrix[N_veff, {veff_decl$trial}] z_trial;
    vector<lower=0, upper=3>[N_veff] sigma_trial;
    cholesky_factor_corr[N_veff] rho_trial;
  ")
    } else {
        ""
    }

    veff_block <- glue::glue_collapse(
        c(id_veff_block, trial_veff_block),
        sep = "\n"
    )

    parameters_block <- glue::glue("
parameters {{
    {distribution_param_declaration}
    {s_param}
    {k_param}
    {f_param}
    {ILVi_param}
    {if (model_type=='full') {ILVs_param} else ''}
    {if (model_type=='full') {ILVm_param} else ''}
    {veff_block}
}}
")

    # Transformed parameters block
    id_veff_transformed_block <- if (N_veff > 0 && is.element("id", veff_type)) {
        glue::glue("matrix[{veff_decl$id},N_veff] v_id;
    v_id = (diag_pre_multiply(sigma_id, rho_id) * z_id)';")
    } else {
        ""
    }
    trial_veff_transformed_block <- if (N_veff > 0 && is.element("trial", veff_type)) {
        glue::glue("matrix[{veff_decl$trial},N_veff] v_trial;
    v_trial = (diag_pre_multiply(sigma_trial, rho_trial) * z_trial)';")
    } else {
        ""
    }
    veff_transformed_block <- glue::glue_collapse(
        c(id_veff_transformed_block, trial_veff_transformed_block),
        sep = "\n"
    )

    transformed_parameters_block <- glue::glue("
transformed parameters {{
   {veff_transformed_block}
   {transformed_params_declaration}
   {distribution_transformed_declaration}
}}
")

    # create string inputs cuz recursion don't work 2 levels down in glue
    if (model_type == "full") {
        i_social_info_statement <- glue::glue("real i_soc = 1.0 * (net_effect{ILVs_variable_effects});")
        # i_lambda_statement <- glue::glue("real i_lambda = {ILVm_variable_effects} * (i_ind + i_soc);")
        i_lambda_statement <- glue::glue("real i_lambda = {ILVm_variable_effects} (i_ind + i_soc);")
        j_social_info_statement <- glue::glue("real j_soc = 1.0 * (net_effect_j{ILVs_variable_effects_j});")
        # j_lambda_statement <- glue::glue("real j_lambda = {ILVm_variable_effects_j} * (j_ind + j_soc);")
        j_lambda_statement <- glue::glue("real j_lambda = {ILVm_variable_effects_j} (j_ind + j_soc);")

        target_increment_statement <- glue::glue("target += log(i_lambda) - log(sum(j_rates));")
        log_lik_statement <- glue::glue("log_lik_matrix[trial, n] = log(i_lambda) - log(sum(j_rates));")
    } else if (model_type == "asocial") {
        i_social_info_statement <- ""
        # i_lambda_statement <- glue::glue("real i_lambda = {ILVm_variable_effects} * i_ind;")
        i_lambda_statement <- glue::glue("real i_lambda = {ILVm_variable_effects} i_ind;")
        j_social_info_statement <- ""
        # j_lambda_statement <- glue::glue("real j_lambda = {ILVm_variable_effects_j} * j_ind;")
        j_lambda_statement <- glue::glue("real j_lambda = {ILVm_variable_effects_j} j_ind;")
        target_increment_statement <- glue::glue("target += log(i_lambda) - log(sum(j_rates));")
        log_lik_statement <- glue::glue("log_lik_matrix[trial, n] = log(i_lambda) - log(sum(j_rates));")
    }

    # Model block
    model_block <- glue::glue("
model {{
    {if (model_type=='full') paste0('log_s_prime_mean ~ ',prior_s,';') else ''}
    {if (model_type=='full') {f_prior} else ''}
    {if (model_type=='full') {k_prior} else ''}
    {ILVi_prior}
    {if (model_type=='full') {ILVs_prior} else ''}
    {if (model_type=='full') {ILVm_prior} else ''}
    {distribution_model_block}

    {veff_priors}

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

    # I'll have to come back to this, if it's even possible. as it is, the rng can pick the same ind twice, and Z won't match the simulated order
    # if this ever gets uncommented, be sure to add back {est_acqTime_code} into the generated quantities glue.
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
    stan_model <- format_stancode(
        gsub("(?m)^[ \\t]*\\n", "", paste0(stan_model, "\n"), perl = TRUE)
    )
    return(stan_model)
}
