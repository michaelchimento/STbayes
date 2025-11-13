#' generate_STb_model_TADA()
#'
#' Helper function that creates stan code for cTADA and dTADA type models.
#'
#' @param STb_data a list of formatted data returned from the STbayes_data() function
#' @param model_type string specifying the model type: "asocial" or "full"
#' @param intrinsic_rate string specifying whether intrinsic rate is static or can change over time ("constant", "weibull")
#' @param transmission_func string specifying transmission function: "standard", "freqdep_f" or "freqdep_k" for complex contagion. Defaults to "standard". Complex contagion with multi-network model is not supported.
#' @param dTADA boolean indicating whether dTADA should be used.
#' @param veff_params Vector of parameter names (string) for which to estimate varying effects. Default is no varying effects.
#' @param veff_type string/vector specifying whether varying effects should be applied to "id", "trial" or both (c("id","trial")). Default is id, applied only if user also gives `veff_params`.
#' @param gq Boolean to indicate whether the generated quantities block is added (incl. ll for WAIC)
#' @param est_acqTime Boolean to indicate whether gq block includes estimates for acquisition time. At the moment this uses 'one weird trick' to accomplish this and does not support estimates for non-integer learning times.
#' @param priors named list with strings containing the prior for log baserate, s, f, k.
#' @return A STAN model (character) that is customized to the input data.

generate_STb_model_TADA <- function(STb_data,
                                    model_type = "full",
                                    intrinsic_rate = "constant",
                                    transmission_func = "standard",
                                    dTADA = F,
                                    veff_params = c(),
                                    veff_type = "id",
                                    gq = TRUE,
                                    est_acqTime = FALSE,
                                    priors = list()) {
    if (!model_type %in% c("asocial", "full")) {
        stop("Invalid model_type. Choose 'asocial' or 'full'.")
    }

    if (est_acqTime == TRUE && min(check_integer(STb_data$time)) == 0) {
        message("WARNING: You have input float times, and unfortunately estimating acquisition times in the GQ block is only possible with integer times at the moment.\nThe model will be created with est_acqTime=F.")
        est_acqTime <- FALSE
    }

    # set var for index string to be used for veff
    veff_idx <- paste0(veff_type, collapse = ",")

    # set var for declaring vectors holding estimates, either pop size or num trials
    veff_decl <- list(id = "P", trial = "K")

    prior_lambda0 <- priors[["log_lambda0"]]
    prior_s <- priors[["log_sprime"]]
    prior_f <- priors[["log_f"]]
    prior_k <- priors[["k_raw"]]
    prior_z_veff <- priors[["z_veff"]]
    prior_sigma_veff <- priors[["sigma_veff"]]
    prior_rho_veff <- priors[["rho_veff"]]
    prior_beta <- priors[["beta_ILV"]]
    prior_gamma <- priors[["gamma"]]

    # check if edgeweights are sampled from posterior distribution
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

    # check if user wants to fit f parameter w complex contagion
    if (intrinsic_rate == "weibull") {
        gamma_param <- "real log_gamma;"
        gamma_prior <- paste0("log_gamma ~ ", prior_gamma, ";")
        if (is.element("gamma", veff_params)) gamma_term <- paste0("gamma[", veff_idx, "]") else gamma_term <- "gamma"
    } else {
        gamma_param <- ""
        gamma_prior <- ""
    }


    # Declare network variables
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
    } else {
        network_term <- ""
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
                gq_block = if (!is.element("lambda_0", veff_params)) "real<lower=0> s = s_prime/lambda_0;" else "real<lower=0> s = s_prime/exp(log_lambda_0_mean);"
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "vector<lower=0>[N_networks] s_prime;",
                    calc = "s_prime = exp(log_s_prime_mean);"
                ),
                gq_block = if (!is.element("lambda_0", veff_params)) "vector<lower=0>[N_networks] s = s_prime ./ lambda_0;" else "vector<lower=0>[N_networks] s = s_prime ./ exp(log_lambda_0_mean);"
            )
        ),
        veff_1 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "vector<lower=0>[{veff_decl[[veff_type]]}] s_prime;",
                    calc = "s_prime = exp(log_s_prime_mean + {v_term});"
                ),
                gq_block = "real<lower=0> s_mean = exp(log_s_prime_mean) / exp(log_lambda_0_mean);"
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "array[N_networks] vector<lower=0>[{veff_decl[[veff_type]]}] s_prime;",
                    calc = "for (n in 1:N_networks) s_prime[n] = exp(log_s_prime_mean[n] + {v_term_multinetwork});"
                ),
                gq_block = "vector[N_networks] s_mean = exp(log_s_prime_mean) / exp(log_lambda_0_mean);"
            )
        ),
        veff_2 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "array[K] vector<lower=0>[P] s_prime;",
                    calc = "for (trial in 1:K) s_prime[trial] = exp(log_s_prime_mean + {v_term});"
                ),
                gq_block = "real<lower=0> s_mean = exp(log_s_prime_mean) / exp(log_lambda_0_mean);"
            ),
            multi_network = list(
                transformed_params_block = list(
                    decl = "array[N_networks, K] vector<lower=0>[P] s_prime;",
                    calc = "for (n in 1:N_networks) for (trial in 1:K) s_prime[n, trial] = exp(log_s_prime_mean[n] + {v_term_multinetwork});"
                ),
                gq_block = "vector[N_networks] s_mean = exp(log_s_prime_mean) / exp(log_lambda_0_mean);"
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


    # INTRINSIC RATE (LAMBDA_0)
    parameter_template_lambda_0 <- list(
        # default is if there are no veff, no multinetwork
        veff_0 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "real<lower=0> lambda_0;",
                    calc = "lambda_0 = exp(log_lambda_0_mean);"
                ),
                gq_block = ""
            )
        ),
        veff_1 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "vector<lower=0>[{veff_decl[[veff_type]]}] lambda_0;",
                    calc = "lambda_0 = exp(log_lambda_0_mean + {v_term});"
                ),
                gq_block = "real<lower=0> lambda_0_mean = exp(log_lambda_0_mean);"
            )
        ),
        veff_2 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "array[K] vector<lower=0>[P] lambda_0;",
                    calc = "for (trial in 1:K) lambda_0[trial] = exp(log_lambda_0_mean + {v_term});"
                ),
                gq_block = "real<lower=0> lambda_0_mean = exp(log_lambda_0_mean);"
            )
        )
    )

    veff_key <- veff_case("lambda_0", veff_params, veff_type)
    strings <- parameter_template_lambda_0[[veff_key]]$single_network
    lambda_0_decl <- glue::glue(strings$transformed_params_block$decl)
    tmp <- glue::glue(strings$transformed_params_block$calc, .envir = environment())
    lambda_0_calc <- glue::glue(tmp, .envir = environment())
    lambda_0_gq <- glue::glue(strings$gq_block)
    if (veff_key != "veff_0") count <- count + 1


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

    parameter_template_gamma <- list(
        veff_0 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "real<lower=0> gamma;",
                    calc = "gamma = exp(log_gamma);"
                ),
                gq_block = ""
            )
        ),
        veff_1 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "vector<lower=0>[{veff_decl[[veff_type]]}] gamma;",
                    calc = "gamma = exp(log_gamma + {v_term});"
                ),
                gq_block = "real<lower=0> gamma_mean = exp(log_gamma);"
            )
        ),
        veff_2 = list(
            single_network = list(
                transformed_params_block = list(
                    decl = "array[K] vector<lower=0>[P] gamma;",
                    calc = "for (trial in 1:K) gamma[trial] = exp(log_gamma + {v_term});"
                ),
                gq_block = "real<lower=0> gamma_mean = exp(log_gamma);"
            )
        )
    )

    if (intrinsic_rate == "weibull") {
        veff_key <- veff_case("gamma", veff_params, veff_type)
        strings <- parameter_template_gamma[[veff_key]]$single_network
        gamma_decl <- glue::glue(strings$transformed_params_block$decl)
        tmp <- glue::glue(strings$transformed_params_block$calc, .envir = environment())
        gamma_calc <- glue::glue(tmp, .envir = environment())
        gamma_gq <- glue::glue(strings$gq_block)
        if (veff_key != "veff_0") count <- count + 1
    } else {
        gamma_decl <- gamma_calc <- gamma_gq <- ""
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
        s_decl, lambda_0_decl, f_decl, k_decl, gamma_decl,
        ilvi_result$transformed_decl,
        ilvs_result$transformed_decl,
        ilvm_result$transformed_decl
    )

    # append calculations
    transformed_params <- c(
        transformed_params, s_calc, lambda_0_calc, f_calc,
        k_calc, gamma_calc, ilvi_result$transformed_calc,
        ilvs_result$transformed_calc,
        ilvm_result$transformed_calc
    )

    # CALCULATIONS IN GQ
    gq_transformed_params <- c(s_gq, lambda_0_gq, f_gq, k_gq, gamma_gq)

    # collapse lists into multiline statements
    transformed_params_declaration <- optimize_transformed_params(transformed_params)

    gq_transformed_params_declaration <- paste0(gq_transformed_params, collapse = "\n")

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
    real log_lambda_0_mean;  // Log baseline learning rate
    {s_param}
    {if (model_type=='full') {f_param} else ''}
    {if (model_type=='full') {k_param} else ''}
    {gamma_param}
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
    gamma_statement <- if (intrinsic_rate == "weibull") glue::glue("* pow(elapsed_time, {gamma_term} - 1)") else ""
    eA_gamma_statement <- if (intrinsic_rate == "weibull") glue::glue("* pow(global_time, {gamma_term} - 1)") else ""
    lambda_var <- if (is.element("lambda_0", veff_params)) glue::glue("lambda_0[{veff_idx}]") else "lambda_0"

    if (model_type == "full") {
        psoc_code <- if (model_type == "full") {
            get_ST_prob_term(
                transmission_func = transmission_func,
                is_distribution = is_distribution,
                separate_s = separate_s,
                veff_params = veff_params,
                veff_idx = veff_idx,
                num_networks = num_networks,
                ILVs_variable_effects = ILVs_variable_effects,
                weibull_term = gamma_statement,
                high_res = F
            )
        }

        social_info_statement <- glue::glue(
            "real soc_term = net_effect{ILVs_variable_effects};"
        )



        lambda_statement <- glue::glue(
            "real lambda = {ILVm_variable_effects} ({lambda_var} * ind_term + soc_term) * D[trial, time_step] {gamma_statement};"
        )

        lambda_statement_estAcq <- glue::glue(
            "real lambda = {ILVm_variable_effects} ({lambda_var} * ind_term + soc_term) {eA_gamma_statement};"
        )

        if (dTADA) {
            target_increment_statement <- glue::glue(
                "target += lambda + log1m_exp(-lambda);"
            )

            log_lik_statement <- glue::glue(
                "// dTADA: probability of learning within interval\n",
                "log_lik_matrix[trial, n] = lambda + log1m_exp(-lambda) - cum_hazard;"
            )
        } else {
            target_increment_statement <- glue::glue(
                "target += log({ILVm_variable_effects} ({lambda_var} * ind_term + soc_term){gamma_statement});"
            )

            log_lik_statement <- glue::glue(
                "log_lik_matrix[trial, n] = log({ILVm_variable_effects} ({lambda_var} * ind_term + soc_term){gamma_statement}) - cum_hazard;"
            )
        }
    } else if (model_type == "asocial") {
        psoc_code <- ""
        social_info_statement <- ""
        lambda_statement <- glue::glue("real lambda =  {lambda_var} * ind_term * D[trial, time_step]{gamma_statement};")
        lambda_statement_estAcq <- glue::glue("real lambda =  {lambda_var} * ind_term{eA_gamma_statement};")
        target_increment_statement <- glue::glue("target += log( {lambda_var} * ind_term{gamma_statement});")
        log_lik_statement <- glue::glue("log_lik_matrix[trial, n] = log({lambda_var} * ind_term{gamma_statement}) - cum_hazard;")
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

    {veff_priors}

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
        for (n in 1:Q) {{
            int id = ind_id[trial, n];
            int learn_time = t[trial, id];
            // if demonstrator, skip simulation
            if (learn_time < 0) {{
                acquisition_time[trial, n] = 0;
                continue;
            }}

            real cum_hazard = 0; //set val before adding
            real threshold = -log(uniform_rng(0, 1));
            int global_time = 1;
            acquisition_time[trial, n] = time_max[trial];

            for (time_step in 1:T[trial]) {{
                for (micro_time in 1:D_int[trial, time_step]){{
                    real ind_term = {ILVi_variable_effects};
                    {network_term}
                    {social_info_statement}
                    {lambda_statement_estAcq}
                    cum_hazard += lambda;
                    if (cum_hazard > threshold) {{
                        acquisition_time[trial, n] = global_time;
                        break;  // exit inner loop
                    }}
                    global_time += 1;
                }}
                if (cum_hazard > threshold) break;  // exit outer loop
             }}
        }}
    }}") else ""

    # Parameters block
    id_veff_block <- if (N_veff > 0 && is.element("id", veff_type)) {
        glue::glue("
    corr_matrix[N_veff] Rho_id;
    Rho_id = multiply_lower_tri_self_transpose(rho_id);
  ")
    } else {
        ""
    }

    trial_veff_block <- if (N_veff > 0 && is.element("trial", veff_type)) {
        glue::glue("
    corr_matrix[N_veff] Rho_trial;
    Rho_trial = multiply_lower_tri_self_transpose(rho_trial);
  ")
    } else {
        ""
    }

    veff_block <- glue::glue_collapse(
        c(id_veff_block, trial_veff_block),
        sep = "\n"
    )
    generated_quantities_block <- glue::glue("
generated quantities {{
    {gq_transformed_params_declaration}
    {veff_block}
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
    stan_model <- format_stancode(
        gsub("(?m)^[ \\t]*\\n", "", paste0(stan_model, "\n"), perl = TRUE)
    )
    return(stan_model)
}
