#' STb_lfo()
#'
#' Estimate leave-future-out expected log predictive density for an STbayes model.
#'
#' The function refits the model at selected points and uses Pareto-smoothed
#' importance sampling between refits to approximate M-step-ahead predictive
#' performance. Cross-validation can be done either within trials in parallel or
#' through the full event sequence.
#'
#' @param lfo_type String specifying how future observations are left out.
#' `"parallel"` predicts the next observations within each trial.
#' `"sequential"` moves through trials sequentially.
#' @param fit A fitted `STbayes` model object to use as the initial reference
#'   fit.
#' @param STb_data A formatted `STbayes` data object returned by [import_user_STb()].
#' @param stan_model Character string containing the Stan model code returned by
#' [generate_STb_model()].
#' @param M Integer of future observations to predict at each step.
#' @param L Integer of observations to condition on before the
#'   first prediction/refit. This should be larger than the number of demonstrators.
#' @param k_thres Numeric. Pareto-\eqn{k} threshold above which the model is
#'   refit. Defaults to `0.7`.
#' @param ... Additional arguments passed to the model refitting step.
#'
#' @return An object with attributes:
#' \describe{
#'   \item{`pointwise`}{A data frame containing pointwise `elpd_lfo` and
#'   `pareto_k` values.}
#'   \item{`estimates`}{A named numeric vector containing summaries `elpd_lfo`,
#'   `se_elpd_lfo`, `lfoic`, and `se_lfoic`.}
#'   \item{`refits`}{Observation indices at which the model was refit.}
#' }
#'
#' @export
STb_lfo <- function(lfo_type = c("parallel", "sequential"), fit, STb_data, stan_model,
                    M, L, k_thres = 0.7,
                    ...) {
    event_data <- attr(STb_data, "df_event_data")
    event_data <- subset(event_data, select = -c(id_numeric, index, tie, seed, censored))
    network_data <- attr(STb_data, "df_networks")
    network_is_distribution <- FALSE

    if (is.list(network_data)) {
        network_is_distribution <- TRUE
        is_external_fit <- vapply(
            network_data,
            function(x) inherits(x, "bison_model") || inherits(x, "STRAND Results Object"),
            TRUE
        )
        if (any(is_external_fit)) {
            stop("LFO does not currently work when networks are bisonr or STRAND fits. Please input networks as user-defined arrays when importing data.")
        }
    }

    if (!network_is_distribution) network_data <- subset(network_data, select = -c(focal_numeric, other_numeric))
    network_type_data <- if (STb_data$directed) "directed" else "undirected"

    # choose yer function
    lfo_type <- match.arg(lfo_type)
    subset_events <- get(paste0("subset_events_", lfo_type), mode = "function")
    subset_networks <- get(paste0("subset_networks_", lfo_type), mode = "function")
    subset_tweights <- get(paste0("subset_t_weights_", lfo_type), mode = "function")
    subset_ILV_tv <- get(paste0("subset_ILV_tv_", lfo_type), mode = "function")
    grab_indices <- get(paste0("grab_indices_", lfo_type), mode = "function")

    ILV_c_data <- attr(STb_data, "df_ILV_c")
    ILV_tv_data <- attr(STb_data, "df_ILV_tv")
    t_weights_data <- attr(STb_data, "df_t_weights")

    ILV_names <- list(STb_data$ILVi_names, STb_data$ILVs_names, STb_data$ILVm_names)
    ILV_names <- lapply(ILV_names, function(x) if (x == "ILVabsent") NULL else x)

    if (lfo_type == "parallel") {
        # we need to shape ELPD output as loglik vector with length N*Q
        N <- STb_data$Q # max number of individuals in each trial
        N_trials <- length(STb_data$N) # number of trials
        length_loglik <- N * N_trials
    } else if (lfo_type == "sequential") {
        # we need to shape ELPD output as loglik vector with length N*Q
        N <- sum(STb_data$N) + sum(STb_data$N_c) # max number of individuals in each trial
        N_trials <- length(STb_data$N) # number of trials
        length_loglik <- N
    }

    # print(paste("N =", N))


    out <- c(1:(length_loglik))
    elpd_output <- rep(NA_real_, length_loglik)
    ks <- rep(NA_real_, length_loglik)
    # conv <- vector("list", max_N)

    refits <- numeric(0)
    refits <- c(refits, L)

    fit_star <- fit
    i_star <- L

    init <- fit_star$draws(format = "list")

    # annoyingly need to get nveff from the stan model
    N_veff <- return_N_veff(stan_model)
    message(paste("Detected N_veff =", N_veff))
    STb_data$N_veff <- N_veff

    # write stan model to temp file and compile once
    stan_file <- tempfile(
        pattern = "STb_model_",
        tmpdir = tempdir(),
        fileext = ".stan"
    )

    writeLines(stan_model, con = stan_file)

    stan_model <- cmdstanr::cmdstan_model(stan_file)

    # begin with i=L, predict L+M
    # print("Initializing with L obs.")

    # first subset past and predicted observations in event_data
    event_data_past <- subset_events(L, event_data)
    event_data_oos <- subset_events(L + M, event_data)
    idx_from_alldata <- grab_indices(L, M, event_data)
    idx_from_oosdata <- grab_indices(L, M, event_data_oos)

    # network data
    network_data_past <- subset_networks(L, network_data, event_data_past)
    network_data_oos <- subset_networks(L + M, network_data, event_data_oos)

    # t_weights
    if (!is.null(t_weights_data)) {
        t_weights_data_past <- subset_t_weights(L, t_weights_data, event_data_past)
        t_weights_data_oos <- subset_t_weights(L + M, t_weights_data, event_data_past)
    } else {
        t_weights_data_past <- t_weights_data_oos <- NULL
    }

    # ILV_tv
    if (!is.null(ILV_tv_data)) {
        ILV_tv_data_past <- subset_t_weights(L, ILV_tv_data, event_data_past)
        ILV_tv_data_oos <- subset_t_weights(L + M, ILV_tv_data, event_data_past)
    } else {
        ILV_tv_data_past <- ILV_tv_data_oos <- NULL
    }

    # then create new STb_data
    standata_past <- suppressMessages(import_user_STb(
        event_data = event_data_past,
        networks = network_data_past,
        network_type = network_type_data,
        ILV_c = ILV_c_data,
        ILV_tv = ILV_tv_data_past,
        ILVi = ILV_names[[1]],
        ILVs = ILV_names[[2]],
        ILVm = ILV_names[[3]],
        t_weights = t_weights_data_past,
        high_res = STb_data$high_res
    ))
    standata_oos <- suppressMessages(import_user_STb(
        event_data = event_data_oos,
        networks = network_data_oos,
        network_type = network_type_data,
        ILV_c = ILV_c_data,
        ILV_tv = ILV_tv_data_oos,
        ILVi = ILV_names[[1]],
        ILVs = ILV_names[[2]],
        ILVm = ILV_names[[3]],
        t_weights = t_weights_data_oos,
        high_res = STb_data$high_res
    ))

    standata_past$N_veff <- standata_oos$N_veff <- N_veff

    # refit with past data
    fit_star <- refit_model(standata_past, stan_model, init)

    # predict with oos data
    predict_oos <- predict_model(standata_oos, stan_model, fit_star)

    # calc elpds
    loglik <- log_lik_matrix(predict_oos)
    elpds <- calc_elpds(loglik, oos_idx = idx_from_oosdata, psis = NULL, M = M)
    elpd_output[idx_from_alldata] <- elpds
    # conv[[L]] <- convergence_summary(fit_star)
    ks[idx_from_alldata] <- 0

    # start from L + 1 as we already handled L above
    for (i in (L + 1):(N - M)) {
        # print(paste("Looping through idx", i))
        oos <- (i + 1):(i + M)

        # subset data to 1:max(oos)
        event_data_oos <- subset_events(max(oos), event_data)
        # print(paste("OOS event data length =", nrow(event_data_oos)))

        # need to calculate 2 sets of idx, one with the entire dataset so
        # ELPD gets assigned to correct locations in the end
        idx_from_alldata <- grab_indices(i, M, event_data)

        # second, we need the indexes of OOS obs from the subsetted data to
        # pass to the logratio function
        idx_from_oosdata <- grab_indices(i, M, event_data_oos)
        idx_4_paretok <- grab_indices(i_star, i - i_star, event_data_oos)

        # remake standata with oos
        network_data_oos <- subset_networks(max(oos), network_data, event_data_oos)

        # t weights
        if (!is.null(t_weights_data)) {
            t_weights_data_oos <- subset_t_weights(max(oos), t_weights_data, event_data_past)
        } else {
            t_weights_data_oos <- NULL
        }

        # t weights
        if (!is.null(ILV_tv_data)) {
            ILV_tv_data_oos <- subset_t_weights(max(oos), ILV_tv_data, event_data_past)
        } else {
            ILV_tv_data_oos <- NULL
        }

        standata_oos <- suppressMessages(import_user_STb(
            event_data = event_data_oos,
            networks = network_data_oos,
            network_type = network_type_data,
            ILV_c = ILV_c_data,
            ILV_tv = ILV_tv_data_oos,
            ILVi = ILV_names[[1]],
            ILVs = ILV_names[[2]],
            ILVm = ILV_names[[3]],
            t_weights = t_weights_data_oos,
            high_res = STb_data$high_res
        ))
        standata_oos$N_veff <- N_veff

        # do the actual predictions. This will create a loglik matrix that is the
        # size of 1:OOS
        predict_oos <- predict_model(standata_oos, stan_model, fit_star)

        # extract loglik matrix size S x Obs
        loglik <- log_lik_matrix(predict_oos)
        # print(paste("dimensions of loglik", dim(loglik)))
        # this creates S sums of log ratios, but should only contain values for predicted idx?
        logratio <- sum_log_ratios(loglik, ids = idx_4_paretok)
        # print(paste("Dimensions of logratio", length(logratio)))
        psis_obj <- suppressWarnings(loo::psis(logratio))
        k <- loo::pareto_k_values(psis_obj)
        ks[idx_from_alldata] <- k

        # if k is bad we need to refit before predicting oos
        if (k > k_thres) {
            print(paste("k exceeds threshold. k =", k))
            i_star <- i
            refits <- c(refits, i)

            event_data_past <- subset_events(i, event_data)
            network_data_past <- subset_networks(i, network_data, event_data_past)

            # t_weights
            if (!is.null(t_weights_data)) {
                t_weights_data_past <- subset_t_weights(i, t_weights_data, event_data_past)
            } else {
                t_weights_data_past <- NULL
            }

            # ILV_tv
            if (!is.null(ILV_tv_data)) {
                ILV_tv_data_past <- subset_t_weights(i, ILV_tv_data, event_data_past)
            } else {
                ILV_tv_data_past <- NULL
            }

            standata_past <- suppressMessages(import_user_STb(
                event_data = event_data_past,
                networks = network_data_past,
                network_type = network_type_data,
                ILV_c = ILV_c_data,
                ILV_tv = ILV_tv_data_past,
                ILVi = ILV_names[[1]],
                ILVs = ILV_names[[2]],
                ILVm = ILV_names[[3]],
                t_weights = t_weights_data_past,
                high_res = STb_data$high_res
            ))
            standata_past$N_veff <- N_veff
            fit_star <- refit_model(standata_past, stan_model, init)
            predict_oos <- predict_model(standata_oos, stan_model, fit_star)
            loglik <- log_lik_matrix(predict_oos)
            # conv[[i]] <- convergence_summary(fit_star)
        }
        elpds <- calc_elpds(loglik,
            oos_idx = idx_from_oosdata,
            psis = if (k > k_thres) NULL else psis_obj,
            M = M
        )
        # print(paste("Writing elpd", elpds, "to ouput index", idx_from_alldata))
        elpd_output[idx_from_alldata] <- elpds
    }

    n_pointwise <- sum(!is.na(elpd_output))

    elpd_lfo <- sum(elpd_output, na.rm = TRUE)
    se_elpd_lfo <- sqrt(n_pointwise) * stats::sd(elpd_output, na.rm = TRUE)

    attr(out, "pointwise") <- data.frame(
        elpd_lfo = elpd_output,
        pareto_k = ks
    )

    attr(out, "estimates") <- c(
        elpd_lfo = elpd_lfo,
        se_elpd_lfo = se_elpd_lfo,
        lfoic = -2 * elpd_lfo,
        se_lfoic = 2 * se_elpd_lfo
    )

    attr(out, "refits") <- refits

    # give em a little sugar
    print(
        paste0(
            "Using a pareto-k threshold of ", k_thres,
            ", model was refit ", length(attr(out, "refits")),
            " times at observations ",
            paste(attr(out, "refits"), collapse = ",")
        )
    )

    return(out)
}

#' STb_compare_lfo()
#'
#' Compare two or more objects returned by [STb_lfo()] using leave-future-out
#' expected log predictive density. The model with the highest total `elpd_lfo`
#' is treated as the reference model. Differences are reported relative to that model.
#'
#' @param ... Two or more LFO objects returned by [STb_lfo()].
#' @param model_names Optional character vector of model names. If `NULL`, names
#'   are taken from `...` when available, otherwise generic names are used.
#' @param n_obs Optional integer giving the number of pointwise observations to
#'   use when computing standard errors. If `NULL`, the number of complete
#'   pointwise observations across models is used.
#' @param digits Integer. Number of digits used for rounding the returned table.
#'   Defaults to `1`.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{`elpd_diff`}{Difference in total ELPD from the best model.}
#'   \item{`se_diff`}{Standard error of the ELPD difference.}
#'   \item{`elpd_lfo`}{Total leave-future-out ELPD.}
#'   \item{`se_elpd_lfo`}{Standard error of total leave-future-out ELPD.}
#'   \item{`lfoic`}{LFO information criterion, computed as `-2 * elpd_lfo`.}
#'   \item{`se_lfoic`}{Standard error of `lfoic`.}
#' }
#'
#' @export
STb_compare_lfo <- function(..., model_names = NULL, n_obs = NULL, digits = 1) {
    models <- list(...)

    if (length(models) < 2) {
        stop("Please provide at least two lfo objects.")
    }

    dots_names <- names(models)

    if (is.null(model_names)) {
        if (!is.null(dots_names) && all(nzchar(dots_names))) {
            model_names <- dots_names
        } else {
            model_names <- paste0("model", seq_along(models))
        }
    }

    if (length(model_names) != length(models)) {
        stop("model_names must have the same length as the number of models")
    }

    get_pointwise_elpd <- function(x, model_name) {
        pointwise <- attr(x, "pointwise")

        if (is.null(pointwise)) {
            stop("model '", model_name, "' is missing attr(x, 'pointwise')")
        }

        if (!is.data.frame(pointwise)) {
            stop("attr(x, 'pointwise') must be a data.frame for model '", model_name, "'")
        }

        if (!"elpd_lfo" %in% names(pointwise)) {
            stop("attr(x, 'pointwise') must contain column 'elpd_lfo' for model '", model_name, "'")
        }

        as.numeric(pointwise$elpd_lfo)
    }

    elpd_list <- Map(get_pointwise_elpd, models, model_names)

    lens <- vapply(elpd_list, length, integer(1))
    if (length(unique(lens)) != 1) {
        stop("all pointwise elpd_lfo vectors must have the same length")
    }

    n_pointwise <- if (is.null(n_obs)) {
        sum(stats::complete.cases(do.call(cbind, elpd_list)))
    } else {
        n_obs
    }

    summarize_one <- function(elpd) {
        elpd_lfo <- sum(elpd, na.rm = TRUE)
        se_elpd_lfo <- sqrt(n_pointwise) * stats::sd(elpd, na.rm = TRUE)

        c(
            elpd_lfo = elpd_lfo,
            se_elpd_lfo = se_elpd_lfo,
            lfoic = -2 * elpd_lfo,
            se_lfoic = 2 * se_elpd_lfo
        )
    }

    estimates <- do.call(rbind, lapply(elpd_list, summarize_one))
    estimates <- as.data.frame(estimates)

    estimates$model <- model_names

    best_idx <- which.max(estimates$elpd_lfo)
    best_elpd <- elpd_list[[best_idx]]

    diff_list <- lapply(elpd_list, function(elpd) {
        elpd - best_elpd
    })

    estimates$elpd_diff <- vapply(diff_list, sum, numeric(1), na.rm = TRUE)

    estimates$se_diff <- vapply(diff_list, function(diff) {
        sqrt(sum(!is.na(diff))) * stats::sd(diff, na.rm = TRUE)
    }, numeric(1))

    estimates <- estimates[order(estimates$elpd_lfo, decreasing = TRUE), ]

    estimates <- estimates[, c(
        "model",
        "elpd_diff",
        "se_diff",
        "elpd_lfo",
        "se_elpd_lfo",
        "lfoic",
        "se_lfoic"
    )]

    rownames(estimates) <- estimates$model
    estimates$model <- NULL

    estimates[] <- lapply(estimates, function(x) round(x, digits))

    class(estimates) <- c("lfo_compare", class(estimates))
    estimates
}

#' plot_ks()
#'
#' Create diagnostic plot of Pareto-\eqn{k} values from LFO-CV.
#'
#' @param ks Numeric vector of Pareto-\eqn{k} values from
#'   `attr(x, "pointwise")$pareto_k` for an object returned by [STb_lfo()].
#' @param L Integer. Initial number of observations used before LFO prediction
#'   begins.
#' @param k_thres Numeric. Pareto-\eqn{k} threshold to draw as a reference line.
#'   Defaults to `0.7`.
#'
#' @return NULL
#'
#' @export
plot_ks <- function(ks, L, k_thres = 0.7) {
    ids <- (L + 1):(length(ks) - 1)
    dat_ks <- data.frame(
        ks = ks[ids],
        ids = ids
    )

    minval <- min(ks, na.rm = T) - .5
    maxval <- max(2, max(ks, na.rm = T) + .5)

    cols <- ifelse(dat_ks$ks > k_thres, "darkblue", "cornflowerblue")

    plot(
        dat_ks$ids,
        dat_ks$ks,
        pch = 3,
        col = cols,
        xlab = "Data point",
        ylab = "Pareto k",
        ylim = c(minval, maxval)
    )

    abline(h = k_thres, lty = 2, col = "red2")
}


#' @noRd
log_sum_exp <- function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
}

#' @noRd
log_mean_exp <- function(x) {
    log_sum_exp(x) - log(length(x))
}

#' Compute log of raw importance ratios
#' @noRd
sum_log_ratios <- function(ll, ids = NULL) {
    # print(ids)
    if (!is.null(ids)) ll <- ll[, ids, drop = FALSE]
    out <- rowSums(ll)
    out
}

#' Extract loglik from cmdstanr fit as matrix
#' @noRd
log_lik_matrix <- function(fit) {
    fit$draws("log_lik", format = "matrix")
}

#' Calculate elpds loglik from loglik matrix for specified indices
#' @noRd
calc_elpds <- function(loglik, oos_idx = NULL, psis = NULL, M = 1, ...) {
    loglik_subset <- loglik[, oos_idx]

    K <- length(oos_idx) / M

    sum_ll <- matrix(NA, nrow = nrow(loglik_subset), ncol = K)

    # for each block of M obs in a trial, we need to sum the samples
    for (k in seq_len(K)) {
        cols <- ((k - 1) * M + 1):(k * M)
        sum_ll[, k] <- rowSums(loglik_subset[, cols, drop = FALSE])
    }

    if (is.null(psis)) {
        # which one is right here? since i could be predicting multiple M, do they all get the same ELPD?
        # or should export by observation number?..
        # out <- log_mean_exp(sum_ll)
        out <- apply(FUN = log_mean_exp, MARGIN = 2, sum_ll)
    } else {
        lw <- weights(psis, normalize = TRUE)[, 1]
        sum_ll <- lw + sum_ll
        # same here, not sure if i export 1 value or one for each oos_idx
        # out <- log_sum_exp(lw + sum_ll)
        out <- apply(FUN = log_sum_exp, MARGIN = 2, lw + sum_ll)
    }

    return(out)
}

#' Calculates which indices should be retained given M-steps ahead
#' @noRd
grab_indices_parallel <- function(i, M, event_data) {
    trial_indices <- split(seq_along(event_data$trial), event_data$trial)
    keep_indices <- unlist(lapply(trial_indices, function(idx) {
        if (length(idx) < i + M) {
            integer(0) # drop trials with fewer than i + M obs
        } else {
            idx[(i + 1):(i + M)]
        }
    }))
    return(keep_indices)
}

#' Calculates which indices should be retained given M-steps ahead
#' @noRd
grab_indices_sequential <- function(i, M, event_data) {
    keep_indices <- (i + 1):(i + M)
    # print(paste("Predicting:", keep_indices))
    return(keep_indices)
}

#' Calculates which indices should be retained given M-steps ahead
#' @noRd
subset_events_parallel <- function(i, data) {
    # split by trial
    split_df <- split(data, data$trial_numeric)

    # process each trial
    out <- lapply(split_df, function(d) {
        # order by time
        ord <- order(d$time)
        d_sorted <- d[ord, ]

        # if i exceeds number of rows
        if (i > nrow(d_sorted)) {
            return(d)
        }

        # otherwise get cutoff time (ith event)
        t_cut <- d_sorted$time[i]

        # update t_end
        d$t_end <- t_cut

        # cap times above cutoff
        d$time[d$time > t_cut] <- t_cut + 1

        return(d)
    })

    # recombine
    do.call(rbind, out)
}

#' @noRd
subset_events_sequential <- function(i, data) {
    # get current trial
    current_trial <- data$trial_numeric[i]

    # only take up through this trial
    d_sub <- data[data$trial_numeric <= current_trial, ]

    if (nrow(d_sub) > i) {
        t_cut <- d_sub$time[i]
        idx <- d_sub$trial_numeric == current_trial
        d_sub$t_end[idx] <- t_cut
        d_sub$time[idx & d_sub$time > t_cut] <- t_cut + 1
    }
    return(d_sub)
}

#' @noRd
elpd_diffs <- function(loo_a, loo_b) {
    pt_a <- loo_a$pointwise
    pt_b <- loo_b$pointwise
    elpd <- grep("^elpd", colnames(pt_a))
    pt_b[, elpd] - pt_a[, elpd]
}

#' @noRd
subset_networks_parallel <- function(i, network_data, event_data) {
    current_time <- event_data$discrete_time[i]

    if (inherits(network_data, "data.frame")) {
        if ("discrete_time" %in% names(network_data)) {
            # need to add 1 here for the interval between last learned and end of obs. period
            network_data <- network_data[network_data$discrete_time <= current_time + 1, ]
            network_data <- network_data[
                order(network_data$trial_numeric, network_data$discrete_time),
            ]
        } else {
            network_data <- network_data[order(network_data$trial_numeric), ]
        }

        return(network_data)
    }

    if (inherits(network_data, "array")) {
        network_data <- list(network_data)
    }

    if (!is.list(network_data)) {
        stop("network_data must be a data.frame, array, or list of arrays.")
    }

    lapply(network_data, function(net) {
        if (!inherits(net, "array")) {
            stop("all list elements in network_data must be arrays.")
        }

        dn <- names(dimnames(net))

        if (is.null(dn) || any(dn == "")) {
            stop("network arrays must have named dimensions.")
        }

        valid_dims <- c("draw", "trial", "time", "focal_ID", "other_ID")
        unknown_dims <- setdiff(dn, valid_dims)

        if (length(unknown_dims) > 0) {
            stop(
                "unknown network array dimension(s): ",
                paste(unknown_dims, collapse = ", ")
            )
        }

        required_dims <- c("draw", "focal_ID", "other_ID")
        missing_dims <- setdiff(required_dims, dn)

        if (length(missing_dims) > 0) {
            stop(
                "network array is missing required dimension(s): ",
                paste(missing_dims, collapse = ", ")
            )
        }

        idx <- setNames(vector("list", length(dn)), dn)
        for (d in dn) idx[[d]] <- TRUE

        if ("time" %in% dn) {
            idx[["time"]] <- idx[["time"]] <- seq_len(max(current_time + 1, 1))
        }

        do.call(`[`, c(list(net), unname(idx), list(drop = FALSE)))
    })
}

#' @noRd
subset_t_weights_parallel <- function(i, t_weights_data, event_data) {
    current_time <- event_data$discrete_time[i]

    if ("discrete_time" %in% names(t_weights_data)) {
        t_weights_data <- t_weights_data[t_weights_data$discrete_time <= current_time, ]
        t_weights_data <- t_weights_data[order(t_weights_data$trial_numeric, t_weights_data$discrete_time), ]
    } else {
        t_weights_data <- t_weights_data[order(t_weights_data$trial_numeric), ]
    }
    return(t_weights_data)
}

#' @noRd
subset_ILV_tv_parallel <- function(i, ILV_tv_data, event_data) {
    current_time <- event_data$discrete_time[i]

    if ("discrete_time" %in% names(ILV_tv_data)) {
        ILV_tv_data <- ILV_tv_data[ILV_tv_data$discrete_time <= current_time, ]
        ILV_tv_data <- ILV_tv_data[order(ILV_tv_data$trial_numeric, ILV_tv_data$discrete_time), ]
    } else {
        ILV_tv_data <- ILV_tv_data[order(ILV_tv_data$trial_numeric), ]
    }
    return(ILV_tv_data)
}

#' @noRd
subset_networks_sequential <- function(i, network_data, event_data) {
    current_trial <- event_data$trial_numeric[i]
    current_time <- event_data$discrete_time[i]

    # old dataframe behavior
    if (inherits(network_data, "data.frame")) {
        d_sub <- network_data[network_data$trial_numeric <= current_trial, ]

        idx_curr <- d_sub$trial_numeric == current_trial
        d_sub <- d_sub[!(idx_curr & d_sub$discrete_time > (current_time + 1)), ]

        return(d_sub)
    }

    # allow a bare array, but return a list for array inputs
    if (inherits(network_data, "array")) {
        network_data <- list(network_data)
    }

    if (!is.list(network_data)) {
        stop("network_data must be a data.frame, array, or list of arrays.")
    }

    lapply(network_data, function(net) {
        if (!inherits(net, "array")) {
            stop("all list elements in network_data must be arrays.")
        }

        dn <- names(dimnames(net))

        if (is.null(dn) || any(dn == "")) {
            stop("network arrays must have named dimensions.")
        }

        valid_dims <- c("draw", "trial", "time", "focal_ID", "other_ID")
        unknown_dims <- setdiff(dn, valid_dims)

        if (length(unknown_dims) > 0) {
            stop(
                "unknown network array dimension(s): ",
                paste(unknown_dims, collapse = ", ")
            )
        }

        required_dims <- c("draw", "focal_ID", "other_ID")
        missing_dims <- setdiff(required_dims, dn)

        if (length(missing_dims) > 0) {
            stop(
                "network array is missing required dimension(s): ",
                paste(missing_dims, collapse = ", ")
            )
        }

        # no trial dimension means there is nothing sequential to subset
        if (!"trial" %in% dn) {
            return(net)
        }

        keep_trials <- seq_len(current_trial)

        # if there is no time dimension, just keep trials up to current_trial
        if (!"time" %in% dn) {
            idx <- setNames(vector("list", length(dn)), dn)
            for (d in dn) idx[[d]] <- TRUE
            idx[["trial"]] <- keep_trials

            return(do.call(`[`, c(list(net), unname(idx), list(drop = FALSE))))
        }

        # for trial-time arrays, later trials are dropped and the current
        # trial is censored after current_time. because arrays need rectangular
        # dimensions, later times in earlier trials are retained.
        idx <- setNames(vector("list", length(dn)), dn)
        for (d in dn) idx[[d]] <- TRUE
        idx[["trial"]] <- keep_trials

        out <- do.call(`[`, c(list(net), unname(idx), list(drop = FALSE)))

        # mask future network states in the current trial
        out_dn <- names(dimnames(out))
        time_dim <- dim(out)[match("time", out_dn)]


        first_future_time <- current_time + 2

        if (first_future_time <= time_dim) {
            idx_future <- setNames(vector("list", length(out_dn)), out_dn)
            for (d in out_dn) idx_future[[d]] <- TRUE

            idx_future[["trial"]] <- length(keep_trials)
            idx_future[["time"]] <- seq.int(first_future_time, time_dim)

            out <- do.call(`[<-`, c(
                list(out),
                unname(idx_future),
                list(value = 0)
            ))
        }
        return(out)
    })
}

#' @noRd
subset_t_weights_sequential <- function(i, t_weights_data, event_data) {
    # identify current trial from global index i
    current_trial <- event_data$trial_numeric[i]
    # print(paste("Subsetting t_weights up to trial", current_trial))

    # keep only up to current trial
    d_sub <- t_weights_data[t_weights_data$trial_numeric <= current_trial, ]

    # work on current trial only
    idx_curr <- d_sub$trial_numeric == current_trial
    d_curr <- d_sub[idx_curr, ]

    # unique sorted times within trial
    current_time <- event_data$discrete_time[i]
    # print(paste("current_time:", current_time))
    # remove rows with times above current time
    d_sub <- d_sub[!(idx_curr & d_sub$discrete_time > current_time), ]

    return(d_sub)
}

#' @noRd
subset_ILV_tv_sequential <- function(i, ILV_tv_data, event_data) {
    # identify current trial from global index i
    current_trial <- event_data$trial_numeric[i]
    print(paste("Subsetting ILV_tv_data up to trial", current_trial))

    # keep only up to current trial
    d_sub <- ILV_tv_data[ILV_tv_data$trial_numeric <= current_trial, ]

    # work on current trial only
    idx_curr <- d_sub$trial_numeric == current_trial
    d_curr <- d_sub[idx_curr, ]

    # unique sorted times within trial
    current_time <- event_data$discrete_time[i]
    print(paste("current_time:", current_time))
    # remove rows with times above current time
    d_sub <- d_sub[!(idx_curr & d_sub$discrete_time > current_time), ]

    return(d_sub)
}

#' @noRd
refit_model <- function(STb_data, stan_model, init = NULL, ...) {
    print("Refitting model.")
    data_list_clean <- Filter(function(x) {
        is.numeric(x) || is.integer(x) || is.logical(x) || is.array(x) || is.matrix(x)
    }, STb_data)

    fit_new <- suppressMessages(suppressWarnings(stan_model$sample(
        data = data_list_clean,
        init = init,
        chains = 4,
        parallel_chains = 4,
        iter_warmup = 1000,
        iter_sampling = 1000,
        refresh = 0,
        show_messages = FALSE,
        output_dir = tempdir(),
        ...
    )))
    return(fit_new)
}

#' @noRd
predict_model <- function(new_data, stan_model, fitted_params, ...) {
    # print("Predicting OOS.")

    data_list_clean <- Filter(function(x) {
        is.numeric(x) || is.integer(x) || is.logical(x) || is.array(x) || is.matrix(x)
    }, new_data)

    predictions <- suppressMessages(suppressWarnings(stan_model$generate_quantities(
        data = data_list_clean,
        fitted_params = fitted_params,
        parallel_chains = 4
    )))
    return(predictions)
}
