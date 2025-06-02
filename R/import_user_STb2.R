#' Import user data 2
#'
#' First step of analysis pipeline: create STbayes data object from user
#' supplied data to be used for generating and fitting models. This function is
#' basically used when *fitting models of complex transmission to high-resolution data.*
#' Rather than pre-process high res data, this will just create a massive
#' data-list to be sent to Stan. Models created with this function will take much longer to run.
#'
#' @param event_data dataframe with columns id, trial, time, t_end
#' @param networks Either a dataframe, a bisonr fit, or a list of bisonr fits.
#' If dataframe: with columns trial, from, to, and one or more columns of edge weights named descriptively.
#' @param network_type "undirected" or "directed".
#' @param ILV_c optional dataframe with columns id, and any constant individual-level variables that might be of interest
#' @param ILV_tv optional dataframe with columns trial, id, time and any time-varying variables. Variable values should summarize the variable for each inter-acquisition period.
#' @param ILVi Optional character vector of column names from ILV metadata to be considered when estimating intrinsic rate. If not specified, all ILV are applied to both.
#' @param ILVs Optional character vector of column names from ILV metadata to be considered when estimating social transmission rate. If not specified, all ILV are applied to both.
#' @param ILVm Optional character vector of column names from ILV metadata to be considered in a multiplicative model.
#' @param t_weights Optional dataframe with columns trial, id, time and t_weight. Transmission rates represent rates of production/relevant cues per inter-event period.
#' @param high_res Boolean indicating whether or not user is providing networks and transmission weights per period duration=1
#' @importFrom Rcpp evalCpp
#' @importFrom stats qlogis median cov aggregate
#' @return A list object containing properly formatted data to run social transmission models.
#' @export
import_user_STb2 <- function(event_data,
                             networks,
                             network_type = c("undirected", "directed"),
                             ILV_c = NULL,
                             ILV_tv = NULL,
                             ILVi = NULL,
                             ILVs = NULL,
                             ILVm = NULL,
                             t_weights = NULL,
                             high_res = FALSE) {
    if (inherits(event_data, "data.frame")) {
        check_required_cols(event_data, required_cols = c("id", "trial", "time", "t_end"), df_name = "event_data")
    } else {
        stop("\U0001F614 Please feed me a dataframe for the event_data argument.")
    }

    if (inherits(networks, "data.frame")) {
        check_required_cols(networks, required_cols = c("trial", "from", "to"), df_name = "networks")
    }

    # other warnings
    if (inherits(networks, "data.frame")) {
        message("User supplied edge weights as point estimates \U0001F4CD")
        is_distribution <- FALSE
    } else if (inherits(networks, "bison_model") || inherits(networks, "STRAND Results Object")) {
        message("User supplied edge weights as posterior distributions \U0001F308")
        networks <- list(networks)
        is_distribution <- TRUE
    } else if (is.list(networks) && all(sapply(networks, function(x) inherits(x, "bison_model") || inherits(x, "STRAND Results Object")))) {
        message("User supplied a list of Bayesian network fits [\U0001F308 , \U0001F308] ")
        is_distribution <- TRUE
    } else {
        stop("\U0001F614 Please feed me a dataframe for the networks argument, or a bisonR model fit, or a list of fits.")
    }

    if (all(is.null(ILVi), is.null(ILVs), is.null(ILVm)) & (!is.null(ILV_c) | !is.null(ILV_tv))) {
        message("\U0001F914 You have provided ILVs, yet did not specify whether they should be additive or multiplicative (missing arguments ILVi, ILVs, ILVm). They will not be included in the model.")
    }

    id_check <- standardize_ids(networks, event_data, ILV_c, ILV_tv, t_weights)
    event_data <- id_check$event_data
    ILV_c <- id_check$ILV_c
    ILV_tv <- id_check$ILV_tv
    t_weights <- id_check$t_weights

    network_type <- match.arg(network_type)

    #### event_data ####
    # event_data should be in format id, trial, time, t_end
    # if time==0, assume to be trained demonstrator, if time>t_end, assume to be censored
    # create numeric variables in case user has supplied strings
    # event_data$id_numeric <- as.numeric(as.factor(event_data$id))
    event_data$trial_numeric <- as.numeric(as.factor(event_data$trial))
    # order in case user has not, assign indexes per trial
    event_data <- event_data[order(event_data$trial_numeric, event_data$time), ]
    event_data$index <- with(event_data, ave(trial_numeric, trial_numeric, FUN = seq_along))

    if (high_res) {
        event_data$discrete_time <- event_data$time
    } else {
        # create discrete time (this should be 0 if ID was a demo/seed)
        event_data$discrete_time <- NA
        event_data$discrete_time[event_data$time != 0] <- with(
            event_data[event_data$time != 0, ],
            ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x)))
        )
        event_data$discrete_time[event_data$time == 0] <- 0
    }

    # not using these now but just in case..
    event_data$tie <- with(event_data, ave(time, interaction(trial_numeric, time), FUN = function(x) length(x) > 1))
    event_data$seed <- ifelse(event_data$time == 0, 1, 0)

    # summarize censored/uncensored for each trial
    event_data$censored <- ifelse(event_data$time <= event_data$t_end, 0, 1)
    N_data <- do.call(rbind, lapply(split(event_data, event_data$trial_numeric), function(df) {
        data.frame(
            trial_numeric = unique(df$trial_numeric),
            num_uncensored = nrow(df) - sum(df$censored),
            num_censored = sum(df$censored),
            max_periods = if (high_res) max(df$t_end) else max(df$discrete_time)
        )
    }))

    # create matrix of where rows = trial and columns = discrete_time, and values = duration
    # summarize by trial and time
    temp_data <- event_data
    temp_data$time[temp_data$time > temp_data$t_end] <- temp_data$t_end[temp_data$time > temp_data$t_end]
    D_data <- unique(temp_data[order(temp_data$trial_numeric, temp_data$time), c("trial_numeric", "time", "discrete_time")])
    D_data <- D_data[D_data$time != 0, ]

    if (high_res) {
        # Expand to full grid of [trial x time] from 1 to max time observed
        max_times <- aggregate(discrete_time ~ trial_numeric, data = D_data, FUN = max)
        D_data <- do.call(rbind, lapply(1:nrow(max_times), function(i) {
            data.frame(
                trial_numeric = max_times$trial_numeric[i],
                discrete_time = seq(1, max_times$discrete_time[i]),
                duration = 1
            )
        }))
    } else {
        # Calculate the duration as time - lag(time) for each grouping
        D_data$duration <- with(D_data, ave(time, trial_numeric, FUN = function(x) c(x[1], diff(x))))
        D_data$duration <- ifelse(is.na(D_data$duration), D_data$time, D_data$duration)
        D_data <- D_data[, !(names(D_data) %in% "time")]
    }

    # pivot wider
    D_data_real <- with(D_data, tapply(duration, list(trial_numeric, discrete_time), FUN = max, default = 0))
    D_data_int <- with(D_data, tapply(duration, list(trial_numeric, as.integer(discrete_time)), FUN = max, default = 0))
    dim(D_data_int) <- dim(D_data_real)
    dimnames(D_data_real) <- NULL
    dimnames(D_data_int) <- NULL

    # create a matrix where rows are trial_numeric, columns are id_numeric, and values are discrete_time
    t_data <- with(event_data, tapply(discrete_time, list(trial_numeric, id_numeric), FUN = max, default = -1))
    # create a matrix where rows are trial_numeric, columns are id_numeric, and values are actual time measurements
    time_data <- with(event_data, tapply(time, list(trial_numeric, id_numeric), FUN = max, default = -1))
    # gets the end of obs period of each trial
    time_max <- with(event_data, tapply(t_end, list(trial_numeric), FUN = max, default = -1))

    # create a matrix where rows are trial_numeric, columns are N + N_c, and values are id_numeric
    # values are then left aligned such that individuals who were in each trial will appear in the first indices
    # i know this is a really weird way of doing a ragged list, in case there are differing numbers of individuals in each trial... but should work fine
    id_data <- with(event_data, tapply(id_numeric, list(trial_numeric, index), FUN = max, default = NA))
    id_data <- t(apply(id_data, 1, function(row) c(stats::na.omit(row), rep(-1, sum(is.na(row))))))

    data_list <- list(
        K = length(unique(event_data$trial_numeric)),
        P = length(unique(event_data$id_numeric)),
        N = structure(N_data$num_uncensored, dim = length(N_data$num_uncensored)),
        N_c = structure(N_data$num_censored, dim = length(N_data$num_censored)),
        T = structure(N_data$max_periods, dim = length(N_data$max_periods)),
        t = if (high_res) time_data else t_data,
        T_max = max(N_data$max_periods),
        time = time_data,
        time_max = time_max,
        Q = max(event_data$index),
        D = if (high_res) D_data_real[, 1:max(N_data$max_periods)] else D_data_real,
        D_int = if (high_res) D_data_int[, 1:max(N_data$max_periods)] else D_data_int,
        ind_id = id_data,
        Zn = create_Z_matrix(event_data, high_res, if (high_res) "real" else "interval"),
        Z = create_Z_matrix(event_data, high_res, if (high_res) "real" else "interval"),
        high_res = high_res,
        directed = if (network_type == "directed") T else F
    )

    dim(data_list$D) <- c(data_list$K, data_list$T_max)
    dim(data_list$D_int) <- c(data_list$K, data_list$T_max)

    if (!is.null(t_weights)) {
        data_list$W <- create_W_matrix(t_weights, data_list$T_max)
        if (!all(dim(data_list$Z) == dim(data_list$W))) {
            stop(paste0(
                "Dimensions of Z (",
                paste(dim(data_list$Z), collapse = ","),
                ") do not match W (",
                paste(dim(data_list$W), collapse = ","), ")."
            ))
        }
        data_list$Z <- data_list$Z * data_list$W
        # data_list$Z <- data_list$W
    }

    #### Constant ILV ####
    ILV_datatypes <- c()
    ILV_names <- c()
    # identify ILVs from ILV_c
    if (!is.null(ILV_c)) {
        message("Constant ILV supplied.")
        # ILV_c$id_numeric <- as.numeric(as.factor(ILV_c$id))
        # order in case user has not
        # ILV_c <- ILV_c[order(ILV_c$id_numeric), ]
        rownames(ILV_c) <- NULL

        # get column names
        exclude_cols <- c("id", "id_numeric")
        ILV_cols <- setdiff(names(ILV_c), exclude_cols)
        ILV_names <- append(ILV_names, ILV_cols)
        # loop through each ILV_c column and add to datalist
        for (col in ILV_cols) {
            ILV_datatypes[[paste0("ILV_", col)]] <- detect_ILV_datatype(ILV_c$col)
            data_list[[paste0("ILV_", col)]] <- ILV_c[[col]]
        }
    }

    #### Time-varying ILV ####
    # identify ILVs from ILV_tv
    if (!is.null(ILV_tv)) {
        message("Time-varying ILV supplied.")
        # convert to numeric if not
        # ILV_tv$id_numeric <- as.numeric(as.factor(ILV_tv$id))
        ILV_tv$trial_numeric <- as.numeric(as.factor(ILV_tv$trial))
        ILV_tv$discrete_time <- with(
            ILV_tv, ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x)))
        )
        # order in case user has not
        # ILV_tv <- ILV_tv[order(ILV_tv$trial_numeric, ILV_tv$id_numeric, ILV_tv$discrete_time), ]
        rownames(ILV_tv) <- NULL
        # Get the ILV column names
        exclude_cols <- c("id", "id_numeric", "time", "discrete_time", "trial", "trial_numeric")
        ILV_cols <- setdiff(names(ILV_tv), exclude_cols)
        ILV_names <- append(ILV_names, ILV_cols)
        # loop through each ILV column and add to datalist
        for (col in ILV_cols) {
            # get data type
            ILV_datatypes[[paste0("ILV_", col)]] <- detect_ILV_datatype(ILV_tv$col)
            # reshape data into matrix
            mat <- with(ILV_tv, tapply(ILV_tv[[col]], list(trial, discrete_time, id), FUN = mean, simplify = TRUE))
            mat[is.na(mat)] <- 0
            data_list[[paste0("ILV_", col)]] <- mat
        }
    }

    # write names
    if (!is.null(ILV_c) | !is.null(ILV_tv)) {
        if (is.null(ILVi)) {
            data_list$ILVi_names <- "ILVabsent"
        } else {
            data_list$ILVi_names <- ILVi
        }
        if (is.null(ILVs)) {
            data_list$ILVs_names <- "ILVabsent"
        } else {
            data_list$ILVs_names <- ILVs
        }
        if (is.null(ILVm)) {
            data_list$ILVm_names <- "ILVabsent"
        }
        # write datatypes
        data_list$ILV_datatypes <- ILV_datatypes
    } else {
        message("No ILV supplied.")
        data_list$ILVi_names <- "ILVabsent"
        data_list$ILVs_names <- "ILVabsent"
        data_list$ILVm_names <- "ILVabsent"
    }


    if (!is_distribution) {
        if (data_list$P != length(unique(c(networks$from, networks$to)))) {
            stop("Networks do not contain the same number of unique individuals as the event data.")
        }
        if (data_list$K != length(unique(networks$trial))) {
            stop("Networks do not contain the same number of trials as the event data.")
        }
        is_dynamic <- "time" %in% names(networks)
        if (!is_dynamic) {
            message("User input indicates static network(s). If dynamic, include 'time' column.")
            networks$time <- 1
        }
        exclude_cols <- c("trial", "time", "from", "to")
        network_cols <- setdiff(names(networks), exclude_cols)

        temp_names <- data.frame(
            name = sort(unique(c(networks$from, networks$to))),
            numeric = seq_along(sort(unique(c(networks$from, networks$to))))
        )
        networks$from_numeric <- as.integer(temp_names$numeric[match(networks$from, temp_names$name)])
        networks$to_numeric <- as.integer(temp_names$numeric[match(networks$to, temp_names$name)])
        networks$trial_numeric <- as.integer(as.factor(networks$trial))
        if (!high_res) {
            networks$discrete_time <- with(networks, ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x))))
        } else {
            networks$discrete_time <- networks$time
        }
        is_symmetric <- nrow(networks[networks$trial_numeric == 1 & networks$discrete_time == min(networks$discrete_time), ]) == data_list$P * (data_list$P - 1)
        max_timesteps <- if (high_res) max(event_data$t_end) else max(data_list$T)

        dims <- c(length(network_cols), data_list$K, max_timesteps, data_list$P, data_list$P)
        A_array <- array(0, dim = dims)

        # flatten ONCE
        A_flat <- as.numeric(aperm(A_array, c(1, 2, 3, 4, 5)))

        # fill in-place
        for (n in seq_along(network_cols)) {
            column <- network_cols[n]
            for (k in 1:data_list$K) {
                temp_df <- networks[networks$trial_numeric == k, ]
                from <- temp_df$from_numeric
                to <- temp_df$to_numeric
                time <- temp_df$discrete_time
                value <- temp_df[[column]]

                fill_array(
                    A_flat, dims,
                    from, to, time, value,
                    n - 1, k - 1,
                    (!is_symmetric && network_type == "undirected")
                )
            }
        }

        # reshape ONCE after all fill_array calls
        A_array <- aperm(array(A_flat, dim = dims), c(1, 2, 3, 4, 5))

        # replicate across time for static networks
        if (!is_dynamic) {
            for (n in seq_along(network_cols)) {
                for (k in 1:data_list$K) {
                    A_array[n, k, , , ] <- A_array[n, k, rep(1, dim(A_array)[3]), , ]
                }
            }
        }

        # zero out diagonals
        for (n in seq_along(network_cols)) {
            for (k in 1:data_list$K) {
                for (t in 1:max(data_list$T)) {
                    diag(A_array[n, k, t, , ]) <- 0
                }
            }
        }
        if (min(A_array) < 0) stop("Edgeweights below zero detected. Rescale so that 0 = no connection.")
        data_list$A <- A_array
        data_list$network_names <- network_cols
        data_list$N_networks <- length(network_cols)
    }
    # IF DISTRIBUTION
    else {
        N_networks <- length(networks)

        for (i in seq_along(networks)) {
            net_obj <- networks[[i]]

            if (inherits(net_obj, "bison_model")) {
                # BISON handling
                dyads <- if (net_obj$directed) {
                    strsplit(net_obj$dyad_names, " -> ")
                } else {
                    strsplit(net_obj$dyad_names, " <-> ")
                }

                dyad_matrix <- do.call(rbind, dyads)
                dyad_df <- data.frame(
                    from = as.integer(dyad_matrix[, 1]),
                    to = as.integer(dyad_matrix[, 2]),
                    stringsAsFactors = FALSE
                )

                # Drop self-loops
                edge_samples <- net_obj$edge_samples[, dyad_df$from != dyad_df$to, drop = FALSE]
                dyad_df <- dyad_df[dyad_df$from != dyad_df$to, ]
                data_list$from_ID <- dyad_df$from
                data_list$to_ID <- dyad_df$to

                if (net_obj$directed) {
                    # Use samples as-is
                    edges <- edge_samples
                } else {
                    # Duplicate edges to simulate directed version
                    # from_ids <- c(dyad_df$from, dyad_df$to)
                    # to_ids   <- c(dyad_df$to, dyad_df$from)
                    # dyad_df <- data.frame(from = from_ids, to = to_ids)
                    # edges <- cbind(edge_samples, edge_samples)
                    edges <- edge_samples
                }
            } else if (inherits(net_obj, "STRAND Results Object")) {
                # STRAND handling
                ass_matrix <- net_obj$samples$predicted_network_sample # [draws, from, to]
                draws <- dim(ass_matrix)[1]
                P <- dim(ass_matrix)[2]

                # Create edge sample matrix: [draws, dyads]
                edge_mat <- matrix(NA, nrow = draws, ncol = P * (P - 1))
                idx <- 1
                for (i_from in 1:P) {
                    for (i_to in 1:P) {
                        if (i_from != i_to) {
                            edge_mat[, idx] <- ass_matrix[, i_from, i_to]
                            idx <- idx + 1
                        }
                    }
                }

                edges <- qlogis(edge_mat) # convert to logit space
            } else {
                stop("\u274C Unrecognized distributional network object. Expected bison_model.")
            }

            if (i == 1) {
                N_dyad <- ncol(edges)
                logit_edge_mu <- matrix(NA, nrow = N_networks, ncol = N_dyad)
                logit_edge_cov <- array(NA, dim = c(N_networks, N_dyad, N_dyad))
            }

            # Get point estimate and covariance
            logit_edge_mu[i, ] <- apply(edges, 2, median)
            logit_edge_cov[i, , ] <- cov(edges)
        }

        # Assume same number of dyads for all networks
        data_list$logit_edge_mu <- logit_edge_mu
        data_list$logit_edge_cov <- logit_edge_cov
        data_list$N_dyad <- N_dyad
        data_list$network_names <- paste0("net", seq_len(N_networks))
        data_list$N_networks <- N_networks
    }

    ## Final sanity check and return ##
    dl_sanity_check(data_list = data_list)
    return(data_list)
}
