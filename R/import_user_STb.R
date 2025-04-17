#' import_user_STb: Create STbayes_data object from user supplied data
#'
#' @param event_data dataframe with columns id, trial, time, t_end
#' @param networks Either a dataframe, a bisonr fit or STRAND fit, or a list of bisonr or STRAND fits. If dataframe: with columns trial, from, to, and one or more columns of edge weights named descriptively. Optionally an integer time column can be provided for dynamic network analysis, although networks must be provided for each time period between transmission events.
#' @param ILV_c optional dataframe with columns id, and any constant individual-level variables that might be of interest
#' @param ILV_tv optional dataframe with columns trial, id, time and any time-varying variables. Variable values should summarize the variable for each inter-acquisition period.
#' @param ILVi Optional character vector of column names from ILV metadata to be considered when estimating intrinsic rate. If not specified, all ILV are applied to both.
#' @param ILVs Optional character vector of column names from ILV metadata to be considered when estimating social transmission rate. If not specified, all ILV are applied to both.
#' @param ILVm Optional character vector of column names from ILV metadata to be considered in a multiplicative model.
#' @param t_weights Optional dataframe with columns trial, id, time and t_weight. Transmission rates represent rates of production/relevant cues per inter-event period.
#' @param high_res Boolean indicating whether or not user is providing networks and transmission weights per period duration=1
#'
#' @return A list object containing properly formatted data to run social transmission models.
#' @export
#'
#' @examples
#' # very mock data
#' event_data <- data.frame(
#'   trial = rep(1:2, each = 3),
#'   id = LETTERS[1:6],
#'   time = c(0, 1, 2, 0, 1, 4),
#'   t_end = c(3, 3, 3, 4, 4, 4)
#' )
#' networks <- data.frame(
#'   trial = rep(1:2, each = 3),
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
#' ILV_tv <- data.frame(
#'   trial = c(rep(1, each = 9), rep(2, each = 9)),
#'   id = c(rep(LETTERS[1:3], each = 3), rep(LETTERS[4:6], each = 3)),
#'   # these times correspond to the inter-acquisition periods
#'   # e.g. 1 is from [t_0 to t_1), 2 is [t_1 to t_2), 3 = [t_2 to t_3 or t_end] if censored inds. present)
#'   time = c(rep(1:3, times = 3), rep(1:3, times = 3)),
#'   # ensure the variable is summarizing these inter-acquisition time periods
#'   dist_from_task = rnorm(18)
#' )
#' t_weights <- data.frame(
#'   trial = c(rep(1, each = 9), rep(2, each = 9)),
#'   id = c(rep(LETTERS[1:3], each = 3), rep(LETTERS[4:6], each = 3)),
#'   time = c(rep(1:3, times = 3), rep(1:3, times = 3)),
#'   t_weight = exp(rnorm(18))
#' )
#' imported_data <- import_user_STb(
#'   event_data = event_data,
#'   networks = networks,
#'   ILV_c = ILV_c,
#'   ILV_tv = ILV_tv,
#'   ILVi = c("age", "dist_from_task"), # Use 'age' and time-varying ILV 'dist_from_task' for asocial learning
#'   ILVs = c("sex"), # Use only 'sex' for social learning
#'   ILVm = c("weight") # Use weight for multiplicative effect on asocial and social learning
#' )
import_user_STb_2 <- function(event_data,
                            networks,
                            network_type = c("undirected", "directed"),
                            ILV_c = NULL,
                            ILV_tv = NULL,
                            ILVi = NULL,
                            ILVs = NULL,
                            ILVm = NULL,
                            t_weights = NULL,
                            high_res = FALSE) {


    # warnings
    if ( all(is.null(ILVi), is.null(ILVs), is.null(ILVm)) & (!is.null(ILV_c) | !is.null(ILV_tv))) {
        message("WARNING: You have provided ILV, yet did not specify whether they should be additive or multiplicative (missing arguments ILVi, ILVs, ILVm). If not specified, they will not be included in the model.")
    }

    if (inherits(networks, "data.frame")) {
        message("User supplied edge weights as point estimates.")
        is_distribution <- FALSE
    } else if (inherits(networks, "bison_model") || inherits(networks, "STRAND Results Object")) {
        message("User supplied a single Bayesian network fit.")
        networks <- list(networks)
        is_distribution <- TRUE
    } else if (is.list(networks) && all(sapply(networks, function(x) inherits(x, "bison_model") || inherits(x, "STRAND Results Object")))) {
        message("User supplied a list of Bayesian network fits (bisonr or STRAND).")
        is_distribution <- TRUE
    } else {
        stop("Please supply a dataframe, a bison_model or STRAND fit object, or a list of them.")
    }

    network_type <- match.arg(network_type)

    #### event_data ####
    # event_data should be in format id, trial, time, t_end
    # if time==0, assume to be trained demonstrator, if time>t_end, assume to be censored
    # create numeric variables in case user has supplied strings
    event_data$id_numeric <- as.numeric(as.factor(event_data$id))
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

    #not using these now but just in case..
    event_data$tie <- with(event_data, ave(time, interaction(trial_numeric, time), FUN = function(x) length(x) > 1))
    event_data$seed <- ifelse(event_data$time == 0, 1, 0)

    # summarize censored/uncensored for each trial
    event_data$censored <- ifelse(event_data$time <= event_data$t_end, 0, 1)
    N_data <- do.call(rbind, lapply(split(event_data, event_data$trial_numeric), function(df) {
        data.frame(
            trial_numeric = unique(df$trial_numeric),
            num_uncensored = nrow(df) - sum(df$censored),
            num_censored = sum(df$censored),
            max_periods = max(df$discrete_time)
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
        full_grid <- do.call(rbind, lapply(1:nrow(max_times), function(i) {
            data.frame(trial_numeric = max_times$trial_numeric[i],
                       discrete_time = seq(1, max_times$discrete_time[i]),
                       duration = 1)
        }))
        D_data <- full_grid
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
        D = D_data_real,
        D_int = D_data_int,
        ind_id = id_data,
        Zn = create_Z_matrix(event_data),
        Z = create_Z_matrix(event_data)
    )

    if (!is.null(t_weights)) {
        data_list$W <- create_W_matrix(t_weights, data_list$T_max)
        if (!all(dim(data_list$Z) == dim(data_list$W))) stop("Dimensions of Z do not match W.")
        data_list$Z <- data_list$Z * data_list$W
    }

    #### Constant ILV ####
    ILV_datatypes <- c()
    ILV_names <- c()
    # identify ILVs from ILV_c
    if (!is.null(ILV_c)) {
        message("Constant ILV supplied.")
        ILV_c$id_numeric <- as.numeric(as.factor(ILV_c$id))
        # order in case user has not
        ILV_c <- ILV_c[order(ILV_c$id_numeric), ]
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
        ILV_tv$id_numeric <- as.numeric(as.factor(ILV_tv$id))
        ILV_tv$trial_numeric <- as.numeric(as.factor(ILV_tv$trial))
        ILV_tv$discrete_time <- with(
            ILV_tv, ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x)))
        )
        # order in case user has not
        ILV_tv <- ILV_tv[order(ILV_tv$trial_numeric, ILV_tv$id_numeric, ILV_tv$discrete_time), ]
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

    ## Network handling for high-resolution data
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
        networks$from_numeric <- temp_names$numeric[match(networks$from, temp_names$name)]
        networks$to_numeric <- temp_names$numeric[match(networks$to, temp_names$name)]
        networks$trial_numeric <- as.numeric(as.factor(networks$trial))
        if (!high_res) {
            networks$discrete_time <- with(networks, ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x))))
        } else {
            networks$discrete_time <- networks$time
        }
        is_symmetric <- nrow(networks[networks$trial_numeric == 1 & networks$discrete_time == min(networks$discrete_time), ]) == data_list$P * (data_list$P - 1)
        max_timesteps <- max(data_list$T)

        for (column in network_cols) {
            dims <- c(data_list$K, max_timesteps, data_list$P, data_list$P)
            A_matrix <- array(0, dim = dims)
            for (k in 1:data_list$K) {
                temp_df <- networks[networks$trial_numeric == k, ]
                for (i in 1:nrow(temp_df)) {
                    from <- as.integer(temp_df[i, "from_numeric"])
                    to <- as.integer(temp_df[i, "to_numeric"])
                    time <- as.integer(temp_df[i, "discrete_time"])
                    value <- as.numeric(temp_df[i, column])
                    A_matrix[k, time, from, to] <- value
                    if (!is_symmetric & network_type == "undirected") A_matrix[k, time, to, from] <- value
                }
                if (!is_dynamic) {
                    A_matrix[k, , , ] <- rep(A_matrix[k, 1, , ], each = dim(A_matrix)[2])
                }
                for (t in 1:max(data_list$T)) {
                    diag(A_matrix[k, t, , ]) <- 0
                }
            }
            if (min(A_matrix) < 0) stop("Edgeweights below zero detected. Rescale so that 0 = no connection.")
            data_list[[paste0("A_", column)]] <- A_matrix
        }
        data_list$network_names <- network_cols
    }
    #IF DISTRIBUTION
    else {
        # create and add network names
        network_names <- paste0("A_", seq_along(networks))
        data_list$network_names <- network_names

        for (i in seq_along(networks)) {
            net_obj <- networks[[i]]
            net_name <- network_names[i]

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
                    to   = as.integer(dyad_matrix[, 2]),
                    stringsAsFactors = FALSE
                )

                # Drop self-loops
                non_self_idx <- which(dyad_df$from != dyad_df$to)
                edges <- net_obj$edge_samples[, non_self_idx]

            } else if (inherits(net_obj, "STRAND Results Object")) {
                # STRAND handling
                ass_matrix <- net_obj$samples$predicted_network_sample  # [draws, from, to]
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

                edges <- qlogis(edge_mat)  # convert to logit space

            } else {
                stop("Unrecognized distributional network object. Expected bison_model or STRAND Results Object.")
            }

            # Get point estimate and covariance
            mu <- apply(edges, 2, median)
            cov_mat <- cov(edges)

            # Add to data list
            data_list[[paste0("logit_edge_mu_", net_name)]] <- mu
            data_list[[paste0("logit_edge_cov_", net_name)]] <- cov_mat
        }

        # Assume same number of dyads for all networks
        data_list$N_dyad <- length(mu)
    }

    ## Final sanity check and return ##
    dl_sanity_check(data_list = data_list)
    return(data_list)
}
