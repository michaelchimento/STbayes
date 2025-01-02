#' import_NBDA_STb: create STbayes_data object from nbda object
#'
#' @param diffusion_data dataframe with columns id, trial, time, max_time
#' @param networks dataframe with columns trial, from, to, and one or more columns of edge weights named descriptively. Optionally an integer time column can be provided for dynamic network analysis, although networks must be provided for each time period between transmission events.
#' @param ILV_metadata dataframe with columns id, and any individual-level variables that might be of interest
#' @param network_names character vector of descriptive names of networks you're importing.
#' @param ILVi Optional character vector of ILVS to be considered when estimating asocial learning rate. If not specified, taken from NBDA object.
#' @param ILVs Optional character vector of ILVS to be considered when estimating social learning rate. If not specified, taken from NBDA object.
#' @param ILVm Optional character vector of ILVS to be considered in a multiplicative model. If not specified, taken from NBDA object.
#'
#' @return A list containing properly formatted data to run social transmission models.
#' @export
#'
#' @examples
#' nbda_object<-nbdaData_cTADA <- STbayes::nbdaData_cTADA
#' data_list <- import_NBDA_STb(nbda_object, network_names="assoc")

import_NBDA_STb <- function(nbda_object, network_names= c("default"), ILVi = NULL, ILVs = NULL, ILVm = NULL) {
    # Initialize list
    data_list <- list()

    #### Validate Input ####
    if (!inherits(nbda_object, "nbdaData")) {
        stop("The provided object is not a valid nbdaData object.")
    }

    if (is.na(sum(nbda_object@timeAcq))){
        message("This NBDA object is likely in OADA format.")
        diffusion_data <- data.frame(
            id = nbda_object@orderAcq,
            trial = 1,  # NBDA assumes a single trial (adjust if multi-trial support is added)
            time = seq_along(1:max(nbda_object@orderAcq)),
            max_time = max(nbda_object@orderAcq)  # Maximum acquisition time
        )
    } else {
        message("This NBDA object is likely in TADA format.")
        diffusion_data <- data.frame(
            id = nbda_object@orderAcq,
            trial = 1,  # NBDA assumes a single trial (adjust if multi-trial support is added)
            time = nbda_object@timeAcq,
            max_time = nbda_object@endTime  # Maximum acquisition time
        )
    }

    #### Diffusion_data ####
    # Extract event times, IDs, and trial information


    # Convert IDs and trials to numeric values
    diffusion_data$id_numeric <- as.numeric(as.factor(diffusion_data$id))
    diffusion_data$trial_numeric <- as.numeric(as.factor(diffusion_data$trial))

    # Sort data and assign index within trials
    diffusion_data <- diffusion_data[order(diffusion_data$trial_numeric, diffusion_data$time), ]
    diffusion_data$index <- with(diffusion_data, ave(trial_numeric, FUN = seq_along))

    # create discrete time (this should be 0 if ID was a demo/seed)
    diffusion_data$discrete_time <- NA
    diffusion_data$discrete_time[diffusion_data$time != 0] <- with( # assign values only where time != 0
        diffusion_data[diffusion_data$time != 0, ],
        ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x)))
    )
    diffusion_data$discrete_time[diffusion_data$time == 0] = 0

    # Identify seeds (demonstrators) and censored individuals
    diffusion_data$seed <- with(diffusion_data, ifelse(time == 0, 1, 0))
    diffusion_data$censored <- with(diffusion_data, ifelse(time < max_time, 0, 1))

    # Summarize trial data
    N_data <- do.call(rbind, lapply(split(diffusion_data, diffusion_data$trial_numeric), function(df) {
        data.frame(
            trial_numeric = unique(df$trial_numeric),
            num_uncensored = nrow(df) - sum(df$censored),
            num_censored = sum(df$censored),
            max_periods = max(df$discrete_time)
        )
    }))

    #### Duration Matrix ####
    D_data <- unique(diffusion_data[order(diffusion_data$trial_numeric, diffusion_data$time), c("trial_numeric", "time", "discrete_time")])
    D_data <- D_data[D_data$time!=0,] # remove demos
    D_data$duration <- with(D_data, ave(time, trial_numeric, FUN = function(x) c(x[1], diff(x))))
    D_data$duration <- ifelse(is.na(D_data$duration), D_data$time, D_data$duration)
    D_data <- D_data[, !(names(D_data) %in% "time")]
    # pivot wider
    D_data_real <- with(D_data, tapply(duration, list(trial_numeric, discrete_time), FUN = max, default = 0))
    D_data_int <- as.integer(with(D_data, tapply(duration, list(trial_numeric, discrete_time), FUN = max, default = 0)))

    #### Discrete Time Matrix ####
    # create a matrix where rows are trial_numeric, columns are id_numeric, and values are discrete_time
    t_data <- with(diffusion_data, tapply(discrete_time, list(trial_numeric, id_numeric), FUN = max, default = -1))
    # create a matrix where rows are trial_numeric, columns are id_numeric, and values are actual time measurements
    time_data <- with(diffusion_data, tapply(time, list(trial_numeric, id_numeric), FUN = max, default = -1))
    # gets the end of obs period of each trial
    time_max <- with(diffusion_data, tapply(max_time, list(trial_numeric), FUN = max, default = -1))
    #t_data[is.na(t_data)] <- -1

    #### Individual IDs Matrix ####
    # create a matrix where rows are trial_numeric, columns are N + N_c, and values are id_numeric
    # values are then left aligned such that individuals who were in each trial will appear in the first indices
    # i know this is a really weird way of doing a ragged list, in case there are differing numbers of individuals in each trial... but should work fine
    # actually not sure whether this will ever be an issue with nbda data objs..
    id_data <- with(diffusion_data, tapply(id_numeric, list(trial_numeric, index), FUN = max, default = NA))
    id_data <- t(apply(id_data, 1, function(row) {
        c(na.omit(row), rep(-1, sum(is.na(row))))
    }))

    #### Populate Data List ####
    data_list$K <- length(unique(diffusion_data$trial_numeric))  # Number of trials
    data_list$Z <- length(unique(diffusion_data$id_numeric))  # Number of individuals
    data_list$N <- N_data$num_uncensored  # Uncensored counts per trial
    dim(data_list$N) = length(data_list$N)
    data_list$N_c <- N_data$num_censored  # Censored counts per trial
    dim(data_list$N_c) = length(data_list$N_c)
    data_list$T <- N_data$max_periods  # Max time periods (discrete)
    dim(data_list$T) = length(data_list$T)
    data_list$t <- t_data  # Discrete time matrix
    data_list$T_max <- max(N_data$max_periods)
    data_list$time <- time_data
    data_list$time_max <- time_max
    data_list$Q <- max(diffusion_data$index)  # Max individuals per trial
    data_list$D <- D_data_real # duration data
    data_list$D_int = D_data_int #duration data in integer format (experimental)
    dim(data_list$D_int) <- dim(data_list$D) #why is r so annoying
    data_list$ind_id <- id_data  # Individual IDs matrix
    data_list$C = create_knowledge_matrix(diffusion_data) # knowledge state matrix

    #### ILV Metadata ####
    ILV_metadata <- data.frame(id = nbda_object@idname)
    if (nbda_object@asoc_ilv != "ILVabsent"){ILV_metadata <- cbind(ILV_metadata, extract_ILV(nbda_object, "asocILVdata"))}
    if (nbda_object@int_ilv != "ILVabsent"){ILV_metadata <- cbind(ILV_metadata, extract_ILV(nbda_object, "intILVdata"))}
    if (nbda_object@multi_ilv != "ILVabsent"){ILV_metadata <- cbind(ILV_metadata, extract_ILV(nbda_object, "multiILVdata"))}
    ILV_metadata <- ILV_metadata[, !(names(ILV_metadata) %in% "id")]
    ILV_metadata = remove_duplicate_columns(ILV_metadata)
    ILV_cols <- names(ILV_metadata)

    for (column in ILV_cols) {
        data_list[[column]] <- ILV_metadata[[column]]
    }

    if (is.null(ILVi)) {ILVi <- nbda_object@asoc_ilv}
    if (is.null(ILVs)) {ILVs <- nbda_object@int_ilv}
    if (is.null(ILVm)) {ILVm <- nbda_object@multi_ilv}

    data_list$ILVi_names <- ILVi
    data_list$ILVs_names <- ILVs
    data_list$ILVm_names <- ILVm

    #### Networks ####
    # extract network from nbdaData object
    networks <- nbda_object@assMatrix

    # Extract the number of networks (k) and time periods (t)
    dim_networks <- dim(networks)
    num_networks <- dim(networks)[3]
    num_time_periods <- dim(networks)[4]

    # Max timesteps
    max_timesteps <- max(data_list$T)

    for (i in seq_len(num_networks)) {
        network_name <- network_names[i]

        # Extract the network for the current edge type
        network <- networks[, , i, ]  # [n, n, 1]

        # if only one time period, replicate the network across all timesteps
        if (num_time_periods == 1) {
            # repeat the same network for each timestep
            network <- array(network, dim = c(dim(network)[1], dim(network)[2], max_timesteps))
            for (t in 2:max_timesteps) {  # Explicitly copy for safety
                network[ , , t] <- network[ , , 1]
            }
            network <- aperm(network, perm = c(3, 1, 2))
        } else {
            # reorder dimensions to [time, row, col]
            network <- aperm(network, perm = c(3, 1, 2))  # [t, n, n]
        }

        # add the first dimension for trials (k = 1)
        network <- array(network, dim = c(1, dim(network)))  # [k = 1, t, n, n]

        data_list[[paste0("A_", network_name)]] <- network
    }

    data_list$network_names <- network_names

    #### Output Messages ####
    dl_sanity_check(data_list=data_list)

    return(data_list)
}
