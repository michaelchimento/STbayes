#' import_NBDA_STb: create STbayes_data object from nbda object
#'
#' @param nbda_object object of NBDAdata class
#' @param network_names character vector of descriptive names of networks you're importing.
#' @param multinetwork_s "separate" or "shared". If supplying more than one network, specify whether each network receives it's own s parameter, or shares a single parameter.
#' @param ILVi Optional character vector of ILVS to be considered when estimating intrinsic rate. If not specified, taken from NBDA object.
#' @param ILVs Optional character vector of ILVS to be considered when estimating social transmission rate. If not specified, taken from NBDA object.
#' @param ILVm Optional character vector of ILVS to be considered in a multiplicative model. If not specified, taken from NBDA object.
#'
#' @return A list containing properly formatted data to run social transmission models.
#' @export
#'
#' @examples
#' nbda_object <- nbdaData_cTADA <- STbayes::nbdaData_cTADA
#' data_list <- import_NBDA_STb(nbda_object, network_names = "assoc")
import_NBDA_STb <- function(nbda_object,
                            network_names = c("default"),
                            multinetwork_s = c("separate", "shared"),
                            ILVi = NULL,
                            ILVs = NULL,
                            ILVm = NULL,
                            high_res=F) {
  if (!.pkg_state$nbda_msg_shown) {
    message("👉 This function is convenient, but there's more flexibility and functionality available using the import_user_STb() function.")
    .pkg_state$nbda_msg_shown = TRUE
  }
  multinetwork_s <- match.arg(multinetwork_s)

  # Initialize list
  data_list <- list()

  #### Validate Input ####
  if (!inherits(nbda_object, "nbdaData")) {
    stop("😔 The provided object is not a valid nbdaData object.")
  }

  if (is.na(sum(nbda_object@timeAcq))) {
    message("🤔 This NBDA object is likely in OADA format.")
    event_data <- data.frame(
      id = nbda_object@orderAcq,
      trial = 1, # NBDA assumes a single trial (adjust if multi-trial support is added)
      time = seq_along(1:max(nbda_object@orderAcq)),
      max_time = max(nbda_object@orderAcq) + 1 # Maximum acquisition time
    )
  } else {
    message("🤔 This NBDA object is likely in TADA format.")
    event_data <- data.frame(
      id = nbda_object@orderAcq,
      trial = 1, # NBDA assumes a single trial (adjust if multi-trial support is added)
      time = nbda_object@timeAcq,
      max_time = nbda_object@endTime # Maximum acquisition time
    )
  }

  # Identify missing IDs
  all_ids <- c(1:dim(nbda_object@assMatrix)[1])
  learned_ids <- nbda_object@orderAcq
  censored_ids <- setdiff(all_ids, learned_ids)

  if (length(censored_ids > 0)) {
    # Create dataframe for censored individuals
    censored_data <- data.frame(
      id = censored_ids,
      trial = 1,
      time = nbda_object@endTime,
      max_time = nbda_object@endTime
    )
    # Combine both dataframes
    event_data <- rbind(event_data, censored_data)
  }



  #### event_data ####
  # Extract event times, IDs, and trial information


  # Convert IDs and trials to numeric values
  event_data$id_numeric <- as.numeric(as.factor(event_data$id))
  event_data$trial_numeric <- as.numeric(as.factor(event_data$trial))

  # Sort data and assign index within trials
  event_data <- event_data[order(event_data$trial_numeric, event_data$time), ]
  event_data$index <- with(event_data, ave(trial_numeric, trial_numeric, FUN = seq_along))

  # create discrete time (this should be 0 if ID was a demo/seed)
  event_data$discrete_time <- NA
  event_data$discrete_time[event_data$time != 0] <- with( # assign values only where time != 0
    event_data[event_data$time != 0, ],
    ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x)))
  )
  event_data$discrete_time[event_data$time == 0] <- 0

  # Identify seeds (demonstrators) and censored individuals
  event_data$seed <- with(event_data, ifelse(time == 0, 1, 0))
  event_data$censored <- with(event_data, ifelse(time < max_time, 0, 1))

  # Summarize trial data
  N_data <- do.call(rbind, lapply(split(event_data, event_data$trial_numeric), function(df) {
    data.frame(
      trial_numeric = unique(df$trial_numeric),
      num_uncensored = nrow(df) - sum(df$censored),
      num_censored = sum(df$censored),
      max_periods = max(df$discrete_time)
    )
  }))

  #### Duration Matrix ####
  D_data <- unique(event_data[order(event_data$trial_numeric, event_data$time), c("trial_numeric", "time", "discrete_time")])
  D_data <- D_data[D_data$time != 0, ] # remove demos
  D_data$duration <- with(D_data, ave(time, trial_numeric, FUN = function(x) c(x[1], diff(x))))
  D_data$duration <- ifelse(is.na(D_data$duration), D_data$time, D_data$duration)
  D_data <- D_data[, !(names(D_data) %in% "time")]
  # pivot wider
  D_data_real <- with(D_data, tapply(duration, list(trial_numeric, discrete_time), FUN = max, default = 0))
  D_data_int <- with(D_data, tapply(as.integer(duration), list(trial_numeric, discrete_time), FUN = max, default = 0))
  dimnames(D_data_real) <- NULL
  dimnames(D_data_int) <- NULL

  #### Discrete Time Matrix ####
  # create a matrix where rows are trial_numeric, columns are id_numeric, and values are discrete_time
  t_data <- with(event_data, tapply(discrete_time, list(trial_numeric, id_numeric), FUN = max, default = -1))
  # create a matrix where rows are trial_numeric, columns are id_numeric, and values are actual time measurements
  time_data <- with(event_data, tapply(time, list(trial_numeric, id_numeric), FUN = max, default = -1))
  # gets the end of obs period of each trial
  time_max <- with(event_data, tapply(max_time, list(trial_numeric), FUN = max, default = -1))

  #### Individual IDs Matrix ####
  # create a matrix where rows are trial_numeric, columns are N + N_c, and values are id_numeric
  # values are then left aligned such that individuals who were in each trial will appear in the first indices
  # i know this is a really weird way of doing a ragged list, in case there are differing numbers of individuals in each trial... but should work fine
  # actually not sure whether this will ever be an issue with nbda data objs..
  id_data <- with(event_data, tapply(id_numeric, list(trial_numeric, index), FUN = max, default = NA))
  id_data <- t(apply(id_data, 1, function(row) {
    c(stats::na.omit(row), rep(-1, sum(is.na(row))))
  }))

  #### Populate Data List ####
  data_list <- list(
    K = length(unique(event_data$trial_numeric)), # Number of trials
    P = length(unique(event_data$id_numeric)), # Number of individuals
    N = N_data$num_uncensored, # Uncensored counts per trial
    N_c = N_data$num_censored, # Censored counts per trial
    T = N_data$max_periods, # Max time periods (discrete)
    t = t_data, # Discrete time matrix
    T_max = max(N_data$max_periods),
    time = time_data,
    time_max = time_max,
    Q = max(event_data$index), # Max individuals per trial
    D = D_data_real, # duration data
    D_int = D_data_int, # duration data in integer format (experimental)
    ind_id = id_data, # Individual IDs matrix
    Zn = create_Z_matrix(event_data, high_res = F), # knowledge state matrix
    multinetwork_s = multinetwork_s,
    prop_k = NULL, #not used by nbda imports
    high_res = F,  #not used by nbda imports
    directed = F  #not used by nbda imports
  )
  dim(data_list$N) <- length(data_list$N)
  dim(data_list$N_c) <- length(data_list$N_c)
  dim(data_list$T) <- length(data_list$T)
  dim(data_list$D_int) <- dim(data_list$D) # why is r so annoying
  data_list$Z <- sweep(data_list$Zn, MARGIN = 3, STATS = nbda_object@weights, FUN = "*") # mult with weights

  #### ILV Metadata ####

  if (nbda_object@asocialTreatment == "constant") {
    ILV_metadata <- data.frame(id = nbda_object@idname)
    if (nbda_object@asoc_ilv != "ILVabsent") {
      ILV_metadata <- cbind(ILV_metadata, extract_ILV(nbda_object, "asocILVdata"))
    }
    if (nbda_object@int_ilv != "ILVabsent") {
      ILV_metadata <- cbind(ILV_metadata, extract_ILV(nbda_object, "intILVdata"))
    }
    if (nbda_object@multi_ilv != "ILVabsent") {
      ILV_metadata <- cbind(ILV_metadata, extract_ILV(nbda_object, "multiILVdata"))
    }

    ILV_metadata <- ILV_metadata[, !(names(ILV_metadata) %in% "id")]
    ILV_metadata <- remove_duplicate_columns(ILV_metadata)
    ILV_cols <- names(ILV_metadata)
    for (col in ILV_cols) {
      data_list[[paste0("ILV_", col)]] <- ILV_metadata[[col]]
    }
  } else {
    #### Time-varying ILV ####
    message("Time-varying ILV supplied.")
    ILV_tv <- data.frame(id = rep(nbda_object@idname, times = max(nbda_object@orderAcq)), time = rep(1:max(nbda_object@orderAcq), each = data_list$P))
    ILV_tv$trial <- 1
    if (!"ILVabsent" %in% nbda_object@asoc_ilv) {
      ILV_tv <- cbind(ILV_tv, extract_tv_ILV(nbda_object, "asocILVdata"))
    }
    if (!"ILVabsent" %in% nbda_object@int_ilv) {
      ILV_tv <- cbind(ILV_tv, extract_tv_ILV(nbda_object, "intILVdata"))
    }
    if (!"ILVabsent" %in% nbda_object@multi_ilv) {
      ILV_tv <- cbind(ILV_tv, extract_tv_ILV(nbda_object, "multiILVdata"))
    }

    ILV_tv <- remove_duplicate_columns(ILV_tv)
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
    # loop through each ILV column and add to datalist
    for (col in ILV_cols) {
      # reshape data into matrix
      mat <- with(ILV_tv, tapply(ILV_tv[[col]], list(trial, discrete_time, id), FUN = mean, simplify = TRUE))
      mat[is.na(mat)] <- 0
      data_list[[paste0("ILV_", col)]] <- mat
    }
  }

  if (is.null(ILVi)) {
    ILVi <- nbda_object@asoc_ilv
  }
  if (is.null(ILVs)) {
    ILVs <- nbda_object@int_ilv
  }
  if (is.null(ILVm)) {
    ILVm <- nbda_object@multi_ilv
  }

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
    network <- networks[, , i, ] # [n, n, 1]

    # if only one time period, replicate the network across all timesteps
    if (num_time_periods == 1) {
      # repeat the same network for each timestep
      network <- array(network, dim = c(dim(network)[1], dim(network)[2], max_timesteps))
      for (t in 2:max_timesteps) { # Explicitly copy for safety
        network[, , t] <- network[, , 1]
      }
    }

    # reorder dimensions to [time, row, col]
    network <- aperm(network, perm = c(3, 1, 2))

    # zero out self loops
    for (t in 1:max(data_list$T)) {
      diag(network[t, , ]) <- 0
    }

    # add the first dimension for trials (k = 1)
    network <- array(network, dim = c(1, dim(network))) # [k = 1, t, n, n]

    data_list[[paste0("A_", network_name)]] <- network
  }

  data_list$network_names <- network_names

  # convert to new user_STb format, don't really want to touch the code above
  # After your current network processing loop
  network_arrays <- lapply(network_names, function(nm) data_list[[paste0("A_", nm)]])
  num_networks <- length(network_arrays)
  num_trials <- dim(network_arrays[[1]])[1]
  num_timesteps <- dim(network_arrays[[1]])[2]
  num_individuals <- dim(network_arrays[[1]])[3]

  # Create unified A array: [network, trial, time, n, n]
  data_list$A <- array(0, dim = c(num_networks, num_trials, num_timesteps, num_individuals, num_individuals))

  for (i in seq_along(network_arrays)) {
    data_list$A[i, , , , ] <- network_arrays[[i]]
  }

  for (nm in network_names) {
    data_list[[paste0("A_", nm)]] <- NULL
  }

  data_list$high_res <- F
  data_list$N_networks <- num_networks

  if (data_list$multinetwork_s == "separate" & data_list$N_networks == 1) {
    data_list$multinetwork_s <- "shared"
  }

  #### Output Messages ####
  dl_sanity_check(data_list = data_list)

  return(data_list)
}
