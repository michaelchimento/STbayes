#' import_user_STb: Create STbayes_data object from user supplied data
#'
#' @param diffusion_data dataframe with columns id, trial, time, max_time
#' @param networks dataframe with columns trial, from, to, and one or more columns of edge weights named descriptively. Optionally an integer time column can be provided for dynamic network analysis, although networks must be provided for each time period between transmission events.
#' @param ILV_c optional dataframe with columns id, and any constant individual-level variables that might be of interest
#' @param ILV_tv optional dataframe with columns trial, id, time and any time-varying variables. Variable values should summarize the variable for each inter-acquisition period.
#' @param ILVi Optional character vector of column names from ILV metadata to be considered when estimating asocial learning rate. If not specified, all ILV are applied to both.
#' @param ILVs Optional character vector of column names from ILV metadata to be considered when estimating social learning rate. If not specified, all ILV are applied to both.
#' @param ILVm Optional character vector of column names from ILV metadata to be considered in a multiplicative model.
#'
#' @return A list object containing properly formatted data to run social transmission models.
#' @export
#'
#' @examples
#' # very mock data
#' diffusion_data <- data.frame(
#'   trial = rep(1:2, each = 3),
#'   id = LETTERS[1:6],
#'   time = c(0, 1, 2, 0, 1, 4),
#'   max_time = c(3, 3, 3, 4, 4, 4)
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
#' imported_data <- import_user_STb(
#'   diffusion_data = diffusion_data,
#'   networks = networks,
#'   ILV_c = ILV_c,
#'   ILV_tv = ILV_tv,
#'   ILVi = c("age", "dist_from_task"), # Use 'age' and time-varying ILV 'dist_from_task' for asocial learning
#'   ILVs = c("sex"), # Use only 'sex' for social learning
#'   ILVm = c("weight") # Use weight for multiplicative effect on asocial and social learning
#' )
#'
import_user_STb <- function(diffusion_data, networks, ILV_c = NULL, ILV_tv = NULL, ILVi = NULL, ILVs = NULL, ILVm = NULL) {
  # warnings
  if (is.null(ILVi) & is.null(ILVs) & is.null(ILVm) & (!is.null(ILV_c) | !is.null(ILV_tv))) {
    message("WARNING: You have provided ILV, yet did not specify whether they should be additive or multiplicative (missing arguments ILVi, ILVs, ILVm). STbayes defaults to unconstrained additive (i.e. each variable's effect on both asocial and social learning will be separately estimated). It is recommended to explicitly define this.")
  }

  if (is.null(ILVm)) {
    ILVm <- "ILVabsent"
  }

  # Initialize list
  data_list <- list()

  #### Diffusion_data ####
  # diffusion_data should be in format id, trial, time, max_time
  # if time==0, assume to be trained demonstrator, if time==max_time, assume to be censored

  # create numeric variables in case user has supplied strings
  diffusion_data$id_numeric <- as.numeric(as.factor(diffusion_data$id))
  diffusion_data$trial_numeric <- as.numeric(as.factor(diffusion_data$trial))

  # order in case user has not, assign indexes per trial
  diffusion_data <- diffusion_data[order(diffusion_data$trial_numeric, diffusion_data$time), ]
  diffusion_data$index <- with(diffusion_data, ave(trial_numeric, trial_numeric, FUN = seq_along))

  # create discrete time (this should be 0 if ID was a demo/seed)
  diffusion_data$discrete_time <- NA
  diffusion_data$discrete_time[diffusion_data$time != 0] <- with( # assign values only where time != 0
    diffusion_data[diffusion_data$time != 0, ],
    ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x)))
  )
  diffusion_data$discrete_time[diffusion_data$time == 0] <- 0

  # group by time/trial and see if there's more than 1
  diffusion_data$tie <- with(diffusion_data, ave(time, interaction(trial_numeric, time), FUN = function(x) length(x) > 1))

  # identify seeds/demonstrators (where time == 0)
  diffusion_data$seed <- with(diffusion_data, ifelse(time == 0, 1, 0))

  # identify censored rows (where time == obs_period)
  diffusion_data$censored <- with(diffusion_data, ifelse(time < max_time, 0, 1))

  # summarize censored/uncensored for each trial
  N_data <- do.call(rbind, lapply(split(diffusion_data, diffusion_data$trial_numeric), function(df) {
    data.frame(
      trial_numeric = unique(df$trial_numeric),
      num_uncensored = nrow(df) - sum(df$censored),
      num_censored = sum(df$censored),
      max_periods = max(df$discrete_time)
    )
  }))

  # create matrix of where rows = trial and columns = discrete_time, and values = duration
  # summarize by trial and time
  D_data <- unique(diffusion_data[order(diffusion_data$trial_numeric, diffusion_data$time), c("trial_numeric", "time", "discrete_time")])
  D_data <- D_data[D_data$time != 0, ] # remove demos
  # Calculate the duration as time - lag(time) for each group
  D_data$duration <- with(D_data, ave(time, trial_numeric, FUN = function(x) c(x[1], diff(x))))
  D_data$duration <- ifelse(is.na(D_data$duration), D_data$time, D_data$duration)
  D_data <- D_data[, !(names(D_data) %in% "time")]
  # pivot wider
  D_data_real <- with(D_data, tapply(duration, list(trial_numeric, discrete_time), FUN = max, default = 0))
  D_data_int <- with(D_data, tapply(duration, list(trial_numeric, as.integer(discrete_time)), FUN = max, default = 0))

  # create a matrix where rows are trial_numeric, columns are id_numeric, and values are discrete_time
  t_data <- with(diffusion_data, tapply(discrete_time, list(trial_numeric, id_numeric), FUN = max, default = -1))
  # create a matrix where rows are trial_numeric, columns are id_numeric, and values are actual time measurements
  time_data <- with(diffusion_data, tapply(time, list(trial_numeric, id_numeric), FUN = max, default = -1))
  # gets the end of obs period of each trial
  time_max <- with(diffusion_data, tapply(max_time, list(trial_numeric), FUN = max, default = -1))
  # Replace NA values with -1 to denote individuals who did not participate in a trial
  # t_data[is.na(t_data)] <- -1

  # create a matrix where rows are trial_numeric, columns are N + N_c, and values are id_numeric
  # values are then left aligned such that individuals who were in each trial will appear in the first indices
  # i know this is a really weird way of doing a ragged list, in case there are differing numbers of individuals in each trial... but should work fine
  id_data <- with(diffusion_data, tapply(id_numeric, list(trial_numeric, index), FUN = max, default = NA))
  id_data <- t(apply(id_data, 1, function(row) {
    c(na.omit(row), rep(-1, sum(is.na(row))))
  }))

  data_list$K <- length(unique(diffusion_data$trial_numeric))
  data_list$Z <- length(unique(diffusion_data$id_numeric)) # number of distinct individuals
  data_list$N <- N_data$num_uncensored # Uncensored counts per trial
  dim(data_list$N) <- length(data_list$N)
  data_list$N_c <- N_data$num_censored # Censored counts per trial
  dim(data_list$N_c) <- length(data_list$N_c)
  data_list$T <- N_data$max_periods # Max time periods (discrete)
  dim(data_list$T) <- length(data_list$T)
  data_list$t <- t_data # vector of times (discrete time periods here) for uncensored
  data_list$T_max <- max(N_data$max_periods)
  data_list$time <- time_data
  data_list$time_max <- time_max
  data_list$Q <- max(diffusion_data$index) # number individuals per trial
  data_list$D <- D_data_real # duration data
  data_list$D_int <- D_data_int # duration data in integer format (experimental)
  dim(data_list$D_int) <- dim(data_list$D) # why is r so annoying
  data_list$ind_id <- id_data # individual id data
  data_list$C <- create_knowledge_matrix(diffusion_data) # knowledge state matrix

  #### Constant ILV ####
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
      # reshape data into matrix
      mat <- with(ILV_tv, tapply(ILV_tv[[col]], list(trial, discrete_time, id), FUN = mean, simplify = TRUE))
      mat[is.na(mat)] <- 0
      data_list[[paste0("ILV_", col)]] <- mat
    }
  }

  if (!is.null(ILV_c) | !is.null(ILV_tv)) {
    # write names
    if (is.null(ILVi)) {
      data_list$ILVi_names <- ILV_names
    } else {
      data_list$ILVi_names <- ILVi
    }
    if (is.null(ILVs)) {
      data_list$ILVs_names <- ILV_names
    } else {
      data_list$ILVs_names <- ILVs
    }
    data_list$ILVm_names <- ILVm
  } else {
    message("No ILV supplied.")
    data_list$ILVi_names <- "ILVabsent"
    data_list$ILVs_names <- "ILVabsent"
    data_list$ILVm_names <- "ILVabsent"
  }

  #### deal with networks ####
  # network_data should be in format trial, time, from, to, network_value1, network_value2, ... etc
  # as many networks as wanted can be supplied. names of networks will be taken from column name.
  # if no time is provided, it will be assumed that this is a static network analysis

  # check if the same number of individuals are included in both datasets
  if (data_list$Z != length(unique(c(networks$from, networks$to)))) {
    message("Networks do not contain the same number of unique individuals as the diffusion data.")
  }

  # check if the same number of individuals are included in both datasets
  if (data_list$K != length(unique(networks$trial))) {
    message("Networks do not contain the same number of trials as the diffusion data.")
  }

  # check if dynamic networks are supplied
  is_dynamic <- "time" %in% names(networks)
  if (!is_dynamic) {
    message("User input indicates a static network. If dynamic network, please include 'time' column.")
    networks$time <- 1
  }

  # check if dynamic networks are supplied
  is_distribution <- "draw" %in% names(networks)
  if (is_distribution) {
    message("User supplied edge weights as distributions. This currently only works when K=1 trials.")
    if (!is_dynamic) networks$time <- 1
  } else {
    message("User supplied edge weights as point estimates.")
  }

  # identify probable networks apart from cols below
  exclude_cols <- c("trial", "time", "from", "to", "draw")
  network_cols <- setdiff(names(networks), exclude_cols)

  temp_names <- data.frame(
    name = sort(unique(c(networks$from, networks$to))),
    numeric = seq_along(sort(unique(c(networks$from, networks$to))))
  )

  networks$from_numeric <- temp_names$numeric[match(networks$from, temp_names$name)]
  networks$to_numeric <- temp_names$numeric[match(networks$to, temp_names$name)]
  networks$trial_numeric <- as.numeric(as.factor(networks$trial))
  networks$discrete_time <- with(networks, ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x))))

  # check if the same number of individuals are included in both datasets
  if (max(data_list$T) != max(networks$discrete_time) & is_dynamic) {
    message("Networks do not contain the same number of discrete timesteps as the diffusion data.")
  }

  if (is_distribution){
    is_symmetric <- nrow(networks[networks$trial_numeric == 1 & networks$discrete_time == 1 & networks$draw == 1, ]) == data_list$Z * (data_list$Z - 1)
  } else {
    is_symmetric <- nrow(networks[networks$trial_numeric == 1 & networks$discrete_time == 1, ]) == data_list$Z * (data_list$Z - 1)
  }
  # precompute constants
  max_timesteps <- max(data_list$T)
  max_draw <- if (is_distribution) max(networks$draw) else NULL

  # for each network
  for (column in network_cols) {
    # initialize + populate matrix
    dims <- if (is_distribution) {
      c(data_list$K, max_timesteps, max_draw, data_list$Z, data_list$Z)
    } else {
      c(data_list$K, max_timesteps, data_list$Z, data_list$Z)
    }
    A_matrix <- array(0, dim = dims)

    for (k in 1:data_list$K) {
      temp_df <- networks[networks$trial_numeric == k, ]

      for (i in 1:nrow(temp_df)) {
        from <- as.integer(temp_df[i, "from_numeric"])
        to <- as.integer(temp_df[i, "to_numeric"])
        time <- as.integer(temp_df[i, "discrete_time"])
        draw <- if (is_distribution) as.integer(temp_df[i, "draw"]) else NULL
        value <- as.numeric(temp_df[i, column])

        # fill in matrix
        if (is_distribution){
          A_matrix[k, time, draw, from, to] <- value
          if (!is_symmetric) A_matrix[k, time, draw, to, from] <- value
        } else{
          A_matrix[k, time, from, to] <- value
          if (!is_symmetric) A_matrix[k, time, to, from] <- value
        }
      }

      if (!is_dynamic) {
          if (is_distribution){
            A_matrix[k, , , , ] <- rep(A_matrix[k, 1, , , ], each = dim(A_matrix)[2])
          } else {
            A_matrix[k,,,] <- rep(A_matrix[k, 1, , ], each = dim(A_matrix)[2])
          }
      }

      # Zero the diagonal for each time step and draw
      if (is_distribution){
        for (t in 1:max(data_list$T)) {
          for (d in 1:max(networks$draw)) {
            diag(A_matrix[k, t, d, , ]) <- 0
          }
        }
      } else {
        for (t in 1:max(data_list$T)) {
          diag(A_matrix[k, t, , ]) <- 0
        }
      }
    }
    if (min(A_matrix) < 0) stop("Edgeweights below zero detected, not allowed. Rescale such that 0 is no connection.")
    data_list[[paste0("A_", column)]] <- A_matrix
  }

  data_list$network_names <- network_cols
  if (is_distribution) data_list$S = max(networks$draw)

  # sanity check
  dl_sanity_check(data_list = data_list)

  return(data_list)
}
