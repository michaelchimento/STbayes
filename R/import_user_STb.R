#' import_user_STb: Create STbayes_data object from user supplied data
#'
#' @param diffusion_data dataframe with columns id, trial, time, max_time
#' @param networks dataframe with columns trial, from, to, and one or more columns of edge weights named descriptively. Optionally an integer time column can be provided for dynamic network analysis, although networks must be provided for each time period between transmission events.
#' @param ILV_metadata optional dataframe with columns id, and any individual-level variables that might be of interest
#' @param ILVi Optional character vector of column names from ILV metadata to be considered when estimating asocial learning rate. If not specified, all ILV are applied to both.
#' @param ILVs Optional character vector of column names from ILV metadata to be considered when estimating social learning rate. If not specified, all ILV are applied to both.
#' @param ILVm Optional character vector of column names from ILV metadata to be considered in a multiplicative model.
#'
#' @return A list object containing properly formatted data to run social transmission models.
#' @export
#'
#' @examples
#' #very mock data
#' diffusion_data <- data.frame(
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
#'   sex = c(0, 1, 1, 0, 1, 0) # Factor ILVs must be input as numeric
#' )
#' imported_data <- import_user_STb(
#'   diffusion_data = diffusion_data,
#'   networks = networks,
#'   ILV_metadata = ILV_metadata,
#'   ILVi = c("age"), # Use only 'age' for asocial learning
#'   ILVs = c("sex") # Use only 'sex' for social learning
#' )
import_user_STb <- function(diffusion_data, networks, ILV_metadata=NULL, ILVi = NULL, ILVs = NULL, ILVm = NULL) {
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
  diffusion_data$index <- with(diffusion_data, ave(trial_numeric, FUN = seq_along))

  # create discrete time (this should be 0 if ID was a demo/seed)
  diffusion_data$discrete_time <- NA
  diffusion_data$discrete_time[diffusion_data$time != 0] <- with( # assign values only where time != 0
      diffusion_data[diffusion_data$time != 0, ],
      ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x)))
  )
  diffusion_data$discrete_time[diffusion_data$time == 0] = 0

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
  D_data <- D_data[D_data$time!=0,] # remove demos
  # Calculate the duration as time - lag(time) for each group
  D_data$duration <- with(D_data, ave(time, trial_numeric, FUN = function(x) c(x[1], diff(x))))
  D_data$duration <- ifelse(is.na(D_data$duration), D_data$time, D_data$duration)
  D_data <- D_data[, !(names(D_data) %in% "time")]
  # pivot wider
  D_data_real <- with(D_data, tapply(duration, list(trial_numeric, discrete_time), FUN = max, default = 0))
  D_data_int <- with(D_data, tapply(duration, list(trial_numeric, discrete_time), FUN = max, default = 0))

  # create a matrix where rows are trial_numeric, columns are id_numeric, and values are discrete_time
  t_data <- with(diffusion_data, tapply(discrete_time, list(trial_numeric, id_numeric), FUN = max, default = NA))
  # create a matrix where rows are trial_numeric, columns are id_numeric, and values are actual time measurements
  time_data <- with(diffusion_data, tapply(time, list(trial_numeric, id_numeric), FUN = max, default = NA))
  # gets the end of obs period of each trial
  time_max <- with(diffusion_data, tapply(max_time, list(trial_numeric), FUN = max, default = NA))
  # Replace NA values with -1 to denote individuals who did not participate in a trial
  #t_data[is.na(t_data)] <- -1

  # create a matrix where rows are trial_numeric, columns are id_numeric, and values are id_numeric
  id_data <- with(diffusion_data, tapply(id_numeric, list(trial_numeric, index), FUN = max, default = NA))

  data_list$K <- length(unique(diffusion_data$trial_numeric))
  data_list$Z <- length(unique(diffusion_data$id_numeric)) # number of distinct individuals
  data_list$N <- N_data$num_uncensored  # Uncensored counts per trial
  dim(data_list$N) = 1
  data_list$N_c <- N_data$num_censored  # Censored counts per trial
  dim(data_list$N_c) = 1
  data_list$T <- N_data$max_periods  # Max time periods (discrete)
  dim(data_list$T) = 1
  data_list$t <- t_data # vector of times (discrete time periods here) for uncensored
  data_list$T_max <- max(N_data$max_periods)
  data_list$time <- time_data
  data_list$time_max <- time_max
  data_list$Q <- max(diffusion_data$index) # number individuals per trial
  data_list$D <- D_data_real # duration data
  data_list$D_int = D_data_int #duration data in integer format (experimental)
  dim(data_list$D_int) <- dim(data_list$D) #why is r so annoying
  data_list$ind_id <- id_data # individual id data
  data_list$C = create_knowledge_matrix(diffusion_data) # knowledge state matrix
  dim(data_list$D)
  #### Metadata ####
  # identify ILVs from ILV_metadata
  if (!is.null(ILV_metadata)){
      message("ILV supplied.")
      ILV_metadata$id_numeric <- as.numeric(as.factor(ILV_metadata$id))
      # order in case user has not, assign indexes per trial
      diffusion_data <- diffusion_data[order(ILV_metadata$id_numeric), ]
      exclude_cols <- c("id", "id_numeric")
      ILV_metadata <- ILV_metadata[, !(names(ILV_metadata) %in% c("id", "id_numeric"))]
      ILV_cols <- names(ILV_metadata)

      for (column in names(ILV_metadata)) {
          data_list[[paste0("ILV_", column)]] <- ILV_metadata[[column]]
      }
      # write names
      if (is.null(ILVi)) {
          ILVi <- ILV_cols
      }
      if (is.null(ILVs)) {
          ILVs <- ILV_cols
      }
      if (is.null(ILVm)) {
          ILVm <- ILV_cols
      }
      data_list$ILVi_names <- ILVi
      data_list$ILVs_names <- ILVs
      data_list$ILVm_names <- ILVm
  } else {
      message("ILV not supplied.")
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

  # identify probable networks apart from cols below
  exclude_cols <- c("trial", "time", "from", "to")
  network_cols <- setdiff(names(networks), exclude_cols)

  networks$from_numeric <- as.numeric(networks$from)
  networks$to_numeric <- as.numeric(networks$to)
  networks$trial_numeric <- as.numeric(as.factor(networks$trial))
  networks$discrete_time <- with(networks, ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x))))

  # check if the same number of individuals are included in both datasets
  if (max(data_list$T) != max(networks$discrete_time) & is_dynamic) {
    message("Networks do not contain the same number of discrete timesteps as the diffusion data.")
  }

  is_symmetric <- nrow(networks[networks$trial_numeric == 1 & networks$discrete_time == 1, ]) == data_list$Z * (data_list$Z - 1)
  max_timesteps <- max(data_list$T)

  #for each network
  for (column in network_cols) {
    # initialize + populate matrix
    A_matrix <- array(0, dim = c(data_list$K, max(data_list$T), data_list$Z, data_list$Z))
    for (k in 1:data_list$K) {
      temp_df <- networks[networks$trial_numeric == k, ]
      for (i in 1:nrow(temp_df)) {
        from <- as.integer(temp_df[i, "from_numeric"])
        to <- as.integer(temp_df[i, "to_numeric"])
        time <- as.integer(temp_df[i, "discrete_time"])

        #fill in matrix
        A_matrix[k, time, from, to] <- as.numeric(temp_df[i, column])

        # if user has supplied non-symmetric data, mirror
        if (!is_symmetric) {
          A_matrix[k, time, to, from] <- as.numeric(temp_df[i, column])
        }

        if (!is_dynamic) {
            # Repeat the same network for each timestep
            for (t in 2:max_timesteps) {
                A_matrix[ k, t, , ] <- A_matrix[ k, 1 , , ]
            }
        }
      }

      # Zero the diagonal for each time step
      for (t in 1:max(data_list$T)) {
        diag(A_matrix[k, t, , ]) <- 0
      }
    }

    data_list[[paste0("A_", column)]] <- A_matrix
  }

  data_list$network_names <- network_cols

  data_list$N_veff = 2

  #sanity check
  dl_sanity_check(data_list=data_list)

  return(data_list)
}
