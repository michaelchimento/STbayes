#' Function to extract estimated acquisition times
#' @param fit CmdStanMCMC fit object
#' @param data_list STb_data object (must include ind_id and time)
#' @param prob interval for CI. defaults to .95
#' @param var_name variable name of acquisition time from GQ block (default = "acquisition_time")
#'
#' @return Dataframe of trial, individual, id, observed_time, mean/median/HPD of estimated time
#' @export
extract_acqTime <- function(fit, data_list, prob = .95, var_name = "acquisition_time") {
  if (!inherits(fit, c("CmdStanMCMC"))) {
    stop(sprintf("Model '%s' must be a CmdStanMCMC model fit.", name))
  }

  posterior_matrix <- fit$draws(variables = var_name, format = "draws_matrix")
  # Dimensions: iterations x chains x Q
  dim_post <- dim(posterior_matrix)
  if (length(dim_post) != 2) {
    stop("Expected array of shape [iterations x chains, individuals].")
  }

  # Get K and Q from data
  K <- dim(data_list$ind_id)[1]
  Q <- dim(data_list$ind_id)[2]

  # Preallocate matrices
  mean_acquisition <- matrix(0, K, Q)
  median_acquisition <- matrix(0, K, Q)
  lower_hpd <- matrix(0, K, Q)
  upper_hpd <- matrix(0, K, Q)

  # Loop through parameters by name
  for (trial in 1:K) {
    for (n in 1:Q) {
      param_name <- sprintf("%s[%d,%d]", var_name, trial, n)
      samples <- posterior_matrix[, param_name]

      mean_acquisition[trial, n] <- mean(samples, na.rm = TRUE)
      median_acquisition[trial, n] <- stats::median(samples, na.rm = TRUE)

      q <- posterior::quantile2(samples, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2))
      lower_hpd[trial, n] <- q[1]
      upper_hpd[trial, n] <- q[2]
    }
  }

  id_all <- as.vector(data_list$ind_id)
  observed_time_all <- as.vector(data_list$time)
  observed_time_ordered <- observed_time_all[match(id_all, seq_along(observed_time_all))]

  df <- data.frame(
    trial = rep(1:K, each = Q),
    individual = rep(1:Q, times = K),
    id = id_all,
    observed_time = observed_time_ordered,
    mean_time = as.vector(mean_acquisition),
    median_time = as.vector(median_acquisition),
    lower_hpd = as.vector(lower_hpd),
    upper_hpd = as.vector(upper_hpd)
  )

  return(df)
}
