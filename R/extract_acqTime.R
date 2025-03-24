#' Function to extract estimated acquisition times
#'
#' @param fit rstan fit object
#' @param data_list STb_data object
#' @param var_name name of where acquisition time is stored in rstan fit
#'
#' @return dataframe of ids, estimated timestep + HPD levels and observed timestep
#' @export
#'
#' @examples
extract_acqTime <- function(fit, data_list, var_name = "acquisition_time") {
    # Extract posterior samples for 'acquisition_time'
    posterior_samples <- rstan::extract(fit, pars = var_name, permuted = TRUE)[[var_name]]

    # Check dimensions
    dim_samples <- dim(posterior_samples) # Should be iterations x K x Q
    if (length(dim_samples) != 3) stop("Variable does not have the expected dimensions (iterations, K, Q).")

    iterations <- dim_samples[1]
    K <- dim_samples[2]
    Q <- dim_samples[3]

    # Flatten and summarize acquisition times
    mean_acquisition <- matrix(0, K, Q)
    median_acquisition <- matrix(0, K, Q)
    lower_hpd <- matrix(0, K, Q)
    upper_hpd <- matrix(0, K, Q)

    # Compute mean and HPD intervals for each trial and individual
    for (trial in 1:K) {
        for (n in 1:Q) {
            id_intrial = rep(n, K)
            samples <- posterior_samples[, trial, n]
            mean_acquisition[trial, n] <- mean(samples, na.rm = TRUE)
            median_acquisition[trial, n] <- median(samples, na.rm = TRUE)
            hpd <- coda::HPDinterval(coda::as.mcmc(samples), prob = 0.95)
            lower_hpd[trial, n] <- hpd[1]
            upper_hpd[trial, n] <- hpd[2]
        }
    }

    id_all <- as.vector(data_list$ind_id)
    observed_time_all <- as.vector(data_list$time)

    # need to reorder observed times to match the order of IDs in ind_id
    observed_time_ordered <- observed_time_all[match(id_all, seq_along(observed_time_all))]

    # Reshape the data for ggplot
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

