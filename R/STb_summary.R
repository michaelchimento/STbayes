#' Return summary table
#'
#' @param fit model fit from fit_STb_model
#' @param depth integer depth of multidimensional parameters to extract
#' @param prob double limits for HPD of estimates
#' @param ignore_params character vector of parameters to ignore
#' @param digits integer of digits to round to
#'
#' @return Summary table
#' @importFrom rstan extract Rhat
#' @importFrom coda as.mcmc HPDinterval effectiveSize
#' @export
#'
#' @examples
#' summary_table <- STb_summary(fit)
STb_summary <- function(fit, depth = 1, prob = 0.95, ignore_params = c("lp__", "idx", "log_lik", "log_lik_matrix", "acquisition_time", "z_ID", "Rho_ID","v_ID"), digits = 3) {
    # extract posterior samples
    samples <- rstan::extract(fit, permuted = TRUE)

    # filter by depth and ignore specific parameters
    param_names <- names(samples)
    param_depths <- sapply(samples, function(x) length(dim(x)))
    selected_params <- param_names[param_depths <= depth & !param_names %in% ignore_params]

    transformed_params = list()

    summary_table <- data.frame(
        Parameter = character(),
        Median = numeric(),
        HPDI_Lower = numeric(),
        HPDI_Upper = numeric(),
        n_eff = numeric(),
        Rhat = numeric(),
        stringsAsFactors = FALSE
    )

    # process each parameter
    for (param in selected_params) {
        if (is.matrix(samples[[param]])) {
            for (d in 1:dim(samples[[param]])[2]) {
                param_samples <- samples[[param]][, d]
                param_samples <- coda::as.mcmc(as.matrix(param_samples))

                # compute statistics
                mean_value <- mean(param_samples)
                median_value <- stats::median(param_samples)
                hpdi <- coda::HPDinterval(param_samples, prob = prob)
                n_eff <- coda::effectiveSize(param_samples)
                rhat <- rstan::Rhat(as.matrix(param_samples))

                # append
                summary_table <- rbind(summary_table, data.frame(
                    Parameter = paste0(param,"[",d,"]"),
                    Mean = mean_value,
                    Median = median_value,
                    HPDI_Lower = hpdi[1, "lower"],
                    HPDI_Upper = hpdi[1, "upper"],
                    n_eff = n_eff,
                    Rhat = rhat,
                    stringsAsFactors = FALSE
                ))
            }
        } else {
            param_samples <- samples[[param]]
            param_samples <- coda::as.mcmc(as.matrix(param_samples))

            # compute statistics
            mean_value <- mean(param_samples)
            median_value <- stats::median(param_samples)
            hpdi <- coda::HPDinterval(param_samples, prob = prob)
            n_eff <- coda::effectiveSize(param_samples)
            rhat <- rstan::Rhat(as.matrix(param_samples))

            # append
            summary_table <- rbind(summary_table, data.frame(
                Parameter = param,
                Mean = mean_value,
                Median = median_value,
                HPDI_Lower = hpdi[1, "lower"],
                HPDI_Upper = hpdi[1, "upper"],
                n_eff = n_eff,
                Rhat = rhat,
                stringsAsFactors = FALSE
            ))
            }
        }

        for (param in names(transformed_params)) {
            param_samples <- transformed_params[[param]]
            param_samples <- coda::as.mcmc(as.matrix(param_samples))

            mean_value <- mean(param_samples)
            median_value <- stats::median(param_samples)
            hpdi <- coda::HPDinterval(param_samples, prob = prob)
            n_eff <- coda::effectiveSize(param_samples)
            rhat <- rstan::Rhat(as.matrix(param_samples))

            summary_table <- rbind(summary_table, data.frame(
                Parameter = param,
                Mean = mean_value,
                Median = median_value,
                HPDI_Lower = hpdi[1, "lower"],
                HPDI_Upper = hpdi[1, "upper"],
                n_eff = n_eff,
                Rhat = rhat,
                stringsAsFactors = FALSE
            ))
    }
    row.names(summary_table) = NULL
    numeric_cols <- sapply(summary_table, is.numeric)
    summary_table[numeric_cols] <- lapply(summary_table[numeric_cols], round, digits = digits)
    return(summary_table)
}

