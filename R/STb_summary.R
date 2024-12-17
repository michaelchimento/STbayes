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
STb_summary <- function(fit, depth = 1, prob = 0.95, ignore_params = c("s", "lambda_0", "lp__"), digits = 3) {
    # Extract posterior samples
    samples <- rstan::extract(fit, permuted = TRUE)

    # Get parameter names and depths
    param_names <- names(samples)
    param_depths <- sapply(strsplit(param_names, "\\["), function(x) length(x))

    # Filter parameters by depth and ignore specific parameters
    selected_params <- param_names[param_depths <= depth & !param_names %in% ignore_params]

    # Initialize summary table
    summary_table <- data.frame(
        Parameter = character(),
        Median = numeric(),
        HPDI_Lower = numeric(),
        HPDI_Upper = numeric(),
        n_eff = numeric(),
        Rhat = numeric(),
        stringsAsFactors = FALSE
    )

    # Process each parameter
        for (param in selected_params) {
            # Extract samples for the parameter
            param_samples <- samples[[param]]

            # Ensure param_samples is treated as a matrix for compatibility with HPDinterval
            param_samples <- coda::as.mcmc(as.matrix(param_samples))

            # Compute statistics
            mean_value <- mean(param_samples)
            median_value <- median(param_samples)
            hpdi <- coda::HPDinterval(param_samples, prob = prob)
            n_eff <- coda::effectiveSize(param_samples)
            rhat <- rstan::Rhat(as.matrix(param_samples))

            # Append to table
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

        # Add transformed parameters
        transformed_params <- list(
            transformed_s = exp(samples$log_s_mean),
            transformed_baserate = 1 / exp(samples$log_lambda_0_mean)
        )

        for (param in names(transformed_params)) {
            param_samples <- transformed_params[[param]]

            # Ensure param_samples is treated as a matrix for compatibility with HPDinterval
            param_samples <- coda::as.mcmc(as.matrix(param_samples))


        # Compute statistics
        mean_value <- mean(param_samples)
        median_value <- median(param_samples)
        hpdi <- coda::HPDinterval(param_samples, prob = prob)
        n_eff <- coda::effectiveSize(param_samples)
        rhat <- rstan::Rhat(as.matrix(param_samples))

        # Append to table
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
