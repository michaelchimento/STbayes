#' Return summary table for CmdStanR fit
#'
#' @param fit CmdStanMCMC model fit
#' @param depth integer depth of multidimensional parameters to extract
#' @param prob double limits for HPD of estimates (default = 0.95)
#' @param ignore_params character vector of parameters to ignore
#' @param digits integer of digits to round to
#'
#' @return Summary table
#' @importFrom posterior as_draws_df summarise_draws rhat ess_bulk
#' @export
STb_summary <- function(fit, depth = 1, prob = 0.95,
                        ignore_params = c("lp__", "idx", "log_lik", "log_lik_matrix",
                                          "acquisition_time", "z_ID", "Rho_ID", "v_ID", ".chain", ".iteration", ".draw"),
                        digits = 3) {

    if (!inherits(fit, c("CmdStanMCMC"))) {
        stop(sprintf("Model '%s' must be a CmdStanMCMC model fit.", name))
    }

    draws_df <- fit$draws(format = "draws_df")

    # Get parameter names and filter
    param_names <- names(draws_df)
    sigma_ID_params <- grep("^sigma_ID\\[", param_names, value = TRUE)
    keep_params <- param_names[!param_names %in% c(ignore_params)]

    # Optionally filter by depth (unchanged)
    depth_filter <- function(param) {
        parts <- strsplit(param, "\\[")[[1]]
        return(length(gsub("]", "", parts)) <= depth)
    }
    keep_params <- Filter(depth_filter, keep_params)
    keep_params <- union(keep_params, sigma_ID_params)
    # Subset
    draws_df <- posterior::subset_draws(draws_df, variable = keep_params)

    # Summarize draws
    summary_stats <- posterior::summarise_draws(
        draws_df,
        ~ stats::median(.x),
        ~ posterior::mad(.x),
        ~ posterior::quantile2(.x, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2)),
        ~ posterior::ess_bulk(.x),
        ~ posterior::rhat(.x)
    )

    # Clean column names and round
    names(summary_stats) <- c("Parameter", "Median", "MAD", "HPDI_Lower", "HPDI_Upper", "n_eff", "Rhat")
    numeric_cols <- sapply(summary_stats, is.numeric)
    summary_stats = as.data.frame(summary_stats)
    summary_stats[numeric_cols] <- lapply(summary_stats[numeric_cols], round, digits = digits)

    return(summary_stats)
}
