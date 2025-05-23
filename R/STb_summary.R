#' Return summary table for CmdStanR fit
#'
#' @param fit CmdStanMCMC model fit
#' @param depth integer depth of multidimensional parameters to extract
#' @param prob double limits for HPD of estimates (default = 0.95)
#' @param ignore_params character vector of parameters to ignore
#' @param digits integer of digits to round to
#' @param CI_method "HPDI" for highest density interval or "PI" for quantiles (equal tails). Defaults to HPDI.
#'
#' @return Summary table
#' @importFrom posterior as_draws_df summarise_draws rhat ess_bulk ess_tail
#' @importFrom bayestestR hdi
#' @export
STb_summary <- function(fit, depth = 1, prob = 0.95,
                        ignore_params = c(
                          "lp__", "idx", "log_lik", "log_lik_matrix", "count_ST", "psocn_sum",
                          "acquisition_time", "z_ID", "Rho_ID", "v_ID", ".chain", ".iteration", ".draw", "s_prime"
                        ),
                        digits = 3,
                        CI_method = c("HPDI", "PI")) {
  if (!inherits(fit, c("CmdStanMCMC"))) {
    stop(sprintf("Model '%s' must be a CmdStanMCMC model fit.", name))
  }

  CI_method <- match.arg(CI_method)

  draws_df <- fit$draws(format = "draws_df")

  # Get parameter names and filter
  param_names <- names(draws_df)
  s_prime_params <- grep("^log_s_prime_mean\\[", param_names, value = TRUE)
  sigma_ID_params <- grep("^sigma_ID\\[", param_names, value = TRUE)
  s_mean_params <- grep("^s_mean\\[", param_names, value = TRUE)
  s_params <- grep("^s\\[", param_names, value = TRUE)
  ST_params <- grep("^percent_ST\\[", param_names, value = TRUE)
  keep_params <- param_names[!param_names %in% c(ignore_params)]

  # Optionally filter by depth (unchanged)
  depth_filter <- function(param) {
    parts <- strsplit(param, "\\[")[[1]]
    return(length(gsub("]", "", parts)) <= depth)
  }
  keep_params <- Filter(depth_filter, keep_params)
  keep_params <- Reduce(union, list(
    keep_params,
    s_prime_params,
    s_mean_params,
    s_params,
    sigma_ID_params,
    ST_params
  ))
  # Subset
  draws_df <- posterior::subset_draws(draws_df, variable = keep_params)

  # Summarize draws
  summary_stats <- posterior::summarise_draws(
    draws_df,
    ~ stats::median(.x),
    ~ posterior::mad(.x),
    ~ {
      if (CI_method == "HPDI") {
        hdi_vals <- bayestestR::hdi(as.numeric(.x), ci = prob)
        c(hdi_vals$CI_low, hdi_vals$CI_high)
      } else {
        posterior::quantile2(.x, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2))
      }
    },
    ~ posterior::ess_bulk(.x),
    ~ posterior::ess_tail(.x),
    ~ posterior::rhat(.x)
  )

  # Clean column names and round
  names(summary_stats) <- c("Parameter", "Median", "MAD", "CI_Lower", "CI_Upper", "ess_bulk", "ess_tail", "Rhat")
  numeric_cols <- sapply(summary_stats, is.numeric)
  summary_stats <- as.data.frame(summary_stats)
  summary_stats[numeric_cols] <- lapply(summary_stats[numeric_cols], round, digits = digits)

  #order
  summary_stats$order_priority <- vapply(summary_stats$Parameter, order_params, FUN.VALUE = numeric(1))
  #natural ordering of indexes
  # Extract numeric index from parameters (e.g., s[10] -> 10), scalar -> NA
  matches <- regmatches(summary_stats$Parameter, gregexpr("(?<=\\[)\\d+(?=\\])", summary_stats$Parameter, perl = TRUE))
  summary_stats$param_index <- sapply(matches, function(x) if (length(x) == 1) as.numeric(x) else NA_real_)
  summary_stats$param_index[is.na(summary_stats$param_index)] <- -1  # ensure scalar params sort first

  # Sort by group and then numeric index
  summary_stats <- summary_stats[order(summary_stats$order_priority, summary_stats$param_index), ]
  summary_stats$param_index <- NULL
  summary_stats$order_priority <- NULL

  return(summary_stats)
}
