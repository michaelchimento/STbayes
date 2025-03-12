#' STb_compare()
#'
#' A function that automates the workflow of loo-psis or waic elpd comparisons. Relies on loo to do this.
#'
#' @param models a list containing model fits for comparison
#' @param model_names an optional list of model names, otherwise taken from object names
#' @param method a string either "loo-psis" or "waic" to indicate the method used for elpd. Defaults to "loo-psis". Pareto diagnostics are not calculated for WAIC.
#'
#' @return list containing loo_objects, comparison, and pareto_diagnostics if using loo-psis.
#' @export
#'
#' @examples
STb_compare <- function(..., model_names = NULL, method = "loo-psis") {
  models <- list(...)

  # double check that at least two models are provided
  if (length(models) < 2) {
    stop("Provide at least two fitted models.")
  }

  # check that method is valid
  if (!method %in% c("loo-psis", "waic", "kfold")) {
      stop("Invalid method. Choose 'loo', 'waic'.")
  }

  # get names from objects if not provided
  if (is.null(model_names)) {
    model_names <- as.character(match.call(expand.dots = FALSE)$`...`)
  }

  # compute LOO-PSIS or WAIC for each model
  message(sprintf("Calculating %s.", toupper(method)))

  loo_results <- lapply(models, function(model) {
    if (!inherits(model, c("stanreg", "brmsfit", "stanfit"))) {
      stop("Each model must be a fitted Bayesian model.")
    }

    if (method == "loo-psis") {
      return(loo::loo(model))
    } else {
        ll = rstan::extract(model, pars = "log_lik")
        WAIC = loo::waic(ll[["log_lik"]])
      return(WAIC)
    }
  })

  # keep identifiable
  names(loo_results) <- model_names

  message("Comparing models.")
  comparison <- loo::loo_compare(loo_results)

  # Extract Pareto k diagnostics (only for LOO)
  if (method == "loo-psis") {
    message("Calculating pareto-k diagnostic (only for loo-psis).")

    pareto_diagnostics <- lapply(seq_along(loo_results), function(i) {
      k_values <- loo_results[[i]]$diagnostics$pareto_k
      problem_indices <- which(k_values > 0.7) # Identify problematic points

      if (length(problem_indices) > 0) {
        warning(sprintf(
          "Model '%s' has %d problematic Pareto k values (k > 0.7).",
          model_names[i], length(problem_indices)
        ))

        if (any(k_values > 1)) {
          cat(sprintf(
            "  - WARNING: %d observations have k > 1, which means LOO is unreliable.\n",
            sum(k_values > 1)
          ))
        }
        cat("  - Problematic observation indices:", paste(problem_indices, collapse = ", "), "\n")
      }

      return(data.frame(
        model = model_names[i],
        observation = seq_along(k_values),
        pareto_k = k_values
      ))
    })

    # combine all Pareto k diagnostics
    pareto_diagnostics_df <- do.call(rbind, pareto_diagnostics)
  } else {
    # If WAIC, no Pareto k diagnostics
    pareto_diagnostics_df <- NULL
  }

  return(list(
    method = method,
    loo_objects = loo_results,
    comparison = comparison,
    pareto_diagnostics = pareto_diagnostics_df
  ))
}
