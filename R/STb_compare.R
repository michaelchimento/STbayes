#' STb_compare() automates the workflow of loo-psis or waic elpd comparisons. Relies on loo to do this.
#'
#' @param ... CmdStanMCMC model fits for comparison
#' @param model_names an optional list of model names, otherwise taken from object names
#' @param method a string either "loo-psis" or "waic" to indicate the method used for elpd. Defaults to "loo-psis". Pareto diagnostics are not calculated for WAIC.
#'
#' @return list containing loo_objects, comparison, and pareto_diagnostics if using loo-psis.
#' @export
STb_compare <- function(..., model_names = NULL, method = "loo-psis") {
  args <- list(...)
  if (is.list(args[[1]])) {
    models <- args[[1]]
    model_names <- names(models)
  } else {
    models <- list(...)
    if (is.null(model_names)) {
      model_names <- as.character(match.call(expand.dots = FALSE)$`...`)
    }
  }

  if (length(models) < 2) {
    stop("Provide at least two fitted models.")
  }

  if (length(models) < 2) {
    stop("Provide at least two fitted models.")
  }

  if (!method %in% c("loo-psis", "waic")) {
    stop("Invalid method. Choose 'loo-psis', 'waic'.")
  }

  if (is.null(model_names)) {
    model_names <- as.character(match.call(expand.dots = FALSE)$`...`)
  }

  message(sprintf("Calculating %s.", toupper(method)))

  loo_results <- mapply(function(model, name) {
    if (!inherits(model, c("CmdStanMCMC"))) {
      stop(sprintf("Model '%s' must be of class CmdStanMCMC.", name))
    }
    log_lik <- model$draws("log_lik", format = "draws_matrix")
    log_lik <- log_lik[, colSums(log_lik) != 0, drop = FALSE] # drop cols that are zero cuz it's just demos

    if (method == "loo-psis") {
      return(loo::loo(log_lik))
    } else {
      return(loo::waic(log_lik))
    }
  }, models, model_names, SIMPLIFY = FALSE)

  names(loo_results) <- model_names

  message("Comparing models.")
  comparison <- loo::loo_compare(loo_results)

  # Pareto k diagnostics (only for loo-psis)
  if (method == "loo-psis") {
    message("Calculating pareto-k diagnostic (only for loo-psis).")

    pareto_diagnostics <- lapply(seq_along(loo_results), function(i) {
      k_values <- loo_results[[i]]$diagnostics$pareto_k
      problem_indices <- which(k_values > 0.7)

      if (length(problem_indices) > 0) {
        cat(sprintf(
          "Model '%s' has %d problematic Pareto k values (k > 0.7).",
          model_names[i], length(problem_indices)
        ))
        if (any(k_values > 1)) {
          cat(sprintf(
            "  - WARNING: %d observations have k > 1 (LOO is unreliable).\n",
            sum(k_values > 1)
          ))
        }
        cat("  - Problematic observation indices:", paste(problem_indices, collapse = ", "), "\n")
      }

      data.frame(
        model = model_names[i],
        observation = seq_along(k_values),
        pareto_k = k_values
      )
    })

    pareto_diagnostics_df <- do.call(rbind, pareto_diagnostics)
  } else {
    pareto_diagnostics_df <- NULL
  }

  return(list(
    method = method,
    loo_objects = loo_results,
    comparison = comparison,
    pareto_diagnostics = pareto_diagnostics_df
  ))
}
