#' Fit STb model using cmdstanr
#'
#' @param data_list list object exported from import_user_STb or import_NBDA_STb
#' @param model_obj Can be either model object exported from generate_STb_model or a filename
#' @param ... Additional arguments passed to mod$sample (e.g., chains, iter_warmup, iter_sampling, adapt_delta)
#'
#' @return A CmdStanMCMC fit object
#' @export
fit_STb <- function(data_list, model_obj, ...) {
  extra_args <- list(...)

  # valid args for $sample() method in cmdstanr
  allowed_sample_args <- c(
    "data", "seed", "chains", "parallel_chains", "threads_per_chain",
    "iter_warmup", "iter_sampling", "save_warmup", "thin", "max_treedepth",
    "adapt_engaged", "adapt_delta", "step_size", "metric", "init", "refresh",
    "sig_figs", "diagnostic_file", "profile_file", "show_messages",
    "opencl_ids", "fixed_param", "output_dir", "validate_csv"
  )
  valid_args <- extra_args[names(extra_args) %in% allowed_sample_args]

  if (!"init" %in% names(valid_args)) {
    if (data_list$multinetwork_s == "separate") {
      valid_args$init <- function(chain_id) list(log_lambda_0_mean = -4, log_s_prime_mean = rep(-4, data_list$N_networks))
    } else {
      valid_args$init <- function(chain_id) list(log_lambda_0_mean = -4, log_s_prime_mean = -4)
    }
  }


  if (!"iter_warmup" %in% names(valid_args) && !"iter_sampling" %in% names(valid_args) && "iter" %in% names(extra_args)) {
    iter <- extra_args$iter
    valid_args$iter_warmup <- floor(iter / 2)
    valid_args$iter_sampling <- ceiling(iter / 2)
  }

  # deal w inline stan code or file
  if (is.character(model_obj) && grepl("\\.stan$", model_obj)) {
    model_code <- readLines(model_obj, warn = FALSE)
    model_code <- glue::glue_collapse(model_code, sep = "\n")
  } else if (is.character(model_obj) && grepl("data \\{", model_obj)) {
    model_code <- model_obj
  } else {
    stop("model_obj must be a file path to an existing .stan file or Stan code as a string (use generate_STb_model()).")
  }

  # Write model to temp file
  temp_file <- tempfile(fileext = ".stan")
  writeLines(model_code, temp_file)

  # Add N_veff
  N_veff <- return_N_veff(model_code)
  message(paste("Detected N_veff =", N_veff))
  data_list$N_veff <- N_veff

  mod <- cmdstanr::cmdstan_model(temp_file)

  message("Sampling...")

  # remove elements of datalist not used... why cmd stan...
  data_list_clean <- Filter(function(x) {
    is.numeric(x) || is.integer(x) || is.logical(x) || is.array(x) || is.matrix(x)
  }, data_list)

  valid_args$data <- data_list_clean

  fit <- do.call(mod$sample, valid_args)

  message("🫴 Use STb_save() to save both the fit and chain csvs in a single, convenient RDS file. Or don't!")

  return(fit)
}
