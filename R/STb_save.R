#' Save STbayes objects easily
#'
#' @param fit model fitted by STbayes
#' @param output_dir where do you want to save it???
#' @param name optional name to use for model
#'
#' @return
#' @export
#'
#' @examples
STb_save <- function(fit, output_dir = "cmdstan_saves", name = NULL) {
  # infer name if not provided
  if (is.null(name)) {
    name <- deparse(substitute(fit))
  }
  # just get rid for nice printing later
  output_dir <- sub("/+$", "", output_dir)

  full_path <- file.path(output_dir)

  # create dir
  dir.create(full_path, showWarnings = FALSE, recursive = TRUE)

  # generate CSV file names
  csv_files <- file.path(full_path, paste0(name, "_chain_", seq_along(fit$output_files()), ".csv"))

  # update the fit object to use new file paths
  fit_newfiles <- cmdstanr::as_cmdstan_fit(fit$output_files())

  # save the R object
  saveRDS(fit_newfiles, file = file.path(full_path, paste0(name, ".rds")))

  message(paste0("Fit & chains successfully saved 💾"))
  message(paste0("You can load fit again with 👉 readRDS('", full_path, "/", name, ".rds')"))
}
