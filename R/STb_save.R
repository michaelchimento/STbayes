#' Save STbayes objects easily
#'
#' @param fit model fitted by STbayes
#' @param name optional name to use for model
#' @param path where do you want to save it???
#'
#' @return
#' @export
#'
#' @examples
STb_save <- function(fit, name = NULL, path = "cmdstan_saves") {
    # infer  name if not provided
    if (is.null(name)) {
        name <- deparse(substitute(fit))
    }

    # create dir
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    # generate CSV file names
    csv_files <- file.path(path, paste0(name, "_chain_", seq_along(fit$output_files()), ".csv"))

    # copy CSVs
    file.copy(fit$output_files(), to = csv_files, overwrite = TRUE)

    # update the fit object to use new file paths
    fit$run_set$output_files_ <- csv_files

    # save the R object
    saveRDS(fit, file = file.path(path, paste0(name, ".rds")))

    message(paste("Fit + chains successfully saved in", path))
}
