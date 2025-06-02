#' STb_save()
#'
#' Save STbayes fits easily, rather than having to save the fit and chain csvs separately.
#' Will provide a message with the command to re-load the fit for convenience.
#'
#' @param fit model fitted by STbayes
#' @param output_dir path of where do you want to save
#' @param name optional name to use for model
#'
#' @return NULL, message if successfully saved.
#' @export
#'
#' @examples
#' \dontrun{
#' data_list <- import_user_STb(STbayes::event_data, STbayes::edge_list)
#' model_obj <- generate_STb_model(data_list, gq = TRUE)
#' fit <- fit_STb(data_list,
#'     model_obj,
#'     parallel_chains = 4,
#'     chains = 4,
#'     cores = 4,
#'     iter = 4000,
#'     refresh = 2000
#' )
#' STb_save(fit, output_dir = "../data/stan_fits", name = "my_fit")
#' }
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

    message(paste0("Fit & chains successfully saved \U0001F4BE"))
    message(paste0("You can load fit again with \U0001F449 readRDS('", full_path, "/", name, ".rds')"))
}
