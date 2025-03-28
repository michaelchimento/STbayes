#' Title
#'
#' @param data_list list object exported from import_user_STb or import_NBDA_STb
#' @param model_obj Can be either model object exported from generate_STb_model or filename
#' @param chains Integer Number of chains to run (default=1)
#' @param cores Integer Number of cores to use (default=1)
#' @param iter Integer Number of iterations to run for
#' @param control List of arguments to pass to control.
#' @param algorithm Defaults to "NUTS". If running asocial OADA, specify "Fixed_param".
#'
#' @return rstan fit
#' @importFrom rstan stan
#' @export
#'
#' @examples
#' data_list_user = import_user_STb(event_data, edge_list)
#' model_obj = generate_STb_model(data_list_user)
#' fit = fit_STb(data_list_user, model_obj, chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.99))
#' #or alternatively
#' fit = fit_STb(data_list_user, "path/to/your/model.stan", chains = 5, cores = 5, iter=2000, control = list(adapt_delta=0.99))
fit_STb <- function(data_list, model_obj, chains=1, cores=1, iter=1000, control=list(), algorithm="NUTS"){

        # if it's a filename
        if (substr(model_obj, nchar(model_obj) - 4, nchar(model_obj)) == ".stan") {
            temp_file = readLines(model_obj, warn = FALSE)
            model_obj = glue::glue_collapse(temp_file, sep = "\n")
        }

        #calculate nveff numberz
        N_veff = return_N_veff(model_obj)
        message(paste("Detected N_veff =", N_veff))
        data_list$N_veff = N_veff

        # check if the model is a character or a file path
        if (is.character(model_obj) && grepl("data \\{", model_obj)) {
            # Write the model code to a temporary file
            temp_file <- tempfile(fileext = ".stan")
            writeLines(model_obj, temp_file)
            model_obj <- temp_file
        }

        message("Compiling model...")
        # fit model
        model <- rstan::stan(
            file = model_obj,
            data = data_list,
            chains = chains,
            cores = cores,
            iter = iter,
            control = control,
            algorithm = algorithm
        )

        return(model)
}
