#' generate_STb_model: Dynamically generate STAN model based on input data
#' @param STb_data a list of formatted data returned from the STbayes_data() function
#' @param model_type string specifying the model type: "asocial" or "full"
#' @param data_type string specifying the type of data you have ("time" or "order"). "time" will generate cTADA specification, "order" will generate OADA specification.
#' @param transmission_func string specifying transmission function: "standard", "freq-dep" for frequency dependent complex contagion. Defaults to "standard".
#' @param veff_ID Parameters for which to estimate varying effects by individuals. Default is no varying effects.
#' @param gq Boolean to indicate whether the generated quantities block is added (incl. ll for WAIC)
#' @param est_acqTime Boolean to indicate whether gq block includes estimates for acquisition time. At the moment this uses 'one weird trick' to accomplish this and does not support estimates for non-integer learning times.
#' @param priors named list with strings containing the prior for log baserate, s, f, k. defaults to list(log_lambda_0 = "uniform(-4, 3)", log_s = "uniform(-4, 3)", log_f = "normal(0,1)", k_raw = "normal(0,1)")
#' @export
#' @return A STAN model (character) that is customized to the input data.
#'
#' @examples
#' #very mock data
#' event_data <- data.frame(
#'   id = c("A", "B", "C", "D", "E", "F"),
#'   trial = c(1, 1, 1, 2, 2, 2),
#'   time = c(0, 1, 2, 0, 1, 4),
#'   max_time = c(3, 3, 3, 4, 4, 4)
#' )
#' networks <- data.frame(
#'   trial = c(1, 1, 1, 2, 2, 2),
#'   from = c("A", "A", "B", "D", "D", "E"),
#'   to = c("B", "C", "C", "E", "F", "F"),
#'   kin = c(1, 0, 1, 0, 1, 1),
#'   inverse_distance = c(0, 1, .5, .25, .1, 0)
#' )
#'  ILV_c<- data.frame(
#'     id = LETTERS[1:6],
#'     age = c(-1, -2, 0, 1, 2), # continuous variables should be normalized
#'     sex = c(0, 1, 1, 0, 1, 0), # Factor ILVs must be input as numeric
#'     weight = c(0.5, .25, .3, 0, -.2, -.4)
#'  )
#' data_list <- import_user_STb(
#'   event_data = event_data,
#'   networks = networks,
#'   ILV_c = ILV_c,
#'   ILVi = c("age"), # Use only 'age' for asocial learning
#'   ILVs = c("sex"), # Use only 'sex' for social learning
#'   ILVm = c("weight") # Use weight for multiplicative effect on asocial and social learning
#' )
#'
#' model = generate_STb_model(data_list) # creates full specification of cTADA model, no varying effects and default priors.
#' model = generate_STb_model(data_list, veff_ID = c("lambda_0", "s")) # estimate varying effects by ID for baseline learning rate and strength of social transmission.
#' model = generate_STb_model(data_list, data="order") # creates OADA specification
#' model = generate_STb_model(data_list, priors=list(log_lambda_0 = "normal(7, 3)", log_s = "uniform(7, 3)")) # adjust priors

#' print(model)
generate_STb_model <- function(STb_data,
                               data_type = c("time", "order"),
                               model_type = c("full","asocial"),
                               transmission_func = c("standard","freq-dep1","freq-dep2"),
                               veff_ID = c(),
                               gq = TRUE,
                               est_acqTime = FALSE,
                               priors = list()) {

    data_type <- match.arg(data_type)
    model_type <- match.arg(model_type)
    transmission_func <- match.arg(transmission_func)

    if (data_type == "order" && "lambda_0" %in% veff_ID) {
        stop('lambda_0 cannot be veff_ID when using OADA.')
    }

    if (transmission_func != "standard" && length(STb_data$network_names)>1) {
        stop('Complex transmission can only be used with single network models. Please use transmission_func="standard"')
    }

    default_priors <- list(log_lambda_0 = "normal(-4, 3)",
                           log_s = "normal(-4, 3)",
                           log_f = "normal(0,1)",
                           k_raw = "normal(0,1)")

    priors <- utils::modifyList(default_priors, priors)

    if (data_type == "time") {
        return(generate_STb_model_cTADA(STb_data = STb_data,
                                        model_type = model_type,
                                        transmission_func = transmission_func,
                                        veff_ID = veff_ID,
                                        gq = gq,
                                        est_acqTime = est_acqTime,
                                        priors = priors))
    } else {
        return(generate_STb_model_OADA(STb_data = STb_data,
                                       model_type = model_type,
                                       transmission_func = transmission_func,
                                       veff_ID = veff_ID,
                                       gq = gq,
                                       priors = priors))
    }
}
