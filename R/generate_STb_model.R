#' generate_STb_model()
#'
#' Second step of analysis pipeline: dynamically generates a STAN model based on
#'  input data. After saving the output to a variable, you can preview a
#'  formatted version in the R console using cat().
#'
#' @param STb_data a list of formatted data returned from the STbayes_data() function
#' @param data_type string specifying the type of data you have ("continuous_time" for cTADA, "discrete_time" for dTADA or "order" for OADA). continuous_time assumes you know precisely when events happened. discrete_time assumes you know roughly when individuals learned (within some discrete period), and order assumes that you have no time information.
#' @param model_type string specifying the model type: "full" or "asocial"
#' @param intrinsic_rate Define shape of intrinsic rate (either "constant" or "weibull"). Weibull fits extra parameter (gamma) that allows for time-varying event rates.
#' @param transmission_func string specifying transmission function: "standard", "freqdep_f" or "freqdep_k" for frequency dependent complex contagion. Defaults to "standard".
#' @param veff_params Vector of parameter names (string) for which to estimate varying effects. Default is no varying effects.
#' @param veff_type string/vector specifying whether varying effects should be applied to "id", "trial" or both (c("id","trial")). Default is id, applied only if user also gives `veff_params`.
#' @param gq Boolean to indicate whether the generated quantities block is added (incl. ll for WAIC)
#' @param est_acqTime Boolean to indicate whether gq block includes estimates for acquisition time. At the moment this uses 'one weird trick' to accomplish this and does not support estimates for non-integer learning times.
#' @param priors named list with strings containing priors.
#' @export
#' @return A STAN model string.
#'
#' @examples
#' # very mock data
#' event_data <- data.frame(
#'     id = LETTERS[1:6],
#'     trial = c(1, 1, 1, 2, 2, 2),
#'     time = c(0, 1, 2, 0, 1, 4),
#'     t_end = c(3, 3, 3, 4, 4, 4)
#' )
#' networks <- data.frame(
#'     trial = c(1, 1, 1, 2, 2, 2),
#'     from = c("A", "A", "B", "D", "D", "E"),
#'     to = c("B", "C", "C", "E", "F", "F"),
#'     kin = c(1, 0, 1, 0, 1, 1),
#'     inverse_distance = c(0, 1, .5, .25, .1, 0)
#' )
#' ILV_c <- data.frame(
#'     id = LETTERS[1:6],
#'     age = c(-1, -2, 0, 1, 2, 3), # continuous variables should be normalized
#'     sex = c(0, 1, 1, 0, 1, 0), # Factor ILVs must be input as numeric
#'     weight = c(0.5, .25, .3, 0, -.2, -.4)
#' )
#' data_list <- import_user_STb(
#'     event_data = event_data,
#'     networks = networks,
#'     ILV_c = ILV_c,
#'     ILVi = c("age"), # Use only 'age' for asocial rate
#'     ILVs = c("sex"), # Use only 'sex' for social rate
#'     ILVm = c("weight") # Use weight for multiplicative effect on both
#' )
#' # creates full specification of cTADA model, no varying effects and default priors.
#' model <- generate_STb_model(data_list)
#' # estimate varying effects by id, trial or both for intrinsic and social rates.
#' model <- generate_STb_model(data_list, veff_params = c("lambda_0", "s"), veff_type = "id")
#' model <- generate_STb_model(data_list,
#'     veff_params = c("lambda_0", "s"),
#'     veff_type = "trial"
#' )
#' model <- generate_STb_model(data_list,
#'     veff_params = c("lambda_0", "s"),
#'     veff_type = c("id", "trial")
#' )
#' # creates OADA specification
#' model <- generate_STb_model(data_list, data = "order")
#' # adjust priors
#' model <- generate_STb_model(data_list, priors = list(
#'     log_lambda0 = "normal(-2, 3)",
#'     log_sprime = "uniform(-7, 3)"
#' ))
#' # quickly inspect model code
#' cat(model)
generate_STb_model <- function(STb_data,
                               data_type = c("continuous_time", "discrete_time", "order"),
                               model_type = c("full", "asocial"),
                               intrinsic_rate = c("constant", "weibull"),
                               transmission_func = c("standard", "freqdep_f", "freqdep_k"),
                               veff_params = c(),
                               veff_type = "id",
                               gq = TRUE,
                               est_acqTime = FALSE,
                               priors = list()) {
    data_type <- match.arg(data_type, choices = c("continuous_time", "discrete_time", "order"))
    model_type <- match.arg(model_type, choices = c("full", "asocial"))
    intrinsic_rate <- match.arg(intrinsic_rate, choices = c("constant", "weibull"))
    transmission_func <- match.arg(transmission_func, choices = c("standard", "freqdep_f", "freqdep_k"))
    veff_type <- check_veff_type(veff_type)
    stopifnot(is.logical(gq), length(gq) == 1)
    stopifnot(is.logical(est_acqTime), length(est_acqTime) == 1)

    if (data_type == "order" && any(c("lambda_0", "gamma") %in% veff_params)) {
        stop("Intrinsic rate parameters (lambda_0, gamma) cannot be a varying effect in an OADA-type model.")
    }

    if (model_type == "asocial" && any(c("s", "f", "k") %in% veff_params)) {
        stop("Social transmission rate parameters (s, f, k) cannot be varying effects in an asocial model.")
    }

    default_priors <- list(
        log_lambda0 = "normal(-4, 2)",
        log_sprime = "normal(-4, 2)",
        beta_ILV = "normal(0,1)",
        log_f = "normal(0,1)",
        k_raw = "normal(0,3)",
        z_veff = "normal(0,1)",
        sigma_veff = "normal(0,1)",
        rho_veff = "lkj_corr_cholesky(3)",
        gamma = "normal(0,1)"
    )

    priors <- utils::modifyList(default_priors, priors)


    if (data_type == "continuous_time") {
        message("Creating cTADA type model with the following default priors:")
        for (p in names(priors)) {
            message(paste(p, "~", priors[[p]]))
        }
        return(generate_STb_model_TADA(
            STb_data = STb_data,
            model_type = model_type,
            intrinsic_rate = intrinsic_rate,
            transmission_func = transmission_func,
            dTADA = F,
            veff_params = veff_params,
            veff_type = veff_type,
            gq = gq,
            est_acqTime = est_acqTime,
            priors = priors
        ))
    } else if (data_type == "discrete_time") {
        message("Creating dTADA type model with the following default priors:")
        for (p in names(priors)) {
            message(paste(p, "~", priors[[p]]))
        }
        return(generate_STb_model_TADA(
            STb_data = STb_data,
            model_type = model_type,
            intrinsic_rate = intrinsic_rate,
            transmission_func = transmission_func,
            dTADA = T,
            veff_params = veff_params,
            veff_type = veff_type,
            gq = gq,
            est_acqTime = est_acqTime,
            priors = priors
        ))
    } else if (data_type == "order") {
        message("Creating OADA type model with the following default priors:")
        for (p in names(priors)) {
            message(paste(p, "~", priors[[p]]))
        }
        return(generate_STb_model_OADA(
            STb_data = STb_data,
            model_type = model_type,
            transmission_func = transmission_func,
            veff_params = veff_params,
            veff_type = veff_type,
            gq = gq,
            priors = priors
        ))
    } else {
        stop("Please choose data_type in 'continuous_time', 'discrete_time', 'order'.")
    }
}
