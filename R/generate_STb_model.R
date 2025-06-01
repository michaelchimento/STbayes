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
#' @param multinetwork_s string specifying how multi-network models are generated. "separate" estimates an s value for each network. "shared" generates model with single s and a vector of weights for each network.
#' @param veff_ID Parameters for which to estimate varying effects by individuals. Default is no varying effects.
#' @param gq Boolean to indicate whether the generated quantities block is added (incl. ll for WAIC)
#' @param est_acqTime Boolean to indicate whether gq block includes estimates for acquisition time. At the moment this uses 'one weird trick' to accomplish this and does not support estimates for non-integer learning times.
#' @param priors named list with strings containing the prior for log(lambda_0), log(s'), log(f), k, z_id, sigma_id, rho_id.
#' @export
#' @return A STAN model string.
#'
#' @examples
#' # very mock data
#' event_data <- data.frame(
#'     id = LETTERS[1:6],
#'     trial = c(1, 1, 1, 2, 2, 2),
#'     time = c(0, 1, 2, 0, 1, 4),
#'     max_time = c(3, 3, 3, 4, 4, 4)
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
#'     age = c(-1, -2, 0, 1, 2), # continuous variables should be normalized
#'     sex = c(0, 1, 1, 0, 1, 0), # Factor ILVs must be input as numeric
#'     weight = c(0.5, .25, .3, 0, -.2, -.4)
#' )
#' data_list <- import_user_STb(
#'     event_data = event_data,
#'     networks = networks,
#'     ILV_c = ILV_c,
#'     ILVi = c("age"), # Use only 'age' for asocial learning
#'     ILVs = c("sex"), # Use only 'sex' for social learning
#'     ILVm = c("weight") # Use weight for multiplicative effect on asocial and social learning
#' )
#'
#' model <- generate_STb_model(data_list) # creates full specification of cTADA model, no varying effects and default priors.
#' model <- generate_STb_model(data_list, veff_ID = c("lambda_0", "s")) # estimate varying effects by ID for baseline learning rate and strength of social transmission.
#' model <- generate_STb_model(data_list, data = "order") # creates OADA specification
#' model <- generate_STb_model(data_list, priors = list(log_lambda_0 = "normal(-2, 3)", log_sprime = "uniform(-7, 3)")) # adjust priors

#' print(model)
generate_STb_model <- function(STb_data,
                               data_type = c("continuous_time", "discrete_time", "order"),
                               model_type = c("full", "asocial"),
                               intrinsic_rate = c("constant", "weibull"),
                               transmission_func = c("standard", "freqdep_f", "freqdep_k"),
                               multinetwork_s = c("separate", "shared"),
                               veff_ID = c(),
                               gq = TRUE,
                               est_acqTime = FALSE,
                               priors = list(),
                               direct_s = F) {
    data_type <- match.arg(data_type)
    model_type <- match.arg(model_type)
    intrinsic_rate <- match.arg(intrinsic_rate)
    transmission_func <- match.arg(transmission_func)
    multinetwork_s <- match.arg(multinetwork_s)


    # used shared s config if only one network
    if (multinetwork_s == "separate" & STb_data$N_networks == 1) {
        multinetwork_s <- "shared"
    }
    STb_data$multinetwork_s <- multinetwork_s

    if (data_type == "order" && "lambda_0" %in% veff_ID) {
        stop("lambda_0 cannot be veff_ID when using OADA.")
    }

    default_priors <- list(
        log_lambda0 = "normal(-4, 2)",
        log_sprime = "normal(-4, 2)",
        beta_ILV = "normal(0,1)",
        log_f = "normal(0,1)",
        k_raw = "normal(0,3)",
        z_ID = "normal(0,1)",
        sigma_ID = "normal(0,1)",
        rho_ID = "lkj_corr_cholesky(3)",
        gamma = "normal(0,1)"
    )

    priors <- utils::modifyList(default_priors, priors)


    if (data_type == "continuous_time" & direct_s == F) {
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
            veff_ID = veff_ID,
            gq = gq,
            est_acqTime = est_acqTime,
            priors = priors
        ))
    } else if (data_type == "continuous_time" & direct_s == T) {
        priors <- utils::modifyList(priors, list(log_sprime = "normal(0,3)"))
        message("Creating cTADA type model (s modelled directly) with the following default priors:")
        for (p in names(priors)) {
            message(paste(p, "~", priors[[p]]))
        }
        return(generate_STb_model_TADA_sdirect(
            STb_data = STb_data,
            model_type = model_type,
            intrinsic_rate = intrinsic_rate,
            transmission_func = transmission_func,
            dTADA = F,
            veff_ID = veff_ID,
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
            veff_ID = veff_ID,
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
            veff_ID = veff_ID,
            gq = gq,
            priors = priors
        ))
    } else {
        stop("Please choose data_type in 'continuous_time', 'discrete_time', 'order'.")
    }
}
