#' Convenience to extract WAIC scores
#'
#' @param fit A STb model fit
#'
#' @return WAIC estimate
#' @export
#'
#' @examples
extract_WAIC <- function(fit) {
    if (inherits(fit, "CmdStanMCMC")) {
        ll <- fit$draws("log_lik", format = "draws_matrix")
    } else {
        stop("Please provide CmdStanMCMC object, rather than ", class(fit)[1])
    }
    WAIC = loo::waic(ll)
    return(WAIC$estimates["waic", "Estimate"])
}
#' Convenience to extract LOOIC scores
#'
#' @param fit A STb model fit
#'
#' @return LOOIC estimate
#' @export
#'
#' @examples
extract_LOOIC <- function(fit) {
    if (inherits(fit, "CmdStanMCMC")) {
        ll <- fit$draws("log_lik", format = "draws_matrix")
    } else {
        stop("please provide CmdStanMCMC object, rather than ", class(fit)[1])
    }
    loo_obj <- loo::loo(ll)
    return(loo_obj$estimates["looic", "Estimate"])
}
