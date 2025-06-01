#' extract_WAIC()
#'
#' Convenience function to extract WAIC scores. You should really be using [STb_compare()].
#'
#' @param fit A STb model fit
#'
#' @return WAIC estimate
#' @importFrom loo waic
#' @export
extract_WAIC <- function(fit) {
    if (inherits(fit, "CmdStanMCMC")) {
        ll <- fit$draws("log_lik", format = "draws_matrix")
    } else {
        stop("Please provide CmdStanMCMC object, rather than ", class(fit)[1])
    }
    WAIC <- loo::waic(ll)
    return(WAIC$estimates["waic", "Estimate"])
}
