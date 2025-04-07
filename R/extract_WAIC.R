#' Convenience to extract WAIC scores
#'
#' @param fit A STb model fit
#'
#' @return Table with WAIC estimates
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
