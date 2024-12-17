#' Convenience to extract WAIC scores
#'
#' @param fit A STb model fit
#'
#' @return Table with WAIC estimates
#' @importFrom rstan extract
#' @importFrom loo waic
#' @export
#'
#' @examples
extract_WAIC <- function(fit) {
    ll = rstan::extract(fit, pars = "log_lik")
    WAIC = loo::waic(ll[["log_lik"]])
    return(WAIC)
}