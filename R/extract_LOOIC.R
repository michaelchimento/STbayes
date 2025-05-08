#' Convenience to extract LOOIC scores
#'
#' @param fit A STb model fit
#'
#' @return LOOIC estimate
#' @importFrom loo loo
#' @export
extract_LOOIC <- function(fit) {
  if (inherits(fit, "CmdStanMCMC")) {
    ll <- fit$draws("log_lik", format = "draws_matrix")
  } else {
    stop("please provide CmdStanMCMC object, rather than ", class(fit)[1])
  }
  loo_obj <- loo::loo(ll)
  return(loo_obj$estimates["looic", "Estimate"])
}
