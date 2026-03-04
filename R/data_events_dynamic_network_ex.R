#' events_dynamic_network_ex
#'
#' @docType data
#' @usage data(events_dynamic_network_ex)
#' @format A data frame with 10 rows and 4 columns:
#' \\describe{
#' \\item{id}{integers denoting identity}
#' \\item{trial}{integer denoting trial (1 trial)}
#' \\item{time}{time of event for each individual. time>t_end indicates censored individuals.}
#' \\item{t_end}{time of end of observation period.}
#' }
#' @source {STbayes} R package
#' @examples
#' data(events_dynamic_network_ex)
#' head(events_dynamic_network_ex)
