#' networks_dynamic_network_ex
#'
#' @docType data
#' @usage data(networks_dynamic_network_ex)
#' @format A data frame with 630 rows and 6 columns:
#' \\describe{
#' \\item{trial}{integer denoting trial (1 trial)}
#' \\item{focal}{integers denoting identity of focal individual (receiving influence)}
#' \\item{other}{integers denoting identity of other individual (exerting influence on focal)}
#' \\item{dist}{binary network indicating whether individuals were in proximity.}
#' \\item{avf}{binary network indicating whether focal was observing other.}
#' \\item{time}{time of each network, can be contiguous integers, or the event cutpoints.}
#' }
#' @source {STbayes} R package
#' @examples
#' data(networks_dynamic_network_ex)
#' head(networks_dynamic_network_ex)
