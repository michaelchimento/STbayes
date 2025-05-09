#' order_params
#'
#' @param p
#'
#' @return
#'
#' @examples
order_params <- function(p) {
    if (grepl("^log_lambda_0", p)) return(1)
    if (grepl("^log_s_prime_mean", p)) return(2)
    if (grepl("^lambda_0_mean", p)) return(3)
    if (grepl("^s_mean", p)) return(4)
    if (grepl("^beta_", p)) return(5)
    if (grepl("^sigma_ID", p)) return(6)
    if (grepl("^percent_ST", p)) return(7)
    return(8)  # everything else
}
