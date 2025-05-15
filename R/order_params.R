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
    if (grepl("^lambda_0", p)) return(4)
    if (grepl("^s_mean$", p)) return(5)
    if (grepl("^s$", p)) return(6)         # scalar s
    if (grepl("^s\\[", p)) return(7)       # vector s[i]
    if (grepl("^sigma_ID", p)) return(8)
    if (grepl("^beta_", p)) return(9)
    if (grepl("^percent_ST", p)) return(10)
    return(11)
}
