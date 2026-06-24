#' Order params
#'
#' Helper function for ordering params for STb_summary()
#'
#' @param p name of parameter
#'
#' @return integer indicating order
order_params <- function(p) {
    if (grepl("^log_lambda_0", p)) {
        return(111)
    }
    if (grepl("^log_s_prime_mean", p)) {
        return(112)
    }
    if (grepl("^lambda_0_mean", p)) {
        return(113)
    }
    if (grepl("^lambda_0", p)) {
        return(114)
    }
    if (grepl("^s_mean$", p)) {
        return(115)
    }
    if (grepl("^s_mean\\[", p)) {
        return(115)
    }
    if (grepl("^s$", p)) {
        return(116)
    } # scalar s
    if (grepl("^s\\[", p)) {
        return(117)
    } # vector s[i]
    if (grepl("^beta_", p)) {
        return(118)
    }
    if (grepl("^percent_ST", p)) {
        return(119)
    }
    if (grepl("^sigma_id", p)) {
        return(120)
    }
    if (grepl("^sigma_trial", p)) {
        return(121)
    }
    return(122)
}
