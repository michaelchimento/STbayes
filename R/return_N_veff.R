#' return_N_veff()
#'
#' Helper function to find the number of parameters with varying effects
#'
#' @param text model object (text)
#' @param num_networks integer, the number of networks (data value of N_networks).
#'   Loop-variable index forms such as \code{v_id[,n]} span this many columns.
#'   Defaults to 1, which leaves literal-digit models unchanged.
#'
#' @return integer of number of varying effects present in model
return_N_veff <- function(text, num_networks = 1) {
    # capture the column (second) index expression of every v_id[...] / v_trial[...]
    # reference: skip the optional first index and comma, then grab up to the ]
    pattern <- "v_(?:id|trial)\\[\\s*[^,\\]]*,\\s*([^\\]]+?)\\s*\\]"
    m <- regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1]]
    exprs <- sub(pattern, "\\1", m, perl = TRUE)

    # evaluate each captured column expression to a max-column value
    values <- vapply(exprs, function(e) {
        e <- trimws(e)
        if (grepl("^\\d+$", e)) {
            as.numeric(e)
        } else if (grepl("^n$", e)) {
            num_networks
        } else if (grepl("^n\\s*\\+\\s*(\\d+)\\s*-\\s*1$", e)) {
            d <- as.numeric(sub("^n\\s*\\+\\s*(\\d+)\\s*-\\s*1$", "\\1", e))
            num_networks + d - 1
        } else {
            NA_real_
        }
    }, numeric(1))

    values <- values[!is.na(values)]

    # return the highest value, or 0 if no matches found
    if (length(values) > 0) {
        return(max(values))
    } else {
        return(0)
    }
}
