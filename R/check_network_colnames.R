#' check_network_colnames()
#'
#' Helper function that makes sure correct colnames are used. Will be removed eventually.
#'
#' @param networks networks dataframe
#'
#' @returns networks dataframe
check_network_colnames <- function(networks) {
    if (all(c("from", "to") %in% names(networks))) {
        warning("\u26A0\ufe0f Columns `from` and `to` are deprecated. Renaming to `focal` and `other`.")
        names(networks)[names(networks) == "from"] <- "focal"
        names(networks)[names(networks) == "to"] <- "other"
    }
    return(networks)
}
