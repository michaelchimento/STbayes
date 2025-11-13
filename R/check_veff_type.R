#' check_veff_type
#'
#' @param veff_type string or character vector user entered
#'
#' @return formatted veffs, ofc
check_veff_type <- function(veff_type) {
    if (!is.character(veff_type)) {
        stop("`veff_type` must be a character vector")
    }

    # standardize
    veff_type <- tolower(veff_type)
    veff_type <- unique(veff_type)

    allowed <- c("id", "trial")

    # validate
    if (!all(veff_type %in% allowed)) {
        stop("`veff_type` must be one or both of: ", paste(allowed, collapse = ", "))
    }

    if (length(veff_type) > 2) {
        stop("`veff_type` must be of length 1 or 2")
    }

    sort(veff_type, decreasing = TRUE)
}
