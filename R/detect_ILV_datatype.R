#' Detect ILV datatype
#'
#' Helper function for importing ILVs
#'
#' @param x column or vector of data
#'
#' @return string with category of variable
detect_ILV_datatype <- function(x) {
    x <- stats::na.omit(x)

    if (is.logical(x)) {
        return("boolean")
    } else if (is.character(x) || is.factor(x)) {
        return("categorical")
    } else if (is.numeric(x)) {
        return("continuous")
    } else {
        # catch-all
        stop("Please input boolean, continuous or categorical variables for ILVs. Categorical variables should be input as factors or characters.")
    }
}
