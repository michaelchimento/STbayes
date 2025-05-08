#' helper function for importing ILVs
#'
#' @param x column or vector of data
#'
#' @return
#'
#' @examples
detect_ILV_datatype <- function(x) {
  x <- stats::na.omit(x)

  if (is.logical(x) || (is.numeric(x) && all(x %in% c(0, 1)))) {
    return("binary")
  }

  # multi-level categorical: integer-like with > 2 unique values
  if (is.numeric(x) && all(x == floor(x)) && length(unique(x)) > 2) {
    return("categorical")
  }

  if (is.numeric(x)) {
    return("continuous")
  }

  # catch-all for factors or characters
  return("categorical")
}
