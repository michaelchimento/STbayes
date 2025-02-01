#' check_integer: private function to check whether or not number is integer or integer-like. can't believe there's not a base function for this!
#'
#' @param x number or vector
#'
#' @return boolean values of whether a number is integer or integer-like
check_integer <- function(x) {
    x == round(x)
}
