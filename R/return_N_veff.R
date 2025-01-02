#' return_N_veff: utility function to find the number of parameters with varying effects
#'
#' @param text model object (text)
#'
#' @return integer of number of varying effects present in model
#'
#' @examples
return_N_veff <- function(text) {
    # regex match v_ID[,*]
    matches <- regmatches(text, gregexpr("v_ID\\[,(\\d+)\\]", text, perl = TRUE))

    # extract values
    numbers <- as.numeric(unlist(regmatches(unlist(matches), gregexpr("\\d+", unlist(matches)))))

    # return the highest number, or 0 if no matches found
    if (length(numbers) > 0) {
        return(max(numbers))
    } else {
        return(0)
    }
}
