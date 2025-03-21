#' Helper function to transform transmission weights into same size as C
#'
#' @param t_weights dataframe with t_weight. all ids across all trials must be present.
#'
#' @return W[k,t,n]
#'
#' @examples
#' t_weights <- data.frame(
#' trial = c(rep(1, each = 9), rep(2, each = 9)),
#' id = c(rep(LETTERS[1:3], each = 3), rep(LETTERS[4:6], each = 3)),
#' time = c(rep(1:3, times = 3), rep(1:3, times = 3)),
#' t_weight = exp(rnorm(18))
#' )
#' create_tw_matrix(t_weights)
create_W_matrix <- function(t_weights){
    t_weights <- t_weights[order(t_weights$trial, t_weights$time), ]
    W <- unclass(with(t_weights, xtabs(t_weight ~ trial + time + id)))
    dimnames(W) <- NULL  # remove dimension names
    return(W)
}


