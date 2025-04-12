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
#' create_W_matrix(t_weights)
create_W_matrix <- function(t_weights, max_time){

    t_weights$trial_numeric <- as.numeric(as.factor(t_weights$trial))
    t_weights$id_numeric <- as.numeric(as.factor(t_weights$id))
    if (!"time" %in% names(t_weights)) {
        # Static case: broadcast static t_weight across time
        message("Detected static transmission weights: broadcasting across time.")

        trial_levels <- sort(unique(t_weights$trial_numeric))
        id_levels <- sort(unique(t_weights$id_numeric))
        max_times <- setNames(rep(max_time, length(trial_levels)), trial_levels)

        # Create full grid of trial x time x id
        rows <- do.call(rbind, lapply(trial_levels, function(tr) {
            do.call(rbind, lapply(1:max_times[as.character(tr)], function(t) {
                cbind(trial_numeric = tr,
                      discrete_time = t,
                      id_numeric = id_levels)
            }))
        }))

        rows <- as.data.frame(rows)
        t_weights_static <- merge(rows, t_weights, by = c("trial_numeric", "id_numeric"), all.x = TRUE)
        W <- with(t_weights_static, xtabs(t_weight ~ trial_numeric + discrete_time + id_numeric))

    }
    else{
        #add discrete time in case user supplied cont. time of events.
        t_weights$discrete_time <- NA
        t_weights$discrete_time[t_weights$time != 0] <- with( # assign values only where time != 0
            t_weights[t_weights$time != 0, ],
            ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x)))
        )
        t_weights$discrete_time[t_weights$time == 0] <- 0
        t_weights <- t_weights[order(t_weights$trial, t_weights$discrete_time), ]
        W <- unclass(with(t_weights, xtabs(t_weight ~ trial_numeric + discrete_time + id_numeric)))
    }

    dimnames(W) <- NULL  # remove dimension names
    return(W)
}


