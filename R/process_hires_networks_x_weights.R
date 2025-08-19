utils::globalVariables(c(
    ":=", ".SD", ".",
    "t_weight", "trial", "time", "duration",
    "focal", "other", "discrete_time"
))
#' process_networks_x_weights_hires()
#'
#' Helper function that pre-processes networks cxn x t_weights for fast high-res models.
#'
#' @param event_data event_data df
#' @param t_weights t_weights df
#' @param networks networks df
#' @param D_data D_data df
#'
#' @importFrom data.table as.data.table setkey merge.data.table
#' @return dataframe where network weights are pre-processed to include cxn x t_weight

process_networks_x_weights_hires <- function(event_data, t_weights, networks, D_data) {
    t_weights$trial <- as.integer(as.factor(t_weights$trial))
    t_weights$id <- as.integer(t_weights$id_numeric)

    names(t_weights)[names(t_weights) == "id"] <- "other"
    names(D_data)[names(D_data) == "trial_numeric"] <- "trial"
    # join with networks

    message("|Hi-Rez|: Pre-processing for standard transmission.")
    networks <- data.table::as.data.table(networks)
    t_weights <- data.table::as.data.table(t_weights)
    temp_net <- data.table::merge.data.table(networks, t_weights, by = c("trial", "time", "other"), all.x = TRUE)
    weight_cols <- setdiff(names(temp_net), c("trial", "time", "other", "focal", "t_weight", "id_numeric"))
    temp_net[, (weight_cols) := lapply(.SD, function(x) x * t_weight), .SDcols = weight_cols]

    temp_net$t_weight <- NULL

    # do the rolling joint
    network_dt <- temp_net
    cutpoints_dt <- data.table::as.data.table(D_data)
    data.table::setkey(cutpoints_dt, trial, time)

    # roll = -Inf ensures that each row gets the most recent discrete time value from D_data
    network_dt <- cutpoints_dt[network_dt, on = c("trial", "time"), roll = -Inf]
    # aggregate, .SDcols = specifies which columns
    network_dt <- network_dt[, c(.(
        duration = unique(duration)
    ), lapply(.SD, sum)),
    by = .(trial, focal, other, discrete_time),
    .SDcols = weight_cols
    ]
    # network_dt$weight = network_dt$weight/(network_dt$duration)
    network_dt[, (weight_cols) := lapply(.SD, function(x) x / duration), .SDcols = weight_cols]

    # remove extra column
    network_dt$time <- network_dt$discrete_time
    network_dt$discrete_time <- NULL
    network_dt$duration <- NULL
    message("|Hi-Rez|: Pre-processing completed.")
    return(as.data.frame(network_dt))
}
