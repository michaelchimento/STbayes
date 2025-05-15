#' Helper function for preprocessing "high resolution" data
#'
#' @param t_weights
#' @param networks
#' @param event_data
#'
#' @importFrom data.table as.data.table setkey merge.data.table
#' @return
#' @export
#'
#' @examples

process_networks_x_weights_hires <- function(event_data, t_weights, networks, D_data) {
    t_weights$trial = as.integer(as.factor(t_weights$trial))
    t_weights$id = as.integer(as.factor(t_weights$id))

    names(t_weights)[names(t_weights) == "id"] <- "to"
    names(D_data)[names(D_data) == "trial_numeric"] <- "trial"
    # join with networks

    message("|Hi-Rez|: Pre-processing for standard transmission.")
    networks <- data.table::as.data.table(networks)
    t_weights <- data.table::as.data.table(t_weights)
    temp_net <- data.table::merge.data.table(networks, t_weights, by = c("trial", "time", "to"), all.x = TRUE)
    weight_cols <- setdiff(names(temp_net), c("trial", "time", "to", "from", "t_weight"))
    temp_net[, (weight_cols) := lapply(.SD, function(x) x * t_weight), .SDcols = weight_cols]

    temp_net$t_weight <- NULL

    # do the rolling joint
    network_dt <- temp_net
    cutpoints_dt <- data.table::as.data.table(D_data)
    data.table::setkey(cutpoints_dt, trial, time)
    #network_dt$time <- ifelse(network_dt$time > 1, network_dt$time - 1, 1)

    # roll = -Inf ensures that each row gets the most recent discrete time value from D_data
    network_dt <- cutpoints_dt[network_dt, on = c("trial", "time"), roll = -Inf]
    # aggregate, .SDcols = specifies which columns
    network_dt <- network_dt[, c(.(
        duration = unique(duration)
    ), lapply(.SD, sum)),
    by = .(trial, from, to, discrete_time),
    .SDcols = weight_cols]
    network_dt$weight = network_dt$weight/(network_dt$duration)

    #remove extra column
    network_dt$time = network_dt$discrete_time
    network_dt$discrete_time = NULL
    network_dt$duration = NULL
    message("|Hi-Rez|: Pre-processing completed.")
    return(as.data.frame(network_dt))
}

# process_networks_x_weights_hires <- function(event_data, t_weights, networks, D_data) {
#     t_weights$trial = as.integer(as.factor(t_weights$trial))
#     t_weights$id = as.integer(as.factor(t_weights$id))
#
#     names(t_weights)[names(t_weights) == "id"] <- "to"
#     names(D_data)[names(D_data) == "trial_numeric"] <- "trial"
#     cutpoints_dt <- data.table::as.data.table(D_data)
#     # join with networks
#
#     message("|Hi-Rez|: Pre-processing for standard transmission.")
#     networks <- data.table::as.data.table(networks)
#     t_weights <- data.table::as.data.table(t_weights)
#     test <- cutpoints_dt[t_weights, on = c("trial", "time"), roll = -Inf]
#     test %>% group_by(discrete_time,to) %>% summarize(sum(t_weight))
#     temp_net <- data.table::merge.data.table(networks, t_weights, by = c("trial", "time", "to"), all.x = TRUE)
#     weight_cols <- setdiff(names(temp_net), c("trial", "time", "to", "from", "t_weight"))
#     temp_net[, (weight_cols) := lapply(.SD, function(x) x * t_weight), .SDcols = weight_cols]
#
#     temp_net$t_weight <- NULL
#
#     # do the rolling joint
#     network_dt <- temp_net
#     data.table::setkey(cutpoints_dt, trial, time)
#     #network_dt$time <- ifelse(network_dt$time > 1, network_dt$time - 1, 1)
#
#     # roll = -Inf ensures that each row gets the most recent discrete time value from D_data
#     # roll = -Inf ensures that each row gets the most recent discrete time value from D_data
#     network_dt <- cutpoints_dt[temp_net, on = c("trial", "time"), roll = -Inf]
#
#     # diagnostic table
#     network_dt <- network_dt[, {
#         sum_weight <- sum(weight, na.rm = TRUE)
#         dur <- unique(duration)
#         .(
#             sum_input = sum_weight,
#             duration = dur,
#             weight = sum_weight / dur
#         )
#     }, by = .(trial, from, to, discrete_time)]
#
#     print(network_dt[, .(trial, from, to, discrete_time, sum_input, duration)])
#
#     #remove extra column
#     network_dt$time = network_dt$discrete_time
#     network_dt$discrete_time = NULL
#     network_dt$duration = NULL
#     network_dt$sum_input = NULL
#     message("|Hi-Rez|: Pre-processing completed.")
#     return(as.data.frame(network_dt))
# }
