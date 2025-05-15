#' Helper function for preprocessing "high resolution" data
#'
#' @param t_weights
#' @param networks
#' @param event_data
#'
#' @importFrom data.table as.data.table setkey
#' @return
#' @export
#'
#' @examples


OLD_process_networks_x_weights_hires <- function(event_data, t_weights, networks, D_data) {
    t_weights$trial = as.integer(as.factor(t_weights$trial))
    t_weights$id = as.integer(as.factor(t_weights$id))
    networks$trial = as.integer(as.factor(networks$trial))
    networks$from = as.integer(as.factor(networks$from))
    networks$to = as.integer(as.factor(networks$to))

    expected_rows <- sum(sapply(split(event_data, event_data$trial_numeric), function(df) {
        P <- length(unique(df$id_numeric))
        T <- max(df$t_end)
        P * (P - 1) * T  # directed dyads × timesteps
    }))
    if (nrow(networks)<expected_rows){
        message("|Hi-Rez|: Completing network grid.")
        full_net_grid <- complete_network_grid(event_data)
        networks <- merge(full_net_grid, networks, by = c("trial", "time", "from", "to"), all.x = TRUE)
        networks[is.na(networks)] <- 0
    }

    names(t_weights)[names(t_weights) == "id"] <- "to"
    names(D_data)[names(D_data) == "trial_numeric"] <- "trial"
    # join with networks

    message("|Hi-Rez|: Pre-processing network x cue data.")
    temp_net <- merge(networks, t_weights, by = c("trial", "time", "to"), all.x = TRUE)
    summary(temp_net)
    # # weight columns
    weight_cols <- setdiff(names(temp_net), c("trial", "from", "to", "time", "t_weight"))
    for (col in weight_cols) {
        temp_net[[col]] <- temp_net[[col]] * temp_net$t_weight
    }
    temp_net$t_weight <- NULL
    #
    # # assign numeric IDs
    # event_data$id_numeric <- as.numeric(as.factor(event_data$id))
    # event_data$trial_numeric <- as.numeric(as.factor(event_data$trial))
    # event_data <- event_data[order(event_data$trial_numeric, event_data$time), ]
    #
    # # create discrete time
    # event_data$discrete_time <- NA
    # event_data$discrete_time[event_data$time != 0] <- with(
    #     event_data[event_data$time != 0, ],
    #     ave(time, trial_numeric, FUN = function(x) as.numeric(as.factor(x)))
    # )
    # event_data$discrete_time[event_data$time == 0] <- 0
    #
    # # clip times above t_end
    # temp_data <- event_data
    # temp_data$time[temp_data$time > temp_data$t_end] <- temp_data$t_end[temp_data$time > temp_data$t_end]
    # D_data <- unique(temp_data[order(temp_data$trial, temp_data$time), c("trial", "time", "discrete_time")])
    # D_data <- D_data[D_data$time != 0, ]

    # do the rolling joint
    network_dt <- data.table::as.data.table(temp_net)
    cutpoints_dt <- data.table::as.data.table(D_data)
    data.table::setkey(cutpoints_dt, trial, time)
    #network_dt$time <- ifelse(network_dt$time > 1, network_dt$time - 1, 1)

    # roll = -Inf ensures that each row gets the most recent discrete time value from D_data
    network_dt <- cutpoints_dt[network_dt, on = c("trial", "time"), roll = -Inf]

    # aggregate, .SDcols = specifies which columns
    network_dt <- network_dt[, lapply(.SD, sum), by = .(trial, from, to, discrete_time),
                             .SDcols = weight_cols]
    network_dt <- merge(network_dt, D_data, by.x = c("trial", "discrete_time"), by.y = c("trial", "discrete_time"))
    for (col in weight_cols) {
        network_dt[[col]] <- network_dt[[col]] / network_dt$duration
    }

    #remove extra column
    network_dt$time = network_dt$discrete_time
    network_dt$discrete_time = NULL
    network_dt$duration = NULL
    message("|Hi-Rez|: Pre-processing completed.")
    return(as.data.frame(network_dt))
}
