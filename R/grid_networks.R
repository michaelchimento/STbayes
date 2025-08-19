utils::globalVariables(c(
    "trial", "from", "other", ".", "t_end", "trial_numeric", "time", "t_max", "dummy"
))
#' grid_networks()
#'
#' Helper function that completes network df with all dyad & trial combinations.
#'
#' @param networks network df
#' @param event_data event df
#'
#' @importFrom data.table as.data.table setkey merge.data.table
#' @return network dataframe
grid_networks <- function(event_data, networks) {
    networks <- data.table::as.data.table(networks)
    event_data <- data.table::as.data.table(event_data)

    # recode trial, from, other as integer factors
    networks[, `:=`(
        trial = as.integer(factor(trial)),
        focal = as.integer(factor(focal)),
        other = as.integer(factor(other))
    )]

    expected_rows <- sum(sapply(split(event_data, event_data$trial_numeric), function(df) {
        P <- length(unique(df$id_numeric))
        T <- max(df$t_end)
        P * (P - 1) * T
    }))

    # merge with full network grid if needed
    if (nrow(networks) < expected_rows) {
        message("|Hi-Rez|: Completing network grid.")
        # get full combinations
        ids <- sort(unique(event_data$id_numeric))
        times <- event_data[, .(time = seq_len(max(t_end))), by = trial_numeric]
        trials <- unique(event_data$trial_numeric)
        trial_id_dt <- data.table::data.table(
            trial_numeric = trials,
            trial = seq_along(trials) # ensure consistent trial ids
        )

        focal_other <- data.table::CJ(focal = ids, other = ids)[focal != other]
        full_grid <- data.table::CJ(trial_numeric = trials, time = 1, sorted = FALSE)[, time := NULL]
        full_grid <- event_data[, .(t_max = max(t_end)), by = trial_numeric][
            , .(time = seq_len(t_max)),
            by = trial_numeric
        ]

        # add dummy column to both
        full_grid[, dummy := 1]
        focal_other[, dummy := 1]

        # Cartesian join on dummy
        full_net_grid <- data.table::merge.data.table(full_grid, focal_other, by = "dummy", allow.cartesian = TRUE)[, dummy := NULL]
        full_net_grid <- data.table::merge.data.table(full_net_grid, trial_id_dt, by = "trial_numeric")
        full_net_grid[, trial_numeric := NULL] # match final format

        # merge with existing network data
        networks <- data.table::merge.data.table(
            full_net_grid,
            networks,
            by = c("trial", "time", "focal", "other"),
            all.x = TRUE
        )
        networks[is.na(networks)] <- 0
    }

    return(as.data.frame(networks))
}
