#' standardize_ids
#'
#' @param event_data dataframe with columns id, trial, time, t_end
#' @param networks networks dataframe with columns trial, from, to, and one or more columns of edge weights named descriptively.
#' @param ILV_c optional dataframe with columns id, and any constant individual-level variables that might be of interest
#' @param ILV_tv optional dataframe with columns trial, id, time and any time-varying variables.
#' @param t_weights Optional dataframe with columns trial, id, time and t_weight.
#'
#' @return list of dataframes with filled out ids
#' @export
#'
#' @examples
standardize_ids <- function(networks, event_data, ILV_c = NULL, ILV_tv = NULL, t_weights = NULL) {
    # Extract all unique ids from networks
    net_ids <- unique(c(as.character(networks$from), as.character(networks$to)))
    id_factor <- factor(net_ids, levels = net_ids)
    id_map <- data.frame(
        id = as.character(id_factor),
        id_numeric = as.integer(id_factor),
        stringsAsFactors = FALSE
    )

    if (!all(event_data$id %in% id_map$id)) {
        missing <- setdiff(event_data$id, id_map$id)
        stop("❌ The following IDs in `event_data` are missing from the `networks`: ", paste(missing, collapse = ", "))
    }
    # Add missing individuals to event_data as censored
    missing_ids <- setdiff(net_ids, unique(as.character(event_data$id)))
    if (length(missing_ids) > 0) {
        warning("⚠️️ The following IDs were present in networks but missing from event_data: ",
                paste(missing_ids, collapse = ", "),
                ". They will be added as censored individuals (t_end+1). You should account for all individuals in all data supplied to STbayes.")
        trials <- unique(event_data$trial)
        t_end_lookup <- aggregate(t_end ~ trial, data = event_data, FUN = max)
        censored_rows <- do.call(rbind, lapply(trials, function(trial) {
            t_end_val <- t_end_lookup$t_end[t_end_lookup$trial == trial]
            data.frame(
                trial = trial,
                id = missing_ids,
                time = t_end_val + 1,
                t_end = t_end_val
            )
        }))
        event_data <- rbind(event_data, censored_rows)
    }

    # Map IDs to numeric
    event_data <- merge(event_data, id_map, by = "id", all.x = TRUE)
    event_data <- event_data[order(event_data$trial, event_data$time), ]

    # ILV_c
    if (!is.null(ILV_c)) {
        if (!all(ILV_c$id %in% id_map$id)) {
            missing <- setdiff(ILV_c$id, id_map$id)
            stop("❌ The following IDs in `ILV_c` are missing from the `networks`: ", paste(missing, collapse = ", "))
        }
        ILV_c <- merge(ILV_c, id_map, by = "id", all.x = TRUE)
        ILV_c <- ILV_c[order(ILV_c$id_numeric), ]
    }

    if (!is.null(ILV_tv)) {
        if (!all(ILV_tv$id %in% id_map$id)) {
            missing <- setdiff(ILV_tv$id, id_map$id)
            stop("❌ The following IDs in `ILV_tv` are missing from the `networks`: ", paste(missing, collapse = ", "))
        }
        missing_ids <- setdiff(id_map$id, unique(ILV_tv$id))
        if (length(missing_ids) > 0) {
            warning("⚠️️ The following IDs were present in `networks` but missing from `ILV_tv`: ",
                    paste(missing_ids, collapse = ", "),
                    ". They will be filled with 0 for all ILV columns.")

            trials <- unique(ILV_tv$trial)
            times <- unique(ILV_tv$time)
            ILV_cols <- setdiff(names(ILV_tv), c("id", "trial", "time", "id_numeric", "trial_numeric", "discrete_time"))

            filler <- expand.grid(
                id = missing_ids,
                trial = trials,
                time = times,
                stringsAsFactors = FALSE
            )
            for (col in ILV_cols) {
                filler[[col]] <- 0
            }

            ILV_tv <- rbind(ILV_tv, filler)
        }

        ILV_tv <- merge(ILV_tv, id_map, by = "id", all.x = TRUE)
        ILV_tv <- ILV_tv[order(ILV_tv$trial, ILV_tv$time), ]
    }


    # t_weights
    if (!is.null(t_weights)) {
        if (!all(t_weights$id %in% id_map$id)) {
            missing <- setdiff(t_weights$id, id_map$id)
            stop("❌ The following IDs in `t_weights` are missing from the `networks`: ", paste(missing, collapse = ", "))
        }
        missing_ids <- setdiff(id_map$id, unique(t_weights$id))
        if (length(missing_ids) > 0) {
            warning("⚠️ The following IDs were present in `networks` but missing from `t_weights`: ",
                    paste(missing_ids, collapse = ", "),
                    ". They will be filled in with default t_weight = 0.")

            trials <- unique(t_weights$trial)
            times <- unique(t_weights$time)

            filler <- expand.grid(
                id = missing_ids,
                trial = trials,
                time = times,
                t_weight = 0
            )
            t_weights <- rbind(t_weights, filler)
        }

        t_weights <- merge(t_weights, id_map, by = "id", all.x = TRUE)
        t_weights <- t_weights[order(t_weights$trial, t_weights$time), ]
    }

    return(list(
        event_data = event_data,
        ILV_c = ILV_c,
        ILV_tv = ILV_tv,
        t_weights = t_weights,
        id_map = id_map
    ))
}
