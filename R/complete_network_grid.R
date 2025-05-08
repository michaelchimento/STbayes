#' pads out networks with all combinations of IDs and trials
#'
#' @param event_data
#'
#' @return a dataframe containing all trials, dyads and times (but not inter-event intervals)
#' @export
#'
#' @examples
complete_network_grid <- function(event_data) {
    trials <- sort(unique(event_data$trial_numeric))

    full_list <- lapply(trials, function(tr) {
        trial_rows <- event_data[event_data$trial_numeric == tr, ]

        # IDs and t_end specific to this trial
        ids <- sort(unique(trial_rows$id_numeric))
        t_max <- max(trial_rows$t_end)
        times <- seq_len(t_max)

        # build per-trial grid
        from_to <- expand.grid(from = ids, to = ids)
        from_to <- from_to[from_to$from != from_to$to, ]  # drop self-loops
        grid <- expand.grid(trial = tr, time = times)
        full_grid <- merge(grid, from_to, by = NULL)

        return(full_grid)
    })

    do.call(rbind, full_list)
}
