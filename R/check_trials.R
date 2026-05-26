#' check_trials()
#'
#' Helper function that ensures all trials are present in networks and events df
#'
#' @param event_data dataframe with columns id, trial, time, t_end
#' @param networks networks dataframe with columns trial, focal, other, and one or more columns of edge weights named descriptively.
#' @param ILV_tv dataframe of time-varying ILVs
#' @param t_weights dataframe of transmission weights, can be timevarying
#'
#' @return error message if trials are missing from dataframes
#' @noRd
check_trials <- function(event_data, networks, ILV_tv = NULL, t_weights = NULL) {
    event_trials <- unique(event_data$trial)

    if (inherits(networks, "data.frame")) {
        network_trials <- unique(networks$trial)

        missing_in_networks <- setdiff(event_trials, network_trials)
        missing_in_events <- setdiff(network_trials, event_trials)

        if (length(missing_in_networks) > 0) {
            stop("\u274C Trials present in `event_data` but missing in `networks`: ", paste(missing_in_networks, collapse = ", "))
        }

        if (length(missing_in_events) > 0) {
            stop("\u274C Trials present in `networks` but missing in `event_data`: ", paste(missing_in_events, collapse = ", "))
        }
        return(NULL)
    }

    if (!is.list(networks)) {
        networks <- list(networks)
    }

    is_external_fit <- vapply(
        networks,
        function(x) inherits(x, "bison_model") || inherits(x, "STRAND Results Object"),
        TRUE
    )

    is_array <- vapply(
        networks,
        function(x) inherits(x, "array"),
        TRUE
    )

    if (any(is_external_fit) && length(event_trials) > 1) {
        warning("\u26A0\ufe0f STbayes does not currently support trial-specific networks when networks are defined as a bisonR or STRAND model fit. The same posterior network will be used for all trials.")
    }

    if (any(is_array)) {
        for (i in which(is_array)) {
            net_obj <- networks[[i]]
            nd <- length(dim(net_obj))

            if (!nd %in% c(3, 4, 5)) {
                stop(
                    "\u274C User-supplied posterior network arrays must be 3D [draw, focal_ID, other_ID], ",
                    "4D [trial, draw, focal_ID, other_ID], or ",
                    "5D [trial, time, draw, focal_ID, other_ID]."
                )
            }

            if (nd %in% c(4, 5)) {
                dim_names <- names(dimnames(net_obj))
                trial_pos <- match("trial", dim_names)

                if (is.na(trial_pos)) {
                    stop(
                        "\u274C Posterior network array ", i,
                        " includes a trial dimension, but that dimension is not named `trial`."
                    )
                }

                trial_dim <- dim(net_obj)[trial_pos]

                if (trial_dim != length(event_trials)) {
                    stop(
                        "\u274C Posterior network array ", i,
                        " has ", trial_dim, " trials, but `event_data` has ",
                        length(event_trials), " trials."
                    )
                }
            }
        }
    }

    if (!is.null(ILV_tv)) {
        ILV_tv_trials <- unique(ILV_tv$trial)
        missing_in_ILV_tv <- setdiff(event_trials, ILV_tv_trials)
        if (length(missing_in_ILV_tv) == 0) {
            NULL
        } else {
            stop("\u274C Trials present in `event_data` but missing in `ILV_tv`: ", paste(missing_in_ILV_tv, collapse = ", "))
        }
    }

    if (!is.null(t_weights)) {
        t_weights_trials <- unique(t_weights$trial)
        missing_in_t_weights <- setdiff(event_trials, t_weights_trials)
        if (length(missing_in_t_weights) == 0) {
            NULL
        } else {
            stop("\u274C Trials present in `event_data` but missing in `t_weights`: ", paste(missing_in_t_weights, collapse = ", "))
        }
    }
}
