#' Title
#'
#' @param event_data
#' @param networks
#' @param t_weights
#' @param D_data
#' @param transmission_type
#'
#'#' @importFrom data.table as.data.table setkey fifelse
#' @return
#' @export
#'
#' @examples
#'
#'
process_hires_complex_k <- function(event_data, networks, t_weights, D_data) {
    message("|Hi-Rez|: Pre-processing for freqdep_k.")

    networks <- as.data.table(networks)
    t_weights <- as.data.table(t_weights)
    D_data <- as.data.table(D_data)
    event_data <- as.data.table(event_data)

    t_weights[, trial := as.integer(as.factor(trial))]
    t_weights[, to := as.integer(as.factor(id))][, id := NULL]
    D_data[, trial := as.integer(as.factor(trial_numeric))][, trial_numeric := NULL]
    event_data[, `:=`(
        trial = as.integer(as.factor(trial)),
        id = as.integer(as.factor(id))
    )]

    # Create binary Z column: whether "to" individual had acquired behavior by each time
    acq_time <- event_data[, .(acquisition_time = time), by = .(trial, id)]
    acq_time[, to := as.integer(as.factor(id))][, id := NULL]
    networks <- data.table::merge.data.table(networks, acq_time, by = c("trial", "to"), all.x = TRUE)
    networks[, z := as.integer(time > acquisition_time)]
    networks[, zn := as.integer(1-z)]

    # Apply transmission weights
    temp_net <- merge(networks, t_weights, by = c("trial", "time", "to"), all.x = TRUE)
    weight_cols <- setdiff(names(temp_net), c("trial", "from", "to", "time", "t_weight", "z", "zn", "acquisition_time"))

    N_networks <- length(weight_cols)
    K <- max(temp_net$trial)
    T_discrete <- max(D_data$discrete_time)
    P <- max(c(temp_net$from, temp_net$to))

    prop_array <- array(0, dim = c(N_networks, K, T_discrete, P))
    setkey(D_data, trial, time)

    for (n in seq_along(weight_cols)) {
        col <- weight_cols[n]
        # Don't apply t_weight to denom
        prop_rows <- temp_net[, .(
            numer = sum(get(col) * t_weight * z),
            denom = sum(get(col)* t_weight * z) + sum(get(col) * zn)
        ), by = .(trial, time, from)][, .(
            trial, time, id = from,
            prop = data.table::fifelse(denom > 0, numer / denom, 0)
        )]

        prop_dt <- D_data[prop_rows, on = .(trial, time), roll = -Inf]
        prop_avg <- prop_dt[, .(prop = sum(prop) / duration), by = .(trial, discrete_time, id)]
        for (i in seq_len(nrow(prop_avg))) {
            k <- prop_avg$trial[i]
            t <- prop_avg$discrete_time[i]
            p <- prop_avg$id[i]
            prop_array[n, k, t, p] <- prop_avg$prop[i]
        }
    }

    return(prop_array)
}
# process_hires_complex_k <- function(event_data, networks, t_weights, Z, D_data) {
#     message("|Hi-Rez|: Pre-processing for freqdep_k.")
#
#     # assume Z, Zn are [K, T_max, P]
#     K <- dim(Z)[1]
#     T_max <- dim(Z)[2]
#     P <- dim(Z)[3]
#
#     # convert to data.table
#     networks <- as.data.table(networks)
#     t_weights <- as.data.table(t_weights)
#     D_data <- as.data.table(D_data)
#
#     # ensure integer trial and ID mapping
#     t_weights[, trial := as.integer(as.factor(trial))]
#     t_weights[, to := as.integer(as.factor(id))][, id := NULL]
#     D_data[, trial := as.integer(as.factor(trial_numeric))][, trial_numeric := NULL]
#
#     # merge t_weights into networks
#     temp_net <- data.table::merge.data.table(networks, t_weights, by = c("trial", "time", "to"), all.x = TRUE)
#
#     weight_cols <- setdiff(names(temp_net), c("trial", "from", "to", "time", "t_weight"))
#     temp_net[, (weight_cols) := lapply(.SD, function(x) x * t_weight), .SDcols = weight_cols]
#     N_networks <- length(weight_cols)
#     T_discrete <- max(D_data$discrete_time)
#
#     prop_array <- array(0, dim = c(N_networks, K, T_discrete, P))
#
#     # precompute discrete time mapping
#     setkey(D_data, trial, time)
#
#     for (n in seq_along(weight_cols)) {
#         col <- weight_cols[n]
#
#         prop_rows <- do.call(rbind, lapply(split(temp_net, list(temp_net$trial, temp_net$time)), function(df) {
#         trial <- unique(df$trial)
#         time <- unique(df$time)
#         out <- lapply(unique(df$from), function(id) {
#             sub <- df[df$from == id, ]
#             z_vals <- Z[trial, time, sub$to]
#             w <- sub[[col]]
#             numer <- sum(w * z_vals)
#             denom <- numer + sum(w * (1 - z_vals))
#             prop <- if (denom > 0) numer / denom else 0
#             data.frame(trial = trial, time = time, id = id, prop = prop)
#         })
#         do.call(rbind, out)
#     }))
#
#         prop_dt <- data.table::as.data.table(prop_rows)
#             cutpoints_dt <- data.table::as.data.table(D_data)
#             data.table::setkey(cutpoints_dt, trial, time)
#             data.table::setkey(prop_dt, trial, time)
#             prop_dt <- cutpoints_dt[prop_dt, roll = -Inf]
#             # Aggregate and normalize
#             prop_avg <- prop_dt[, .(prop = sum(prop) / duration), by = .(trial, discrete_time, id)]
#             for (i in seq_len(nrow(prop_avg))) {
#                 k <- prop_avg$trial[i]
#                 t <- prop_avg$discrete_time[i]
#                 p <- prop_avg$id[i]
#                 prop_array[n, k, t, p] <- prop_avg$prop[i]
#             }
#     }
#
#     return(prop_array)
# }

# process_hires_complex_k <- function(event_data, networks, t_weights, Z, D_data) {
#     message("|Hi-Rez|: Pre-processing for freqdep_k v2.")
#
#     t_weights$trial = as.integer(as.factor(t_weights$trial))
#     t_weights$id = as.integer(as.factor(t_weights$id))
#
#     # assumes: Z and Zn are [K, T_max, P]
#     K <- dim(Z)[1]
#     T_max <- dim(Z)[2]
#     P <- dim(Z)[3]
#
#     names(t_weights)[names(t_weights) == "id"] <- "to"
#     names(D_data)[names(D_data) == "trial_numeric"] <- "trial"
#     # Multiply edge weights by t_weight
#     temp_net <- data.table::merge.data.table(networks, t_weights, by = c("trial", "time", "to"), all.x = TRUE)
#     weight_cols <- setdiff(names(temp_net), c("trial", "from", "to", "time", "t_weight"))
#     N_networks <- length(weight_cols)
#     T_discrete <- max(D_data$discrete_time)
#
#     prop_array <- array(0, dim = c(N_networks, K, T_discrete, P))
#
#     for (n in seq_along(weight_cols)) {
#         col <- weight_cols[n]
#         temp_net[[col]] <- temp_net[[col]] * temp_net$t_weight
#
#         prop_rows <- do.call(rbind, lapply(split(temp_net, list(temp_net$trial, temp_net$time)), function(df) {
#             trial <- unique(df$trial)
#             time <- unique(df$time)
#             out <- lapply(unique(df$from), function(id) {
#                 sub <- df[df$from == id, ]
#                 z_vals <- Z[trial, time, sub$to]
#                 w <- sub[[col]]
#                 numer <- sum(w * z_vals)
#                 denom <- numer + sum(w * (1 - z_vals))
#                 prop <- if (denom > 0) numer / denom else 0
#                 data.frame(trial = trial, time = time, id = id, prop = prop)
#             })
#             do.call(rbind, out)
#         }))
#
#
#     # Rolling join to discrete_time
#     prop_dt <- data.table::as.data.table(prop_rows)
#     cutpoints_dt <- data.table::as.data.table(D_data)
#     data.table::setkey(cutpoints_dt, trial, time)
#     data.table::setkey(prop_dt, trial, time)
#     prop_dt <- cutpoints_dt[prop_dt, roll = -Inf]
#     # Aggregate and normalize
#     prop_avg <- prop_dt[, .(prop = sum(prop) / duration), by = .(trial, discrete_time, id)]
#     for (i in seq_len(nrow(prop_avg))) {
#         k <- prop_avg$trial[i]
#         t <- prop_avg$discrete_time[i]
#         p <- prop_avg$id[i]
#         prop_array[n, k, t, p] <- prop_avg$prop[i]
#     }
#     }
#     return(prop_array)
# }
