#' create_Z_matrix: private function to create the knowledge state matrix Z[k,t,n]
#'
#' @param event_data dataframe (event_data) in import_NBDA_STb.R / import_user_STb.R
#' @export
#' @return matrix [k,t,n] where value=0 if naive at t, otherwise 1
create_Z_matrix <- function(event_data, high_res, time_type=c("interval", "real")) {
  time_type = match.arg(time_type)
  trials <- unique(event_data$trial_numeric)
  # extract unique IDs and acquisition times
  ids <- unique(event_data$id_numeric)
  times <- if (time_type=="real") max(event_data$t_end) else unique(event_data$discrete_time) # 0 = seed, max_time = censored
  Z <- array(0, dim = c(length(trials), max(times), length(ids))) # create C matrix
  for (k in trials) {
    temp_df <- event_data[event_data$trial_numeric == k, ]
    rownames(temp_df) <- NULL
    trial_max_time <- max(times)
    # fill Z based on acquisition times
    for (i in 1:nrow(temp_df)) {
      temp_id <- temp_df$id_numeric[i]
      time <- if (time_type=="real") temp_df$time[i] else temp_df$discrete_time[i]
      if (time <= trial_max_time) {
        if(high_res & time < trial_max_time){
          Z[k, (time+1):max(times), temp_id] <- 1
        } else if (high_res & time == trial_max_time){
          next
        }
        else if (!high_res & time < trial_max_time) {
          Z[k, (time + 1):max(times), temp_id] <- 1 # add time +1 to account for demos/censored.
        }
        else if (!high_res & time == trial_max_time) {
          next
        }
      }
    }
  }
  return(Z)
}
