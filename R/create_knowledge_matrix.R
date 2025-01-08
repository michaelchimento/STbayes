#' create_knowledge_matrix: private function to create the knowledge state matrix C[t,n]
#'
#' @param df dataframe (diffusion_data) in import_NBDA_data.R / create_STb_data.R
#'
#' @return matrix [t,n] where value=0 if naive at t, otherwise 1
create_knowledge_matrix <- function(diffusion_data) {
    trials = unique(diffusion_data$trial_numeric)
    # extract unique IDs and acquisition times
    ids <- unique(diffusion_data$id_numeric)
    times <- unique(diffusion_data$discrete_time) #0 = seed, max_time = censored
    C <- array(0, dim = c(length(trials), max(times), length(ids))) #create C matrix
    for (k in trials) {
        temp_df = diffusion_data[diffusion_data$trial_numeric==k,]
        rownames(temp_df) = NULL
        trial_max_time = max(temp_df$discrete_time)
        # fill C based on acquisition times
        for (i in 1:nrow(temp_df)) {
            temp_id <- temp_df$id_numeric[i]
            time <- temp_df$discrete_time[i]
            if (time != trial_max_time) {
                C[k,(time+1):max(times), temp_id] <- 1 #add time +1 to account for demos/censored.
            }
        }
    }
    return(C)
}

