#' create_knowledge_matrix: private function to create the knowledge state matrix C[t,n]
#'
#' @param df dataframe (diffusion_data) in import_NBDA_data.R / create_STb_data.R
#'
#' @return matrix [t,n] where value=0 if naive at t, otherwise 1
create_knowledge_matrix <- function(df) {
    trials = unique(df$trial_numeric)
    # extract unique IDs and acquisition times
    ids <- unique(df$id_numeric)
    times <- unique(df$discrete_time) #0 = seed, max_time = censored
    C <- array(0, dim = c(length(trials), max(times), length(ids))) #create C matrix
    for (k in trials) {
        # fill C based on acquisition times
        for (i in ids) {
            temp_id <- df$id_numeric[i]
            time <- df$discrete_time[i]
            if (time != max(times)) {
                C[k,(time+1):max(times), temp_id] <- 1 #add time +1 to account for demos/censored.
            }
        }
    }
    return(C)
}

