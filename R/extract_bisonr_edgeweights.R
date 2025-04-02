#' extract_bisonr_edgeweights
#'
#' Creates networks dataframe for STb from bisonr fit with a given number of draws from the posterior distribution of edge weights.
#'
#' @param bisonr_fit object of class "bison_model"
#' @param draws number of draws from posterior, defaults to 100
#'
#' @return dataframe with cols draw, from, to, value
#' @export
#'
#' @examples
extract_bisonr_edgeweights <- function(bisonr_fit, draws=100){

    if (!inherits(bisonr_fit, "bison_model")) stop("Please supply a bisonr fit object.")

    #yoinked code from bisonr so as not to depend
    samples <- bisonr_fit$edge_samples
    draw_ids <- sample(1:nrow(samples), draws, replace = FALSE)
    edge_samples <- samples[draw_ids, ]

    dyads <- strsplit(bisonr_fit$dyad_names, " <-> ") # split "from-to"
    dyad_matrix <- do.call(rbind, dyads)          # Convert to matrix

    # Extract 'from' and 'to' columns
    from <- dyad_matrix[, 1]
    to <- dyad_matrix[, 2]
    from_long <- rep(from, times=draws)
    to_long <- rep(to, times=draws)

    # create draw index for each row in edgelist
    draw_index <- rep(1:draws, each=bisonr_fit$num_dyads)

    # Flatten edgelist into long format, transpose to preserve order, then unlist
    values <- as.vector(t(edge_samples))

    # Create final dataframe
    long_df <- data.frame(
        draw = draw_index,
        trial = 1, # assume single trial
        from = as.integer(from_long),
        to = as.integer(to_long),
        value = values
    )

    return(long_df)
}
