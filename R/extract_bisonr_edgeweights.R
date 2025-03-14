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

    if (class(bisonr_fit) != "bison_model") stop("Please supply a bisonr fit object.")

    #yoinked code from bisonr so as not to depend
    samples <- matrix(as.numeric(bisonr_fit$edge_samples), ncol=ncol(bisonr_fit$edge_samples))
    draw_ids <- sample(1:dim(samples)[1], draws)
    edge_samples <- samples[draw_ids, ]

    dyads <- strsplit(bisonr_fit$dyad_names, " <-> ") # split "from-to"
    dyad_matrix <- do.call(rbind, dyads)          # Convert to matrix

    # Extract 'from' and 'to' columns
    from <- dyad_matrix[, 1]
    to <- dyad_matrix[, 2]

    num_dyads <- ncol(edge_samples)

    # Create a 'draw' index for each row in edgelist
    draw_index <- rep(1:draws, each=bisonr_fit$num_dyads)

    # Flatten edgelist into long format
    values <- as.vector(t(edge_samples))  # Transpose to preserve order, then unlist

    # Replicate 'from' and 'to' for each draw
    from_long <- rep(from, times=draws)
    to_long <- rep(to, times=draws)

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
