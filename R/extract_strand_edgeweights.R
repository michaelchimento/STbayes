#' extract_bisonr_edgeweights
#'
#' Creates networks dataframe for STb from STRAND results obj with a given number of draws from the posterior distribution of edge weights.
#'
#' @param strand_obj object of class "STRAND Results Object"
#' @param draws number of draws from posterior, defaults to 100
#'
#' @return dataframe with cols draw, from, to, value
#' @export
#'
#' @examples
extract_strand_edgeweights <- function(strand_obj, draws=100){

    if (class(strand_obj) != "STRAND Results Object") stop("Please supply a STRAND results object, created using STRAND::summarize_strand_results().")

    #association matrix [draws, from, to]
    ass_matrix <- strand_results_obj$samples$predicted_network_sample

    total_draws <- dim(ass_matrix)[1]
    num_nodes <- dim(ass_matrix)[2]

    # Sample the requested number of draws
    draw_ids <- sample(1:total_draws, draws, replace = FALSE)

    # Create a matrix of from-to pairs (excluding self-loops)
    from_to_pairs <- expand.grid(from = 1:num_nodes, to = 1:num_nodes)
    from_to_pairs <- from_to_pairs[from_to_pairs$from != from_to_pairs$to, ]  # Remove self-loops

    num_dyads <- nrow(from_to_pairs)

    # preallocate vectors
    draw_long <- rep(draw_ids, each = num_dyads)
    from_long <- rep(from_to_pairs$from, times = draws)
    to_long <- rep(from_to_pairs$to, times = draws)

    values <- unlist(lapply(draw_ids, function(d) as.vector(ass_matrix[d, , ])[!diag(num_nodes)]))

    # Create final dataframe
    long_df <- data.frame(
        draw = draw_long,
        from = as.integer(from_long),
        to = as.integer(to_long),
        value = values
    )

    # Create dataframe
    long_df <- data.frame(
        draw = draw_long,
        trial = 1, # assume single trial
        from = as.integer(from_long),
        to = as.integer(to_long),
        value = values
    )

    return(long_df)
}
