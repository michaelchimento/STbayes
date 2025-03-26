
# Function to convert edge_list to networks dataframe
edge_list_to_networks <- function(edge_list) {
    networks <- edge_list
    names(networks) <- c("from", "to", "trial", "assoc")
    return(networks)
}

test_that("Network structure is consistent using point estimate edgeweights.", {
    set.seed(42)  # for reproducibility

    # Load diffusion data
    diffusion_data <- STbayes::event_data

    # Create mock edge_list
    edge_list <- STbayes::edge_list
    g = igraph::graph_from_edgelist(as.matrix(STbayes::edge_list[1:2]), directed = FALSE)
    adj_matrix = igraph::as_adjacency_matrix(g,attr=NULL, sparse = FALSE)

    # Convert edge_list to networks format
    networks <- edge_list_to_networks(edge_list)

    # Import data
    data_imported <- import_user_STb(diffusion_data, networks)

    # Extract A_assoc matrix
    A_assoc <- data_imported$A_assoc

    # Verify that the values in A_assoc match the original edge_list associations
    for (i in seq_len(nrow(edge_list))) {
        from <- edge_list$from[i]
        to <- edge_list$to[i]
        trial <- edge_list$trial[i]
        assoc_value <- edge_list$assoc[i]
        testthat::expect_equal(A_assoc[trial, 1, from, to], adj_matrix[from,to])
    }
    # Ensure the network is properly repeated for all timesteps if static
    for (t in seq_len(data_imported$T_max)) {
        testthat::expect_equal(A_assoc[trial, t, from, to], assoc_value)
    }
})

