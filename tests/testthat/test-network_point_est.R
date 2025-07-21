test_that("Network structure is consistent using point estimate edgeweights.", {
    # Load diffusion data
    event_data <- STbayes::event_data

    # Create mock edge_list
    networks <- STbayes::edge_list
    g <- igraph::graph_from_edgelist(as.matrix(networks[1:2]), directed = FALSE)
    adj_matrix <- igraph::as_adjacency_matrix(g, attr = NULL, sparse = FALSE)


    # Import data
    data_imported <- STbayes::import_user_STb(event_data, networks)

    # Extract A_assoc matrix
    A_assoc <- data_imported$A

    # Verify that the values in A_assoc match the original edge_list associations
    for (i in seq_len(nrow(networks))) {
        focal <- networks$focal[i]
        other <- networks$other[i]
        trial <- networks$trial[i]
        assoc_value <- networks$assoc[i]
        testthat::expect_equal(A_assoc[1, trial, 1, focal, other], adj_matrix[focal, other])
    }
    # Ensure the network is properly repeated for all timesteps if static
    for (t in seq_len(data_imported$T_max)) {
        testthat::expect_equal(A_assoc[1, trial, t, focal, other], assoc_value)
    }
})
