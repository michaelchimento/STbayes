
# Function to convert edge_list to networks dataframe
edge_list_to_networks <- function(edge_list) {
    networks <- edge_list
    names(networks) <- c("from", "to", "trial", "assoc")
    return(networks)
}

test_that("Network structure is consistent using distribution edgeweights.", {
    set.seed(42)  # for reproducibility

    bisonr_fit = STbayes::bisonr_fit
    networks = extract_bisonr_edgeweights(bisonr_fit, draws=100)
    networks$value = scales::rescale(networks$value) #networks can now be used in import_user_STb following normal workflow

    #network has 10 individuals, create mock diffusion data
    diffusion_data <- data.frame(
        trial = 1,
        id = c(1:10),
        time = sample(1:101, 10, replace = FALSE),
        t_end = 100
    )


    # Import data
    data_imported <- import_user_STb(diffusion_data, networks)

    # Extract A_assoc matrix
    A_assoc <- data_imported$A_value

    # Verify that the values in A_assoc match the original edge_list associations
    for (i in seq_len(nrow(networks))) {
        from <- networks$from[i]
        to <- networks$to[i]
        trial <- networks$trial[i]
        draw <- networks$draw[i]
        assoc_value <- networks$value[i]
        testthat::expect_equal(A_assoc[trial, 1, draw, from, to], assoc_value)
    }
    # Ensure the network is properly repeated for all timesteps if static
    for (t in seq_len(data_imported$T_max)) {
        testthat::expect_equal(A_assoc[trial, t, draw, from, to], assoc_value)
    }
})

