
# Function to convert edge_list to networks dataframe
edge_list_to_networks <- function(edge_list) {
    networks <- edge_list
    names(networks) <- c("from", "to", "trial", "assoc")
    return(networks)
}

test_that("Network structure is consistent using distribution edgeweights.", {
    set.seed(42)  # for reproducibility

    bisonr_fit = STbayes::bisonr_fit

    #network has 10 individuals, create mock diffusion data
    diffusion_data <- data.frame(
        trial = 1,
        id = c(1:10),
        time = sample(1:101, 10, replace = FALSE),
        t_end = 100
    )


    # Import data
    data_imported <- import_user_STb(diffusion_data, bisonr_fit)

    testthat::expect_equal(data_imported$network_names, "net1")
    testthat::expect_equal(length(data_imported$logit_edge_mu), 90)
    testthat::expect_equal(dim(data_imported$logit_edge_cov), c(1,90,90))

    # Import data
    data_imported <- import_user_STb(diffusion_data, list(bisonr_fit, bisonr_fit))
    testthat::expect_equal(data_imported$network_names, c("net1", "net2"))
    testthat::expect_equal(data_imported$logit_edge_mu[1,], data_imported$logit_edge_mu[2,])

    # Import data
    data_imported <- import_user_STb(diffusion_data, STbayes::strand_results_obj)
    testthat::expect_equal(length(data_imported$logit_edge_mu), 90)
    testthat::expect_equal(dim(data_imported$logit_edge_cov), c(1,90,90))
})

