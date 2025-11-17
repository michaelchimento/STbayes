testthat::test_that("Importing single and multi-network distribution edge weights works.", {
    # network has 10 individuals, create mock diffusion data
    event_data <- data.frame(
        trial = 1,
        # bison stores ids as characters, ids should match across dataframes
        id = as.character(c(1:10)),
        time = sample(1:101, 10, replace = FALSE),
        t_end = 100
    )

    # mock 1000 posterior draws.
    association <- array(rnorm(1000 * 10 * 10, mean = 0, sd = 1),
        dim = c(1000, 10, 10)
    )

    data_list <- import_user_STb(event_data, networks = association)
    model <- generate_STb_model(data_list)

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(model)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )

    # it's also possible to do multi-network models
    visual_contact <- array(rnorm(1000 * 10 * 10, mean = 0, sd = 1),
        dim = c(1000, 10, 10)
    )

    data_list <- import_user_STb(event_data,
        networks = list(association, visual_contact)
    )
    model <- generate_STb_model(data_list)
    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(model)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

testthat::test_that("Network structure is consistent using distribution edgeweights.", {
    networks <- STbayes::bisonr_fit

    # network has 10 individuals, create mock diffusion data
    events <- data.frame(
        trial = 1,
        id = c(1:10),
        time = sample(1:101, 10, replace = FALSE),
        t_end = 100
    )

    events$id <- as.character(events$id)

    # Import data
    data_imported <- import_user_STb(events, networks, network_type = "directed")

    testthat::expect_equal(data_imported$network_names, "net1")
    testthat::expect_equal(length(data_imported$logit_edge_mu), 90)
    testthat::expect_equal(dim(data_imported$logit_edge_cov), c(1, 90, 90))

    # Import data
    data_imported <- import_user_STb(events, list(networks, networks))
    testthat::expect_equal(data_imported$network_names, c("net1", "net2"))
    testthat::expect_equal(data_imported$logit_edge_mu[1, ], data_imported$logit_edge_mu[2, ])

    # Import data
    events$id <- as.numeric(events$id)
    data_imported <- import_user_STb(event_data = events, networks = STbayes::strand_results_obj)
    testthat::expect_equal(length(data_imported$logit_edge_mu), 90)
    testthat::expect_equal(dim(data_imported$logit_edge_cov), c(1, 90, 90))
})
