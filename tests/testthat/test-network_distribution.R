testthat::test_that("Importing single 3, 4, and 5 dimensional network arrays works.", {
    # network has 10 individuals, create mock diffusion data
    event_data <- data.frame(
        trial = 1,
        # bison stores ids as characters, ids should match across dataframes
        id = as.character(c(1:10)),
        time = sample(1:101, 10, replace = FALSE),
        t_end = 100
    )

    # mock 1000 posterior draws.
    network_3D <- array(rnorm(1000 * 10 * 10, mean = 0, sd = 1),
        dim = c(1000, 10, 10),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    network_4D <- array(
        rnorm(1 * 1000 * 10 * 10),
        dim = c(1, 1000, 10, 10),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    network_5D <- array(
        rnorm(1 * 10 * 1000 * 10 * 10),
        dim = c(1, 10, 1000, 10, 10),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    data_list <- import_user_STb(event_data,
        networks = network_3D
    )
    model <- generate_STb_model(data_list)
    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(model)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )

    data_list <- import_user_STb(event_data,
        networks = network_4D
    )
    model <- generate_STb_model(data_list)
    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(model)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )

    data_list <- import_user_STb(event_data,
        networks = network_5D
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

testthat::test_that("Importing 3D multi-network distribution edge weights works.", {
    # network has 10 individuals, create mock diffusion data
    event_data <- data.frame(
        trial = 1,
        # bison stores ids as characters, ids should match across dataframes
        id = as.character(c(1:10)),
        time = sample(1:101, 10, replace = FALSE),
        t_end = 100
    )

    # mock 1000 posterior draws.
    network_3D1 <- array(rnorm(1000 * 10 * 10, mean = 0, sd = 1),
        dim = c(1000, 10, 10),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    network_3D2 <- array(rnorm(1000 * 10 * 10, mean = 0, sd = 1),
        dim = c(1000, 10, 10),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    data_list <- import_user_STb(event_data,
        networks = list(network_3D1, network_3D2)
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

testthat::test_that("Importing 4D multi-network distribution edge weights works.", {
    # network has 10 individuals, create mock diffusion data
    event_data <- data.frame(
        trial = 1,
        # bison stores ids as characters, ids should match across dataframes
        id = as.character(c(1:10)),
        time = sample(1:101, 10, replace = FALSE),
        t_end = 100
    )

    # mock 1000 posterior draws.
    network_4D1 <- array(
        rnorm(1 * 1000 * 10 * 10),
        dim = c(1, 1000, 10, 10),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    network_4D2 <- array(
        rnorm(1 * 1000 * 10 * 10),
        dim = c(1, 1000, 10, 10),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    data_list <- import_user_STb(event_data,
        networks = list(network_4D1, network_4D2)
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

testthat::test_that("Importing 4D multi-network distribution edge weights works.", {
    # network has 10 individuals, create mock diffusion data
    event_data <- data.frame(
        trial = 1,
        # bison stores ids as characters, ids should match across dataframes
        id = as.character(c(1:10)),
        time = sample(1:101, 10, replace = FALSE),
        t_end = 100
    )

    # mock 1000 posterior draws.
    network_5D1 <- array(
        rnorm(1 * 10 * 1000 * 10 * 10),
        dim = c(1, 10, 1000, 10, 10),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    network_5D2 <- array(
        rnorm(1 * 10 * 1000 * 10 * 10),
        dim = c(1, 10, 1000, 10, 10),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    data_list <- import_user_STb(event_data,
        networks = list(network_5D1, network_5D2)
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

testthat::test_that("Network structure is consistent using bisonr fit.", {
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
