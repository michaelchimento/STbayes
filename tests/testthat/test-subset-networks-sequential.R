test_that("subset_networks_sequential subsets data.frame networks up to current trial and time", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 2, 2, 2),
        discrete_time = c(1, 2, 1, 2, 3)
    )

    network_data <- data.frame(
        trial_numeric = c(1, 1, 1, 2, 2, 2, 3),
        discrete_time = c(1, 2, 3, 1, 2, 3, 1),
        focal = LETTERS[1:7],
        other = LETTERS[2:8],
        weight = seq_len(7)
    )

    out <- subset_networks_sequential(
        i = 4,
        network_data = network_data,
        event_data = event_data
    )

    expected <- network_data[
        network_data$trial_numeric <= 2 &
            !(network_data$trial_numeric == 2 & network_data$discrete_time > 3),
    ]

    expect_s3_class(out, "data.frame")
    expect_equal(out, expected)
    expect_true(all(out$trial_numeric <= 2))
    expect_false(any(out$trial_numeric == 2 & out$discrete_time > 3))
})


test_that("subset_networks_sequential keeps all earlier trial network rows", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 2, 2),
        discrete_time = c(1, 2, 1, 1)
    )

    network_data <- data.frame(
        trial_numeric = c(1, 1, 1, 2, 2),
        discrete_time = c(1, 2, 3, 1, 2),
        focal = LETTERS[1:5],
        other = LETTERS[2:6],
        weight = seq_len(5)
    )

    out <- subset_networks_sequential(
        i = 3,
        network_data = network_data,
        event_data = event_data
    )

    expect_true(all(network_data$trial_numeric[network_data$trial_numeric == 1] %in% out$trial_numeric))
    expect_equal(out[out$trial_numeric == 1, ], network_data[network_data$trial_numeric == 1, ])
    expect_false(any(out$trial_numeric == 2 & out$discrete_time > 2))
})


test_that("subset_networks_sequential leaves 3D posterior arrays unchanged", {
    draws <- 4
    P <- 3

    net <- array(
        seq_len(draws * P * P),
        dim = c(draws, P, P),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    event_data <- data.frame(
        trial_numeric = c(1, 1, 2, 2),
        discrete_time = c(1, 2, 1, 2)
    )

    out <- subset_networks_sequential(
        i = 3,
        network_data = net,
        event_data = event_data
    )

    expect_type(out, "list")
    expect_length(out, 1)
    expect_equal(out[[1]], net)
    expect_equal(dim(out[[1]]), dim(net))
    expect_equal(names(dimnames(out[[1]])), c("draw", "focal_ID", "other_ID"))
})


test_that("subset_networks_sequential subsets 4D trial-varying arrays up to current trial", {
    K <- 4
    draws <- 3
    P <- 3

    net <- array(
        seq_len(K * draws * P * P),
        dim = c(K, draws, P, P),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    event_data <- data.frame(
        trial_numeric = c(1, 1, 2, 2, 3),
        discrete_time = c(1, 2, 1, 2, 1)
    )

    out <- subset_networks_sequential(
        i = 4,
        network_data = net,
        event_data = event_data
    )

    expected <- net[1:2, , , , drop = FALSE]

    expect_type(out, "list")
    expect_length(out, 1)
    expect_equal(out[[1]], expected)
    expect_equal(dim(out[[1]]), c(2, draws, P, P))
    expect_equal(names(dimnames(out[[1]])), c("trial", "draw", "focal_ID", "other_ID"))
})


test_that("subset_networks_sequential subsets 5D arrays and zero-fills future times in current trial", {
    K <- 4
    T_max <- 5
    draws <- 3
    P <- 3

    net <- array(
        seq_len(K * T_max * draws * P * P),
        dim = c(K, T_max, draws, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    event_data <- data.frame(
        trial_numeric = c(1, 1, 2, 2, 3),
        discrete_time = c(1, 2, 1, 2, 1)
    )

    out <- subset_networks_sequential(
        i = 4,
        network_data = net,
        event_data = event_data
    )

    expected <- net[1:2, , , , , drop = FALSE]
    expected[2, 4:5, , , ] <- 0

    expect_type(out, "list")
    expect_length(out, 1)
    expect_equal(out[[1]], expected)
    expect_equal(dim(out[[1]]), c(2, T_max, draws, P, P))
    expect_equal(names(dimnames(out[[1]])), c("trial", "time", "draw", "focal_ID", "other_ID"))
})


test_that("subset_networks_sequential does not zero-fill 5D arrays when current time is final time", {
    K <- 3
    T_max <- 4
    draws <- 3
    P <- 3

    net <- array(
        seq_len(K * T_max * draws * P * P),
        dim = c(K, T_max, draws, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    event_data <- data.frame(
        trial_numeric = c(1, 1, 2, 2),
        discrete_time = c(1, 2, 3, 4)
    )

    out <- subset_networks_sequential(
        i = 4,
        network_data = net,
        event_data = event_data
    )

    expected <- net[1:2, , , , , drop = FALSE]

    expect_equal(out[[1]], expected)
    expect_equal(dim(out[[1]]), c(2, T_max, draws, P, P))
})


test_that("subset_networks_sequential zero-fills all current-trial time slices when current time is zero", {
    K <- 3
    T_max <- 4
    draws <- 3
    P <- 3

    net <- array(
        seq_len(K * T_max * draws * P * P),
        dim = c(K, T_max, draws, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    event_data <- data.frame(
        trial_numeric = c(1, 1, 2, 2),
        discrete_time = c(0, 1, 0, 1)
    )

    out <- subset_networks_sequential(
        i = 3,
        network_data = net,
        event_data = event_data
    )

    expected <- net[1:2, , , , , drop = FALSE]
    expected[2, 2:4, , , ] <- 0

    expect_equal(out[[1]], expected)
})


test_that("subset_networks_sequential works with a list of 3D, 4D, and 5D arrays", {
    K <- 4
    T_max <- 5
    draws <- 3
    P <- 3

    net_3d <- array(
        seq_len(draws * P * P),
        dim = c(draws, P, P),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_4d <- array(
        seq_len(K * draws * P * P),
        dim = c(K, draws, P, P),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_5d <- array(
        seq_len(K * T_max * draws * P * P),
        dim = c(K, T_max, draws, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    event_data <- data.frame(
        trial_numeric = c(1, 1, 2, 2, 3),
        discrete_time = c(1, 2, 1, 2, 1)
    )

    out <- subset_networks_sequential(
        i = 4,
        network_data = list(net_3d, net_4d, net_5d),
        event_data = event_data
    )

    expected_5d <- net_5d[1:2, , , , , drop = FALSE]
    expected_5d[2, 4:5, , , ] <- 0

    expect_type(out, "list")
    expect_length(out, 3)

    expect_equal(out[[1]], net_3d)
    expect_equal(out[[2]], net_4d[1:2, , , , drop = FALSE])
    expect_equal(out[[3]], expected_5d)
})


test_that("subset_networks_sequential works when dimensions are in a valid non-standard order", {
    K <- 4
    T_max <- 5
    draws <- 3
    P <- 3

    net <- array(
        seq_len(draws * K * P * T_max * P),
        dim = c(draws, K, P, T_max, P),
        dimnames = list(
            draw = NULL,
            trial = NULL,
            focal_ID = NULL,
            time = NULL,
            other_ID = NULL
        )
    )

    event_data <- data.frame(
        trial_numeric = c(1, 1, 2, 2, 3),
        discrete_time = c(1, 2, 1, 2, 1)
    )

    out <- subset_networks_sequential(
        i = 4,
        network_data = net,
        event_data = event_data
    )

    expected <- net[, 1:2, , , , drop = FALSE]
    expected[, 2, , 4:5, ] <- 0

    expect_equal(out[[1]], expected)
    expect_equal(dim(out[[1]]), c(draws, 2, P, T_max, P))
    expect_equal(names(dimnames(out[[1]])), c("draw", "trial", "focal_ID", "time", "other_ID"))
})


test_that("subset_networks_sequential rejects unsupported network_data inputs", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 2),
        discrete_time = c(1, 2, 1)
    )

    expect_error(
        subset_networks_sequential(
            i = 2,
            network_data = "not a network",
            event_data = event_data
        ),
        "network_data must be a data.frame, array, or list of arrays",
        fixed = TRUE
    )
})


test_that("subset_networks_sequential rejects lists containing non-arrays", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 2),
        discrete_time = c(1, 2, 1)
    )

    net <- array(
        0,
        dim = c(2, 3, 3),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    expect_error(
        subset_networks_sequential(
            i = 2,
            network_data = list(net, data.frame(x = 1)),
            event_data = event_data
        ),
        "all list elements in network_data must be arrays",
        fixed = TRUE
    )
})


test_that("subset_networks_sequential rejects arrays without named dimensions", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 2),
        discrete_time = c(1, 2, 1)
    )

    net <- array(0, dim = c(2, 3, 3))

    expect_error(
        subset_networks_sequential(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "network arrays must have named dimensions",
        fixed = TRUE
    )
})


test_that("subset_networks_sequential rejects arrays with unknown dimensions", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 2),
        discrete_time = c(1, 2, 1)
    )

    net <- array(
        0,
        dim = c(2, 3, 3, 4),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL,
            wave = NULL
        )
    )

    expect_error(
        subset_networks_sequential(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "unknown network array dimension\\(s\\): wave"
    )
})


test_that("subset_networks_sequential rejects arrays missing draw dimension", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 2),
        discrete_time = c(1, 2, 1)
    )

    net <- array(
        0,
        dim = c(2, 3, 3),
        dimnames = list(
            trial = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    expect_error(
        subset_networks_sequential(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "network array is missing required dimension\\(s\\): draw"
    )
})


test_that("subset_networks_sequential rejects arrays missing focal_ID dimension", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 2),
        discrete_time = c(1, 2, 1)
    )

    net <- array(
        0,
        dim = c(2, 2, 3),
        dimnames = list(
            draw = NULL,
            trial = NULL,
            other_ID = NULL
        )
    )

    expect_error(
        subset_networks_sequential(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "network array is missing required dimension\\(s\\): focal_ID"
    )
})


test_that("subset_networks_sequential rejects arrays missing other_ID dimension", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 2),
        discrete_time = c(1, 2, 1)
    )

    net <- array(
        0,
        dim = c(2, 2, 3),
        dimnames = list(
            draw = NULL,
            trial = NULL,
            focal_ID = NULL
        )
    )

    expect_error(
        subset_networks_sequential(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "network array is missing required dimension\\(s\\): other_ID"
    )
})


test_that("subset_networks_sequential reports multiple missing required dimensions together", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 2),
        discrete_time = c(1, 2, 1)
    )

    net <- array(
        0,
        dim = c(2, 3),
        dimnames = list(
            trial = NULL,
            time = NULL
        )
    )

    expect_error(
        subset_networks_sequential(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "network array is missing required dimension\\(s\\): draw, focal_ID, other_ID"
    )
})
