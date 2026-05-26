test_that("subset_networks_parallel subsets dynamic data.frame networks by current discrete time", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
    )

    network_data <- data.frame(
        trial_numeric = c(2, 1, 1, 2, 1),
        discrete_time = c(3, 2, 1, 1, 3),
        focal = c("A", "A", "B", "C", "C"),
        other = c("B", "C", "A", "A", "B"),
        weight = seq_len(5)
    )

    out <- subset_networks_parallel(
        i = 2,
        network_data = network_data,
        event_data = event_data
    )

    expect_s3_class(out, "data.frame")
    expect_true(all(out$discrete_time <= 2))
    expect_equal(out$trial_numeric, c(1, 1, 2))
    expect_equal(out$discrete_time, c(1, 2, 1))
})


test_that("subset_networks_parallel sorts static data.frame networks by trial", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
    )

    network_data <- data.frame(
        trial_numeric = c(3, 1, 2, 1),
        focal = c("A", "B", "C", "D"),
        other = c("B", "C", "D", "A"),
        weight = seq_len(4)
    )

    out <- subset_networks_parallel(
        i = 2,
        network_data = network_data,
        event_data = event_data
    )

    expect_equal(out$trial_numeric, c(1, 1, 2, 3))
    expect_false("discrete_time" %in% names(out))
})


test_that("subset_networks_parallel leaves 3D posterior arrays unchanged", {
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
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
    )

    out <- subset_networks_parallel(
        i = 2,
        network_data = net,
        event_data = event_data
    )

    expect_type(out, "list")
    expect_length(out, 1)
    expect_equal(out[[1]], net)
    expect_equal(dim(out[[1]]), dim(net))
    expect_equal(names(dimnames(out[[1]])), c("draw", "focal_ID", "other_ID"))
})


test_that("subset_networks_parallel leaves 4D trial-varying posterior arrays unchanged", {
    K <- 2
    draws <- 4
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
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
    )

    out <- subset_networks_parallel(
        i = 2,
        network_data = net,
        event_data = event_data
    )

    expect_type(out, "list")
    expect_length(out, 1)
    expect_equal(out[[1]], net)
    expect_equal(dim(out[[1]]), dim(net))
    expect_equal(names(dimnames(out[[1]])), c("trial", "draw", "focal_ID", "other_ID"))
})


test_that("subset_networks_parallel subsets 5D trial-time posterior arrays up to current time", {
    K <- 2
    T_max <- 5
    draws <- 4
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
        trial_numeric = c(1, 1, 1, 2),
        discrete_time = c(1, 2, 3, 4)
    )

    out <- subset_networks_parallel(
        i = 3,
        network_data = net,
        event_data = event_data
    )

    expected <- net[, 1:3, , , , drop = FALSE]

    expect_type(out, "list")
    expect_length(out, 1)
    expect_equal(out[[1]], expected)
    expect_equal(dim(out[[1]]), c(K, 3, draws, P, P))
    expect_equal(
        names(dimnames(out[[1]])),
        c("trial", "time", "draw", "focal_ID", "other_ID")
    )
})


test_that("subset_networks_parallel keeps one time slice when current discrete time is zero", {
    K <- 2
    T_max <- 5
    draws <- 4
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
        trial_numeric = c(1, 1, 1),
        discrete_time = c(0, 1, 2)
    )

    out <- subset_networks_parallel(
        i = 1,
        network_data = net,
        event_data = event_data
    )

    expected <- net[, 1, , , , drop = FALSE]

    expect_equal(out[[1]], expected)
    expect_equal(dim(out[[1]])[2], 1)
    expect_equal(names(dimnames(out[[1]])), names(dimnames(net)))
})


test_that("subset_networks_parallel works with a list of 3D, 4D, and 5D arrays", {
    K <- 2
    T_max <- 4
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
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
    )

    out <- subset_networks_parallel(
        i = 2,
        network_data = list(net_3d, net_4d, net_5d),
        event_data = event_data
    )

    expect_type(out, "list")
    expect_length(out, 3)

    expect_equal(out[[1]], net_3d)
    expect_equal(out[[2]], net_4d)
    expect_equal(out[[3]], net_5d[, 1:2, , , , drop = FALSE])

    expect_equal(names(dimnames(out[[1]])), c("draw", "focal_ID", "other_ID"))
    expect_equal(names(dimnames(out[[2]])), c("trial", "draw", "focal_ID", "other_ID"))
    expect_equal(names(dimnames(out[[3]])), c("trial", "time", "draw", "focal_ID", "other_ID"))
})


test_that("subset_networks_parallel rejects unsupported network_data inputs", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
    )

    expect_error(
        subset_networks_parallel(
            i = 2,
            network_data = "not a network",
            event_data = event_data
        ),
        "network_data must be a data.frame, array, or list of arrays",
        fixed = TRUE
    )
})


test_that("subset_networks_parallel rejects lists containing non-arrays", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
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
        subset_networks_parallel(
            i = 2,
            network_data = list(net, data.frame(x = 1)),
            event_data = event_data
        ),
        "all list elements in network_data must be arrays",
        fixed = TRUE
    )
})


test_that("subset_networks_parallel rejects arrays without named dimensions", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
    )

    net <- array(0, dim = c(2, 3, 3))

    expect_error(
        subset_networks_parallel(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "network arrays must have named dimensions",
        fixed = TRUE
    )
})


test_that("subset_networks_parallel rejects arrays with unknown named dimensions", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
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
        subset_networks_parallel(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "unknown network array dimension\\(s\\): wave"
    )
})


test_that("subset_networks_parallel rejects arrays missing draw dimension", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
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
        subset_networks_parallel(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "network array is missing required dimension\\(s\\): draw"
    )
})


test_that("subset_networks_parallel rejects arrays missing focal_ID dimension", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
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
        subset_networks_parallel(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "network array is missing required dimension\\(s\\): focal_ID"
    )
})


test_that("subset_networks_parallel rejects arrays missing other_ID dimension", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
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
        subset_networks_parallel(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "network array is missing required dimension\\(s\\): other_ID"
    )
})


test_that("subset_networks_parallel reports multiple missing required dimensions together", {
    event_data <- data.frame(
        trial_numeric = c(1, 1, 1),
        discrete_time = c(1, 2, 3)
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
        subset_networks_parallel(
            i = 2,
            network_data = net,
            event_data = event_data
        ),
        "network array is missing required dimension\\(s\\): draw, focal_ID, other_ID"
    )
})
