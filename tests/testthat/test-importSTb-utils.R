test_that("array_to_edge_matrix extracts 3D posterior arrays in focal-other edge order", {
    P <- 3
    draws <- 4

    net_obj <- array(
        NA_real_,
        dim = c(draws, P, P),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    for (d in seq_len(draws)) {
        for (i_focal in seq_len(P)) {
            for (i_other in seq_len(P)) {
                net_obj[d, i_focal, i_other] <- 100 * d + 10 * i_focal + i_other
            }
        }
    }

    out <- array_to_edge_matrix(net_obj)

    expected <- cbind(
        net_obj[, 1, 2],
        net_obj[, 1, 3],
        net_obj[, 2, 1],
        net_obj[, 2, 3],
        net_obj[, 3, 1],
        net_obj[, 3, 2]
    )

    expect_type(out, "double")
    expect_equal(dim(out), c(draws, P * (P - 1)))
    expect_equal(out, expected)
})


test_that("array_to_edge_matrix extracts 4D posterior arrays for the requested trial", {
    K <- 3
    P <- 3
    draws <- 4
    trial <- 2

    net_obj <- array(
        NA_real_,
        dim = c(K, draws, P, P),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    for (k in seq_len(K)) {
        for (d in seq_len(draws)) {
            for (i_focal in seq_len(P)) {
                for (i_other in seq_len(P)) {
                    net_obj[k, d, i_focal, i_other] <-
                        1000 * k + 100 * d + 10 * i_focal + i_other
                }
            }
        }
    }

    out <- array_to_edge_matrix(net_obj, trial = trial)

    expected <- cbind(
        net_obj[trial, , 1, 2],
        net_obj[trial, , 1, 3],
        net_obj[trial, , 2, 1],
        net_obj[trial, , 2, 3],
        net_obj[trial, , 3, 1],
        net_obj[trial, , 3, 2]
    )

    expect_equal(dim(out), c(draws, P * (P - 1)))
    expect_equal(out, expected)
})


test_that("array_to_edge_matrix extracts 5D posterior arrays for requested trial and time", {
    K <- 3
    T_max <- 2
    P <- 3
    draws <- 4
    trial <- 3
    time <- 2

    net_obj <- array(
        NA_real_,
        dim = c(K, T_max, draws, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    for (k in seq_len(K)) {
        for (tt in seq_len(T_max)) {
            for (d in seq_len(draws)) {
                for (i_focal in seq_len(P)) {
                    for (i_other in seq_len(P)) {
                        net_obj[k, tt, d, i_focal, i_other] <-
                            10000 * k + 1000 * tt + 100 * d + 10 * i_focal + i_other
                    }
                }
            }
        }
    }

    out <- array_to_edge_matrix(net_obj, trial = trial, time = time)

    expected <- cbind(
        net_obj[trial, time, , 1, 2],
        net_obj[trial, time, , 1, 3],
        net_obj[trial, time, , 2, 1],
        net_obj[trial, time, , 2, 3],
        net_obj[trial, time, , 3, 1],
        net_obj[trial, time, , 3, 2]
    )

    expect_equal(dim(out), c(draws, P * (P - 1)))
    expect_equal(out, expected)
})


test_that("array_to_edge_matrix rejects arrays that are not 3D, 4D, or 5D", {
    bad_2d <- array(1, dim = c(2, 3))
    bad_6d <- array(1, dim = c(1, 2, 3, 4, 5, 6))

    expect_error(
        array_to_edge_matrix(bad_2d),
        "posterior network arrays must be 3D, 4D, or 5D",
        fixed = TRUE
    )

    expect_error(
        array_to_edge_matrix(bad_6d),
        "posterior network arrays must be 3D, 4D, or 5D",
        fixed = TRUE
    )
})


test_that("validate_posterior_array ignores non-array network objects", {
    expect_silent(
        validate_posterior_array(
            net_obj = list(a = 1),
            network_index = 1,
            K = 3,
            T_max = 4,
            P = 5
        )
    )
})


test_that("validate_posterior_array accepts valid 3D, 4D, and 5D posterior arrays", {
    K <- 3
    T_max <- 4
    P <- 5
    draws <- 6

    net_3d <- array(
        0,
        dim = c(draws, P, P),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_4d <- array(
        0,
        dim = c(K, draws, P, P),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_5d <- array(
        0,
        dim = c(K, T_max, draws, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    expect_silent(validate_posterior_array(net_3d, 1, K, T_max, P))
    expect_silent(validate_posterior_array(net_4d, 2, K, T_max, P))
    expect_silent(validate_posterior_array(net_5d, 3, K, T_max, P))
})


test_that("validate_posterior_array rejects arrays with invalid dimensionality", {
    bad_2d <- array(
        0,
        dim = c(3, 3),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL
        )
    )

    bad_6d <- array(
        0,
        dim = c(2, 3, 4, 5, 6, 7)
    )

    expect_error(
        validate_posterior_array(bad_2d, 1, K = 2, T_max = 3, P = 3),
        "must be 3D \\[draw, focal_ID, other_ID\\]"
    )

    expect_error(
        validate_posterior_array(bad_6d, 1, K = 2, T_max = 3, P = 3),
        "must be 3D \\[draw, focal_ID, other_ID\\]"
    )
})


test_that("validate_posterior_array requires correct dimension names for 3D arrays", {
    P <- 3
    draws <- 4

    unnamed <- array(0, dim = c(draws, P, P))

    wrong_order <- array(
        0,
        dim = c(draws, P, P),
        dimnames = list(
            focal_ID = NULL,
            draw = NULL,
            other_ID = NULL
        )
    )

    wrong_names <- array(
        0,
        dim = c(draws, P, P),
        dimnames = list(
            draw = NULL,
            focal = NULL,
            other = NULL
        )
    )

    expect_error(
        validate_posterior_array(unnamed, 1, K = 2, T_max = 3, P = P),
        "must have named dimensions: draw, focal_ID, other_ID",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(wrong_order, 1, K = 2, T_max = 3, P = P),
        "must have named dimensions: draw, focal_ID, other_ID",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(wrong_names, 1, K = 2, T_max = 3, P = P),
        "must have named dimensions: draw, focal_ID, other_ID",
        fixed = TRUE
    )
})


test_that("validate_posterior_array requires correct dimension names for 4D arrays", {
    K <- 2
    P <- 3
    draws <- 4

    wrong_order <- array(
        0,
        dim = c(K, draws, P, P),
        dimnames = list(
            draw = NULL,
            trial = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    expect_error(
        validate_posterior_array(wrong_order, 7, K = K, T_max = 3, P = P),
        "Posterior network array 7 must have named dimensions: trial, draw, focal_ID, other_ID",
        fixed = TRUE
    )
})


test_that("validate_posterior_array requires correct dimension names for 5D arrays", {
    K <- 2
    T_max <- 3
    P <- 4
    draws <- 5

    wrong_order <- array(
        0,
        dim = c(K, T_max, draws, P, P),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            time = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    expect_error(
        validate_posterior_array(wrong_order, 8, K = K, T_max = T_max, P = P),
        "Posterior network array 8 must have named dimensions: trial, time, draw, focal_ID, other_ID",
        fixed = TRUE
    )
})


test_that("validate_posterior_array checks non-draw dimensions for 3D arrays", {
    draws <- 4
    P <- 3

    wrong_focal <- array(
        0,
        dim = c(draws, P + 1, P),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    wrong_other <- array(
        0,
        dim = c(draws, P, P + 1),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    expect_error(
        validate_posterior_array(wrong_focal, 1, K = 2, T_max = 3, P = P),
        "has dimensions 4 x 4 x 3, but expected draws x 3 x 3",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(wrong_other, 1, K = 2, T_max = 3, P = P),
        "has dimensions 4 x 3 x 4, but expected draws x 3 x 3",
        fixed = TRUE
    )
})


test_that("validate_posterior_array checks non-draw dimensions for 4D arrays", {
    K <- 2
    draws <- 4
    P <- 3

    wrong_trial <- array(
        0,
        dim = c(K + 1, draws, P, P),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    wrong_focal <- array(
        0,
        dim = c(K, draws, P + 1, P),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    wrong_other <- array(
        0,
        dim = c(K, draws, P, P + 1),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    expect_error(
        validate_posterior_array(wrong_trial, 1, K = K, T_max = 3, P = P),
        "has dimensions 3 x 4 x 3 x 3, but expected 2 x draws x 3 x 3",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(wrong_focal, 1, K = K, T_max = 3, P = P),
        "has dimensions 2 x 4 x 4 x 3, but expected 2 x draws x 3 x 3",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(wrong_other, 1, K = K, T_max = 3, P = P),
        "has dimensions 2 x 4 x 3 x 4, but expected 2 x draws x 3 x 3",
        fixed = TRUE
    )
})


test_that("validate_posterior_array checks non-draw dimensions for 5D arrays", {
    K <- 2
    T_max <- 3
    draws <- 4
    P <- 5

    wrong_trial <- array(
        0,
        dim = c(K + 1, T_max, draws, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    wrong_time <- array(
        0,
        dim = c(K, T_max + 1, draws, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    wrong_focal <- array(
        0,
        dim = c(K, T_max, draws, P + 1, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    wrong_other <- array(
        0,
        dim = c(K, T_max, draws, P, P + 1),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    expect_error(
        validate_posterior_array(wrong_trial, 1, K = K, T_max = T_max, P = P),
        "has dimensions 3 x 3 x 4 x 5 x 5, but expected 2 x 3 x draws x 5 x 5",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(wrong_time, 1, K = K, T_max = T_max, P = P),
        "has dimensions 2 x 4 x 4 x 5 x 5, but expected 2 x 3 x draws x 5 x 5",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(wrong_focal, 1, K = K, T_max = T_max, P = P),
        "has dimensions 2 x 3 x 4 x 6 x 5, but expected 2 x 3 x draws x 5 x 5",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(wrong_other, 1, K = K, T_max = T_max, P = P),
        "has dimensions 2 x 3 x 4 x 5 x 6, but expected 2 x 3 x draws x 5 x 5",
        fixed = TRUE
    )
})


test_that("validate_posterior_array allows any draw count greater than or equal to 2", {
    K <- 2
    T_max <- 3
    P <- 4

    net_3d <- array(
        0,
        dim = c(2, P, P),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_4d <- array(
        0,
        dim = c(K, 10, P, P),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_5d <- array(
        0,
        dim = c(K, T_max, 100, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    expect_silent(validate_posterior_array(net_3d, 1, K, T_max, P))
    expect_silent(validate_posterior_array(net_4d, 1, K, T_max, P))
    expect_silent(validate_posterior_array(net_5d, 1, K, T_max, P))
})


test_that("validate_posterior_array rejects posterior arrays with fewer than 2 draws", {
    K <- 2
    T_max <- 3
    P <- 4

    net_3d <- array(
        0,
        dim = c(1, P, P),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_4d <- array(
        0,
        dim = c(K, 1, P, P),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_5d <- array(
        0,
        dim = c(K, T_max, 1, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    expect_error(
        validate_posterior_array(net_3d, 1, K, T_max, P),
        "at least 2 draws",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(net_4d, 1, K, T_max, P),
        "at least 2 draws",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(net_5d, 1, K, T_max, P),
        "at least 2 draws",
        fixed = TRUE
    )
})


test_that("validate_posterior_array rejects missing values in posterior arrays", {
    K <- 2
    T_max <- 3
    P <- 4
    draws <- 5

    net_3d <- array(
        0,
        dim = c(draws, P, P),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_4d <- array(
        0,
        dim = c(K, draws, P, P),
        dimnames = list(
            trial = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_5d <- array(
        0,
        dim = c(K, T_max, draws, P, P),
        dimnames = list(
            trial = NULL,
            time = NULL,
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_3d[1, 1, 1] <- NA_real_
    net_4d[1, 1, 1, 1] <- NA_real_
    net_5d[1, 1, 1, 1, 1] <- NA_real_

    expect_error(
        validate_posterior_array(net_3d, 1, K, T_max, P),
        "Missing values detected in posterior network array 1",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(net_4d, 2, K, T_max, P),
        "Missing values detected in posterior network array 2",
        fixed = TRUE
    )

    expect_error(
        validate_posterior_array(net_5d, 3, K, T_max, P),
        "Missing values detected in posterior network array 3",
        fixed = TRUE
    )
})


test_that("validate_posterior_array checks dimensions before missing values", {
    P <- 3

    net_obj <- array(
        0,
        dim = c(4, P + 1, P),
        dimnames = list(
            draw = NULL,
            focal_ID = NULL,
            other_ID = NULL
        )
    )

    net_obj[1, 1, 1] <- NA_real_

    expect_error(
        validate_posterior_array(net_obj, 1, K = 2, T_max = 3, P = P),
        "has dimensions 4 x 4 x 3, but expected draws x 3 x 3",
        fixed = TRUE
    )
})
