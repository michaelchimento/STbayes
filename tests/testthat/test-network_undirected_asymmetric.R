# Regression tests for undirected networks supplied with asymmetric
# (one-directionally stored, zero-padded reciprocal) edge lists.
#
# Bug: such networks were silently collapsed toward all-zero because the
# is_symmetric row-count heuristic + a destructive C++ mirror let zero-padded
# reciprocal rows clobber the real weights. See fill_array.cpp (max mirror) and
# import_user_STb.R (always-symmetrize + asymmetry warning).

# Build a small single-trial diffusion with 4 individuals. Intended undirected
# weights (symmetric): A-B = 2, A-C = 1, B-D = 3.
make_event_data <- function() {
    data.frame(
        id = c("A", "B", "C", "D"),
        trial = 1,
        time = c(0, 1, 2, 3),
        t_end = 4,
        stringsAsFactors = FALSE
    )
}

# Expected symmetric adjacency (sorted id order A=1..D=4), diagonal zeroed.
expected_A <- function() {
    E <- matrix(0, 4, 4)
    E[1, 2] <- E[2, 1] <- 2
    E[1, 3] <- E[3, 1] <- 1
    E[2, 4] <- E[4, 2] <- 3
    E
}

# Full ordered-pair edge list with weights stored in ONE direction only and
# every reciprocal / non-edge padded with assoc = 0 (Sonja's storage shape).
make_asymmetric_edges <- function(include_self_loops) {
    ids <- c("A", "B", "C", "D")
    grid <- expand.grid(focal = ids, other = ids, stringsAsFactors = FALSE)
    if (!include_self_loops) {
        grid <- grid[grid$focal != grid$other, ]
    }
    grid$trial <- 1
    grid$assoc <- 0
    set_w <- function(f, o, w) grid$assoc[grid$focal == f & grid$other == o] <<- w
    set_w("A", "B", 2)
    set_w("A", "C", 1)
    set_w("B", "D", 3)
    rownames(grid) <- NULL
    grid[, c("trial", "focal", "other", "assoc")]
}

test_that("undirected network with self-loops + zero-padded reciprocals is symmetrized (mirror branch)", {
    event_data <- make_event_data()
    el <- make_asymmetric_edges(include_self_loops = TRUE)

    expect_warning(
        data_list <- suppressMessages(
            import_user_STb(event_data, el, network_type = "undirected")
        ),
        regexp = "asymmetric"
    )

    A <- data_list$A
    expect_gt(sum(A), 0)

    slice <- matrix(A[1, 1, 1, , ], 4, 4)
    expect_equal(slice, t(slice)) # symmetric
    expect_equal(slice, expected_A()) # weights == pmax(W, t(W))
})

test_that("undirected network without self-loops (P*(P-1) rows) is still symmetrized (previously-skipped branch)", {
    event_data <- make_event_data()
    el <- make_asymmetric_edges(include_self_loops = FALSE)

    # This is exactly P*(P-1) rows, which the old heuristic misread as
    # "already symmetric" and skipped mirroring on.
    expect_equal(nrow(el), 4 * (4 - 1))

    expect_warning(
        data_list <- suppressMessages(
            import_user_STb(event_data, el, network_type = "undirected")
        ),
        regexp = "asymmetric"
    )

    A <- data_list$A
    expect_gt(sum(A), 0)

    slice <- matrix(A[1, 1, 1, , ], 4, 4)
    expect_equal(slice, t(slice))
    expect_equal(slice, expected_A())
})

test_that("truly symmetric undirected input does not warn and is unchanged", {
    event_data <- make_event_data()
    el <- make_asymmetric_edges(include_self_loops = FALSE)
    # symmetrize the stored weights so the input is genuinely symmetric
    key <- paste(el$focal, el$other)
    recip <- setNames(el$assoc, paste(el$other, el$focal))[key]
    recip[is.na(recip)] <- 0
    el$assoc <- pmax(el$assoc, recip)

    expect_no_warning(
        data_list <- suppressMessages(
            import_user_STb(event_data, el, network_type = "undirected")
        )
    )

    slice <- matrix(data_list$A[1, 1, 1, , ], 4, 4)
    expect_equal(slice, expected_A())
})
