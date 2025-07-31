test_that("generate_X_matrix works with categorical input and reference level", {
    x <- c("a", "b", "a", "c")
    x <- as.numeric(as.factor(x))
    result <- generate_X_matrix(x, ilv_name = "foo", n_levels = 3)

    # expect 2 columns (b and c)
    expect_equal(ncol(result), 2)
    expect_equal(colnames(result), c("foo_2", "foo_3"))

    # first two rows correspond to levels a and b
    expect_equivalent(result[1, ], c(0, 0)) # a (ref)
    expect_equivalent(result[2, ], c(1, 0)) # b
})

test_that("generate_X_matrix handles single-level input", {
    x <- as.numeric(as.factor(rep("only", 4)))
    result <- generate_X_matrix(x, ilv_name = "single", n_levels = 1)
    expect_equal(dim(result), c(4, 0)) # model.matrix drops intercept and nothing else remains
})
