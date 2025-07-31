test_that("detect_ILV_datatype correctly detects boolean vectors", {
    expect_equal(detect_ILV_datatype(c(TRUE, FALSE, NA)), "boolean")
})

test_that("detect_ILV_datatype correctly detects categorical variables", {
    expect_equal(detect_ILV_datatype(factor(c("a", "b", "a"))), "categorical")
    expect_equal(detect_ILV_datatype(c("a", "b", "a")), "categorical") # integers ≥ 0
})

test_that("detect_ILV_datatype correctly detects continuous variables", {
    expect_equal(detect_ILV_datatype(c(1.1, 2.5, 3.7, NA)), "continuous")
    expect_equal(detect_ILV_datatype(c(-1, 0, 1, NA)), "continuous")
})

test_that("detect_ILV_datatype throws error for unsupported input", {
    expect_error(detect_ILV_datatype(list(1, 2, 3)), "Please input boolean")
})
