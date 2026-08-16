test_that("it returns 0 when there are no v_id matches", {
    model_text <- "some random model code without v_id"
    expect_equal(return_N_veff(model_text), 0)
})

test_that("it correctly detects a single v_id index", {
    model_text <- "something v_id[,3] something else"
    expect_equal(return_N_veff(model_text), 3)
})

test_that("it returns the max index when multiple v_ids are present", {
    model_text <- "v_id[,2] more text v_id[,5] v_id[,3]"
    expect_equal(return_N_veff(model_text), 5)
})

test_that("it handles messy spacing and newlines", {
    model_text <- "v_id[,1]\n v_id[, 4 ] \n more v_id[,7]"
    expect_equal(return_N_veff(model_text), 7)
})

test_that("it handles empty input", {
    model_text <- ""
    expect_equal(return_N_veff(model_text), 0)
})

test_that("shared v_id[,n] resolves to num_networks (single network = 1)", {
    model_text <- "target += ... v_id[,n] ..."
    expect_equal(return_N_veff(model_text, num_networks = 1), 1)
})

test_that("shared v_id[,n] resolves to num_networks (three networks = 3)", {
    model_text <- "target += ... v_id[,n] ..."
    expect_equal(return_N_veff(model_text, num_networks = 3), 3)
})

test_that("shared v_id[,n] plus literal v_id[,4] takes the max", {
    model_text <- "v_id[,n] ... v_id[,4]"
    expect_equal(return_N_veff(model_text, num_networks = 3), 4)
})

test_that("offset form v_id[,n+4-1] resolves to num_networks + d - 1", {
    model_text <- "v_id[,n+4-1]"
    expect_equal(return_N_veff(model_text, num_networks = 3), 6)
})

test_that("it handles v_trial loop-variable and offset forms", {
    expect_equal(return_N_veff("v_trial[trial,n]", num_networks = 2), 2)
    expect_equal(return_N_veff("v_trial[trial,n+2-1]", num_networks = 3), 4)
})
