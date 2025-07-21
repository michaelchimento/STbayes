test_that("check_trials passes when trials match", {
    event_data <- data.frame(trial = c(1, 2, 2))
    networks <- data.frame(trial = c(1, 2), from = c("A", "B"), to = c("B", "A"))

    expect_silent(check_trials(event_data, networks))
})

test_that("check_trials errors when trial is missing in networks", {
    event_data <- data.frame(trial = c(1, 2, 3))
    networks <- data.frame(trial = c(1, 2), from = c("A", "B"), to = c("B", "A"))

    expect_error(
        check_trials(event_data, networks)
    )
})

test_that("check_trials errors when trial is missing in event_data", {
    event_data <- data.frame(trial = c(1, 2))
    networks <- data.frame(trial = c(1, 2, 3), from = c("A", "B", "C"), to = c("B", "A", "A"))

    expect_error(
        check_trials(event_data, networks)
    )
})

test_that("check_trials warns for Bayesian input with multiple trials", {
    mock_net <- structure(list(), class = "bison_model") # mimic bison_model
    event_data <- data.frame(trial = c(1, 2, 2))

    expect_warning(
        check_trials(event_data, mock_net)
    )
})
