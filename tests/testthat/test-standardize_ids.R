test_that("IDs are properly mapped and filled", {
    networks <- data.frame(
        trial = c(1, 1, 1),
        focal = c("A", "A", "B"),
        other = c("B", "C", "C")
    )

    event_data <- data.frame(
        id = c("A", "B"),
        trial = c(1, 1),
        time = c(0, 1),
        t_end = c(3, 3)
    )

    ILV_c <- data.frame(
        id = c("A", "B", "C"),
        age = c(1, 2, 3)
    )

    ILV_tv <- data.frame(
        id = c("A", "B"),
        trial = c(1, 1),
        time = c(1, 1),
        ilv_var = c(10, 20)
    )

    t_weights <- data.frame(
        id = c("A", "B"),
        trial = c(1, 1),
        time = c(1, 1),
        t_weight = c(0.5, 1.0)
    )

    result <- suppressWarnings(standardize_ids(networks, event_data, ILV_c, ILV_tv, t_weights))

    # check that id_map was created correctly
    expect_equal(result$id_map$id, c("A", "B", "C"))
    expect_equal(result$id_map$id_numeric, 1:3)

    # check that C was added to event_data as censored
    added_c <- result$event_data[result$event_data$id == "C", ]
    expect_equal(nrow(added_c), 1)
    expect_equal(added_c$time, 4)
    expect_equal(added_c$t_end, 3)

    # check that ILV_tv was padded for C
    padded_c <- result$ILV_tv[result$ILV_tv$id == "C", ]
    expect_true(all(padded_c$ilv_var == 0))

    # check that t_weights was padded for C
    padded_tw <- result$t_weights[result$t_weights$id == "C", ]
    expect_true(all(padded_tw$t_weight == 0))

    # check that all outputs include id_numeric
    expect_true("id_numeric" %in% colnames(result$event_data))
    expect_true("id_numeric" %in% colnames(result$ILV_c))
    expect_true("id_numeric" %in% colnames(result$ILV_tv))
    expect_true("id_numeric" %in% colnames(result$t_weights))
})


test_that("Error if event_data contains unknown IDs", {
    networks <- data.frame(trial = 1, from = "A", to = "B")
    event_data <- data.frame(id = "Z", trial = 1, time = 0, t_end = 3)

    expect_error(standardize_ids(networks, event_data),
        regexp = "missing from `networks`"
    )
})

test_that("Error if ILV_c contains unknown IDs", {
    networks <- data.frame(trial = 1, from = "A", to = "B")
    event_data <- data.frame(id = c("A", "B"), trial = 1, time = 0:1, t_end = 3)
    ILV_c <- data.frame(id = c("A", "Z"), age = c(1, 2))

    expect_error(standardize_ids(networks, event_data, ILV_c = ILV_c),
        regexp = "missing from `networks`"
    )
})

test_that("Error if ILV_tv contains unknown IDs", {
    networks <- data.frame(trial = 1, from = "A", to = "B")
    event_data <- data.frame(id = c("A", "B"), trial = 1, time = 0:1, t_end = 3)
    ILV_tv <- data.frame(id = c("A", "Z"), trial = 1, time = 1, ilv_var = 99)

    expect_error(standardize_ids(networks, event_data, ILV_tv = ILV))
})
