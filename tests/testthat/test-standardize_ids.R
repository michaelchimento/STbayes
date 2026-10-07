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

test_that("standardize_ids extracts character IDs from array dimnames", {
    bird_ids <- c("0700ED8BF7", "0700ED8C01", "0700ED8C02")

    net_array <- array(
        rnorm(100 * 3 * 3),
        dim = c(100, 3, 3),
        dimnames = list(draw = NULL, focal_ID = bird_ids, other_ID = bird_ids)
    )

    event_data <- data.frame(
        id = bird_ids, trial = 1, time = c(0, 1, 2), t_end = 3,
        stringsAsFactors = FALSE
    )

    result <- standardize_ids(list(net_array), event_data)

    expect_equal(sort(result$id_map$id), sort(bird_ids))
    expect_equal(nrow(result$id_map), 3)
    expect_true(all(event_data$id %in% result$id_map$id))
    expect_equal(result$id_map$id[result$id_map$id_numeric == 1], bird_ids[1])
    expect_equal(result$id_map$id[result$id_map$id_numeric == 2], bird_ids[2])
    expect_equal(result$id_map$id[result$id_map$id_numeric == 3], bird_ids[3])
})

test_that("standardize_ids falls back to 1:P for arrays with NULL dimnames values", {
    net_array <- array(
        rnorm(100 * 3 * 3),
        dim = c(100, 3, 3),
        dimnames = list(draw = NULL, focal_ID = NULL, other_ID = NULL)
    )

    event_data <- data.frame(
        id = as.character(1:3), trial = 1, time = c(0, 1, 2), t_end = 3,
        stringsAsFactors = FALSE
    )

    result <- standardize_ids(list(net_array), event_data)

    expect_equal(result$id_map$id, c("1", "2", "3"))
    expect_equal(result$id_map$id_numeric, 1:3)
})

test_that("standardize_ids preserves positional order for character IDs", {
    ids <- c("Zebra", "Ant", "Mole")

    net_array <- array(
        rnorm(50 * 3 * 3),
        dim = c(50, 3, 3),
        dimnames = list(draw = NULL, focal_ID = ids, other_ID = ids)
    )

    event_data <- data.frame(
        id = ids, trial = 1, time = c(0, 1, 2), t_end = 3,
        stringsAsFactors = FALSE
    )

    result <- standardize_ids(list(net_array), event_data)

    expect_equal(result$id_map$id_numeric[result$id_map$id == "Zebra"], 1)
    expect_equal(result$id_map$id_numeric[result$id_map$id == "Ant"], 2)
    expect_equal(result$id_map$id_numeric[result$id_map$id == "Mole"], 3)
})
