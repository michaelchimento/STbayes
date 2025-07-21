test_that("complicated data gets imported correctly (pigeons)", {
    # Load test data from testdata/
    load(test_path("testdata", "pigeon_events.rda")) # loads pigeon_events
    pigeon_events$trial <- as.integer(as.factor(pigeon_events$trial))
    load(test_path("testdata", "pigeon_networks.rda")) # loads pigeon_networks
    load(test_path("testdata", "pigeon_tweights.rda")) # loads pigeon_tweights
    load(test_path("testdata", "pigeon_ILVtv.rda")) # loads pigeon_ILVtv
    load(test_path("testdata", "pigeon_data_list_TRUE.rda")) # loads pigeon_ILVtv

    # Use in import_user_STb
    data_list <- import_user_STb(
        event_data = pigeon_events,
        networks = pigeon_networks,
        network_type = "directed",
        t_weights = pigeon_tweights,
        ILV_tv = pigeon_ILVtv,
        ILVi = c("avf_mon")
    )
    expect_identical(data_list$K, pigeon_data_list_TRUE$K)
    expect_identical(data_list$P, pigeon_data_list_TRUE$Z)
    expect_equivalent(data_list$N, pigeon_data_list_TRUE$N)
    expect_equivalent(data_list$N_c, pigeon_data_list_TRUE$N_c)
    expect_equivalent(data_list$T, pigeon_data_list_TRUE$T)
    # expect_equivalent(data_list$t, pigeon_data_list_TRUE$t)
    expect_equivalent(data_list$Q, pigeon_data_list_TRUE$Q)
    expect_equivalent(data_list$D, pigeon_data_list_TRUE$D * 100)
    # expect_equivalent(data_list$ind_id, pigeon_data_list_TRUE$bird_id)
    expect_equivalent(data_list$Z, pigeon_data_list_TRUE$C_fovmon)
    expect_equivalent(data_list$W, pigeon_data_list_TRUE$C_fovmon)
    expect_equivalent(data_list$ILV_avf_mon, pigeon_data_list_TRUE$M_avf)
    expect_equivalent(data_list$A[1, , , , ], pigeon_data_list_TRUE$A_avf)
    expect_equivalent(data_list$A[2, , , , ], pigeon_data_list_TRUE$A_dist)
    expect_equivalent(data_list$A[3, , , , ], pigeon_data_list_TRUE$A_vor)
    expect_equivalent(data_list$A[4, , , , ], pigeon_data_list_TRUE$A_fov)
})
