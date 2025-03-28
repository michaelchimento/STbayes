test_that("complicated data gets imported correctly (pigeons)", {
    data_list = import_user_STb(event_data = STbayes::pigeon_events,
                                networks = STbayes::pigeon_networks,
                                network_type = "directed",
                                t_weights = STbayes::pigeon_tweights,
                                ILV_tv = STbayes::pigeon_ILVtv,
                                ILVi = c("avf_mon"))
    expect_identical(data_list$K, STbayes::pigeon_data_list_TRUE$K)
    expect_identical(data_list$P, STbayes::pigeon_data_list_TRUE$Z)
    expect_equivalent(data_list$N, STbayes::pigeon_data_list_TRUE$N)
    expect_equivalent(data_list$N_c, STbayes::pigeon_data_list_TRUE$N_c)
    expect_equivalent(data_list$T, STbayes::pigeon_data_list_TRUE$T)
    #expect_equivalent(data_list$t, STbayes::pigeon_data_list_TRUE$t)
    expect_equivalent(data_list$Q, STbayes::pigeon_data_list_TRUE$Q)
    expect_equivalent(data_list$D, STbayes::pigeon_data_list_TRUE$D*100)
    #expect_equivalent(data_list$ind_id, STbayes::pigeon_data_list_TRUE$bird_id)
    expect_equivalent(data_list$Z, STbayes::pigeon_data_list_TRUE$C_fovmon)
    expect_equivalent(data_list$W, STbayes::pigeon_data_list_TRUE$C_fovmon)
    expect_equivalent(data_list$ILV_avf_mon, STbayes::pigeon_data_list_TRUE$M_avf)
    expect_equivalent(data_list$A_avf, STbayes::pigeon_data_list_TRUE$A_avf)
    expect_equivalent(data_list$A_dist, STbayes::pigeon_data_list_TRUE$A_dist)
    expect_equivalent(data_list$A_vor, STbayes::pigeon_data_list_TRUE$A_vor)
    expect_equivalent(data_list$A_fov, STbayes::pigeon_data_list_TRUE$A_fov)
})
