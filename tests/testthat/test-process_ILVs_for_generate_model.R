test_that("process_ILVs produces identical output to original ILV handling", {
    # Mock inputs
    STb_data <- list(
        ILV_age = matrix(0.5, nrow = 2, ncol = 3),  # time-varying
        ILV_sex = rep(1, 3),                        # constant
        ILV_weight = rep(0.25, 3)                   # constant
    )

    veff_ID <- c("age", "sex")  # only age and sex have varying effects
    count_start <- 1

    # Run refactored version
    result_ILVi <- process_ILVs(ilv_vars = c("age"),
                                ilv_vars_clean = c("age"),
                                veff_ID = veff_ID,
                                suffix = "i",
                                STb_data = STb_data,
                                count_start = count_start,
                                prior_beta = "normal(0,1)")

    result_ILVs <- process_ILVs(ilv_vars = c("sex"),
                                ilv_vars_clean = c("sex"),
                                veff_ID = veff_ID,
                                suffix = "s",
                                STb_data = STb_data,
                                count_start = result_ILVi$count,
                                prior_beta = "normal(0,1)")

    result_ILVm <- process_ILVs(ilv_vars = c("weight"),
                                ilv_vars_clean = c("weight"),
                                veff_ID = veff_ID,
                                suffix = "m",
                                STb_data = STb_data,
                                count_start = result_ILVs$count,
                                prior_beta = "normal(0,1)")

    # Hardcoded expected outputs
    expect_equal(result_ILVi$param, "real beta_ILVi_age;")
    expect_equal(result_ILVi$prior, "beta_ILVi_age ~ normal(0, 1);")
    expect_true(any(grepl("vector\\[P\\] age_i =", result_ILVi$transformed)))
    expect_true(grepl("exp\\(age_i\\[id\\] \\* ILV_age\\[trial,time_step,id\\]\\)", result_ILVi$term))

    expect_equal(result_ILVs$param, "real beta_ILVs_sex;")
    expect_equal(result_ILVs$prior, "beta_ILVs_sex ~ normal(0, 1);")
    expect_true(any(grepl("vector\\[P\\] sex_s =", result_ILVs$transformed)))
    expect_true(grepl("\\* exp\\(sex_s\\[id\\] \\* ILV_sex\\[id\\]\\)", result_ILVs$term))

    expect_equal(result_ILVm$param, "real beta_ILVm_weight;")
    expect_equal(result_ILVm$prior, "beta_ILVm_weight ~ normal(0, 1);")
    expect_equal(result_ILVm$transformed, NULL)  # not in veff_ID
    expect_true(grepl("exp\\(beta_ILVm_weight \\* ILV_weight\\[id\\]\\) \\*", result_ILVm$term))
})
