test_that("process_ILVs produces identical output to original ILV handling", {
    # Mock inputs
    STb_data <- list(
        ILV_age = matrix(0.5, nrow = 2, ncol = 3), # time-varying
        ILV_sex = c(1, 2, 1), # constant
        ILV_weight = rep(0.25, 3), # constant
        ILV_timevarying = c(TRUE, FALSE, FALSE),
        ILV_datatypes = c("continuous", "categorical", "continuous"),
        ILV_n_levels = c(2)
    )

    names(STb_data$ILV_timevarying) <- c("ILV_age", "ILV_sex", "ILV_weight")
    names(STb_data$ILV_datatypes) <- c("ILV_age", "ILV_sex", "ILV_weight")
    names(STb_data$ILV_n_levels) <- c("ILV_sex")

    veff_params <- c("age", "sex") # only age and sex have varying effects
    count_start <- 1

    # Run refactored version
    result_ILVi <- process_ILVs(
        ilv_vars = c("age"),
        ilv_vars_clean = c("age"),
        ilv_datatypes = STb_data$ILV_datatypes,
        ilv_n_levels = STb_data$ILV_n_levels,
        ilv_timevarying = STb_data$ILV_timevarying,
        veff_params = veff_params,
        veff_type = "id",
        suffix = "i",
        STb_data = STb_data,
        count_start = count_start,
        prior_beta = "normal(0, 1)"
    )

    result_ILVi$prior

    result_ILVs <- process_ILVs(
        ilv_vars = c("sex"),
        ilv_vars_clean = c("sex"),
        ilv_datatypes = STb_data$ILV_datatypes,
        ilv_n_levels = STb_data$ILV_n_levels,
        ilv_timevarying = STb_data$ILV_timevarying,
        veff_params = veff_params,
        veff_type = "id",
        suffix = "s",
        STb_data = STb_data,
        count_start = result_ILVi$count,
        prior_beta = "normal(0, 1)"
    )

    result_ILVm <- process_ILVs(
        ilv_vars = c("weight"),
        ilv_vars_clean = c("weight"),
        ilv_datatypes = STb_data$ILV_datatypes,
        ilv_n_levels = STb_data$ILV_n_levels,
        ilv_timevarying = STb_data$ILV_timevarying,
        veff_params = veff_params,
        veff_type = "id",
        suffix = "m",
        STb_data = STb_data,
        count_start = result_ILVs$count,
        prior_beta = "normal(0, 1)"
    )

    # Hardcoded expected outputs
    expect_equal(result_ILVi$param, "real beta_ILVi_age;")
    expect_equal(result_ILVi$prior, "beta_ILVi_age ~ normal(0, 1);")
    expect_equal(result_ILVi$transformed_decl, "array[K, T_max] vector[P] age_i;")
    expect_equal(result_ILVi$transformed_calc, "for (trial in 1:K) for (timestep in 1:T_max) age_i[trial][timestep] = ILV_age[trial][timestep] * beta_ILVi_age + v_id[,1];")
    expect_equal(result_ILVi$term, "exp(age_i[trial,time_step,id])")

    expect_equal(result_ILVs$param, "vector[1] beta_ILVs_sex;")
    expect_equal(result_ILVs$prior, "beta_ILVs_sex ~ normal(0, 1);")
    expect_equal(result_ILVs$transformed_decl, "vector[P] sex_s;")
    expect_equal(result_ILVs$transformed_calc, "sex_s = ILV_sex * beta_ILVs_sex + v_id[,2];")
    expect_equal(result_ILVs$term, "* exp(sex_s[id])")

    expect_equal(result_ILVm$param, "real beta_ILVm_weight;")
    expect_equal(result_ILVm$prior, "beta_ILVm_weight ~ normal(0, 1);")
    expect_equal(result_ILVm$transformed_decl, "vector[P] weight_m;")
    expect_equal(result_ILVm$transformed_calc, "weight_m = ILV_weight * beta_ILVm_weight;")
    expect_equal(result_ILVm$term, "exp(weight_m[id]) *")
})
