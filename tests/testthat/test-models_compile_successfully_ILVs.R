library(dplyr)
# parameters
N <- 30
Tmax <- 30 # number of inter-event intervals

# create all (id, time) combinations
event_data <- data.frame(
    trial = 1,
    id = 1:N,
    time = 1:N,
    t_end = Tmax
)

# create all (id, time) combinations
networks <- data.frame(
    trial = 1,
    focal = 1:N,
    other = sample(1:N),
    weight = 1
)

# create all (id, time) combinations
ILVtimevarying <- expand.grid(
    trial = 1,
    id = 1:N,
    time = 1:Tmax
)

# simulate 3 types of time-varying ILVs
ILV_tv <- ILVtimevarying %>%
    mutate(
        bool_ILV = as.logical(rbinom(n(), 1, 0.5)), # binary
        cont_ILV = rnorm(n(), mean = 0, sd = 1), # continuous
        cat_ILV = factor(sample(1:4, size = n(), replace = TRUE)) # categorical (4 levels)
    )

# create all (id, time) combinations
ILVconstant <- expand.grid(
    trial = 1,
    id = 1:N
)

# simulate 3 types of time-varying ILVs
ILV_c <- ILVconstant %>%
    mutate(
        bool_ILV = as.logical(rbinom(n(), 1, 0.5)), # binary
        cont_ILV = rnorm(n(), mean = 0, sd = 1), # continuous
        cat_ILV = factor(sample(1:4, size = n(), replace = TRUE)) # categorical (4 levels)
    )

#### BOOLEAN ILV ####
test_that("Constant boolean ILV cTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks,
        ILV_c = ILV_c,
        ILVi = c("bool_ILV"),
        ILVs = c("cont_ILV"),
        ILVm = c("cat_ILV")
    )

    stan_code <- generate_STb_model(data_list)

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})
test_that("Timevarying boolean ILV cTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks,
        ILV_tv = ILV_tv,
        ILVi = c("bool_ILV")
    )

    stan_code <- generate_STb_model(data_list)

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("Constant ILV cTADA veffs compile", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks,
        ILV_c = ILV_c,
        ILVi = c("bool_ILV")
    )

    # id type
    stan_code <- generate_STb_model(data_list,
        veff_params = c("bool_ILV"),
        veff_type = "id"
    )
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )

    # trial type
    stan_code <- generate_STb_model(data_list,
        veff_params = c("bool_ILV"),
        veff_type = "trial"
    )
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )

    # both type
    stan_code <- generate_STb_model(data_list,
        veff_params = c("bool_ILV"),
        veff_type = c("id", "trial")
    )
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("Timevarying Boolean ILV cTADA veffs compile", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks,
        ILV_tv = ILV_tv,
        ILVi = c("bool_ILV")
    )

    # id type
    stan_code <- generate_STb_model(data_list,
        veff_params = c("bool_ILV"),
        veff_type = "id"
    )
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )

    # trial type
    stan_code <- generate_STb_model(data_list,
        veff_params = c("bool_ILV"),
        veff_type = "trial"
    )
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )

    # trial type
    stan_code <- generate_STb_model(data_list,
        veff_params = c("bool_ILV"),
        veff_type = c("id", "trial")
    )
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("Boolean ILV OADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks,
        ILV_tv = ILV_tv,
        ILVi = c("bool_ILV")
    )

    stan_code <- generate_STb_model(data_list, data_type = "order")

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("Boolean ILV OADA veffID compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks,
        ILV_tv = ILV_tv,
        ILVi = c("bool_ILV")
    )

    stan_code <- generate_STb_model(data_list,
        data_type = "order",
        veff_params = c("bool_ILV"),
        veff_type = "ID"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
    stan_code <- generate_STb_model(data_list,
        data_type = "order",
        veff_params = c("bool_ILV"),
        veff_type = "trial"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

#### Continuous ILV ####
test_that("Continuous ILV cTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks,
        ILV_tv = ILV_tv,
        ILVi = c("cont_ILV")
    )

    stan_code <- generate_STb_model(data_list)

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("Continuous ILV OADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks,
        ILV_tv = ILV_tv,
        ILVi = c("cont_ILV")
    )

    stan_code <- generate_STb_model(data_list, data_type = "order")

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("Continuous ILV cTADA veffID compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks,
        ILV_tv = ILV_tv,
        ILVi = c("cont_ILV")
    )

    stan_code <- generate_STb_model(data_list,
        veff_params = c("lambda_0", "s", "cont_ILV")
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("Continuous ILV OADA veffID compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks,
        ILV_tv = ILV_tv,
        ILVi = c("cont_ILV")
    )

    stan_code <- generate_STb_model(data_list,
        data_type = "order",
        veff_params = c("cont_ILV")
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})
