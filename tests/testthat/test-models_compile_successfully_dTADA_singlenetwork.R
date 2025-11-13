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


#### BOOLEAN ILV ####
test_that("No veff dTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        data_type = "discrete_time"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("lambda_0 veff_type=id dTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        veff_params = "lambda_0",
        veff_type = "id",
        data_type = "discrete_time"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("lambda_0 veff_type=trial dTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        veff_params = "lambda_0",
        veff_type = "trial",
        data_type = "discrete_time"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("lambda_0 veff_type=c(id,trial) dTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        veff_params = "lambda_0",
        veff_type = c("id", "trial"),
        data_type = "discrete_time"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("veff_params=c(lambda_0,s) veff_type=c(id,trial) dTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        veff_params = c("lambda_0", "s"),
        veff_type = c("id", "trial"),
        data_type = "discrete_time"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("veff_params=c(lambda_0,s,gamma) veff_type=c(id,trial) dTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        intrinsic_rate = "weibull",
        veff_params = c("lambda_0", "s", "gamma"),
        veff_type = c("id", "trial"),
        data_type = "discrete_time"
    )
    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("veff_params=c(lambda_0,s,gamma, f) veff_type=c(id,trial) dTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        intrinsic_rate = "weibull",
        transmission_func = "freqdep_f",
        veff_params = c("lambda_0", "s", "gamma", "f"),
        veff_type = c("id", "trial"),
        data_type = "discrete_time"
    )
    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("veff_params=c(lambda_0,s,gamma,k) veff_type=c(id,trial) dTADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        intrinsic_rate = "weibull",
        transmission_func = "freqdep_k",
        veff_params = c("lambda_0", "s", "gamma", "k"),
        veff_type = c("id", "trial"),
        data_type = "discrete_time"
    )
    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})
