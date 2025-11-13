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
    weight = 1,
    kin = 1
)


#### BOOLEAN ILV ####
test_that("No veff OADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        data_type = "order"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("s veff_type=id OADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        veff_params = "s",
        veff_type = "id",
        data_type = "order"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("s veff_type=trial OADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        veff_params = "s",
        veff_type = "trial",
        data_type = "order"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("s veff_type=c(id,trial) OADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        veff_params = "s",
        veff_type = c("id", "trial"),
        data_type = "order"
    )

    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})



test_that("veff_params=c(s, f) veff_type=c(id,trial) OADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        intrinsic_rate = "weibull",
        transmission_func = "freqdep_f",
        veff_params = c("s", "f"),
        veff_type = c("id", "trial"),
        data_type = "order"
    )
    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})

test_that("veff_params=c(s,k) veff_type=c(id,trial) OADA compiles", {
    data_list <- import_user_STb(
        event_data = event_data,
        networks = networks
    )

    stan_code <- generate_STb_model(data_list,
        transmission_func = "freqdep_k",
        veff_params = c("s", "k"),
        veff_type = c("id", "trial"),
        data_type = "order"
    )
    # write to a temporary .stan file
    tmp_stan_path <- cmdstanr::write_stan_file(stan_code)

    # try to compile
    expect_error(
        cmdstanr::cmdstan_model(tmp_stan_path, compile = TRUE),
        regexp = NA # expect no error
    )
})
