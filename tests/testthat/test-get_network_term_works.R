test_that("get_network_term generates correct Stan code across conditions", {
    #### NO DISTRIBUTION, NO VEFF ####
    for (func in c("standard", "freqdep_k", "freqdep_f")) {
        network_term <- get_network_term(
            transmission_func = func,
            is_distribution = FALSE,
            num_networks = 2,
            veff_params = c()
        )
        expected <- switch(func,
            standard = "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime[network] * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
}",
            freqdep_k = "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape[network]);
  net_effect += s_prime[network] * dini_transformed;
}",
            freqdep_f = "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f[network]) / (pow(active, f[network]) + pow(inactive, f[network]));
  }
  net_effect += s_prime[network] * frac;
}"
        )
        testthat::expect_equal(network_term, expected)
    }

    #### IS DISTRIBUTION, NO VEFF ####
    for (func in c("standard", "freqdep_k", "freqdep_f")) {
        network_term <- get_network_term(
            transmission_func = func,
            is_distribution = TRUE,
            num_networks = 2,
            veff_params = c()
        )
        expected <- switch(func,
            standard = "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime[network] * dot_product(A[network][id, ],Z[trial][time_step, ]);
}",
            freqdep_k = "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape[network]);
  net_effect += s_prime[network] * dini_transformed;
}",
            freqdep_f = "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f[network]) / (pow(active, f[network]) + pow(inactive, f[network]));
  }
  net_effect += s_prime[network] * frac;
}"
        )
        testthat::expect_equal(network_term, expected)
    }

    #### NO DISTRIBUTION, VEFF (s) ####
    for (func in c("standard", "freqdep_k", "freqdep_f")) {
        network_term <- get_network_term(
            transmission_func = func,
            is_distribution = FALSE,
            num_networks = 2,
            veff_params = c("s")
        )
        expected <- switch(func,
            standard = "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime[network,id] * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
}",
            freqdep_k = "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape[network]);
  net_effect += s_prime[network,id] * dini_transformed;
}",
            freqdep_f = "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f[network]) / (pow(active, f[network]) + pow(inactive, f[network]));
  }
  net_effect += s_prime[network,id] * frac;
}"
        )
        testthat::expect_equal(network_term, expected)
    }

    #### IS DISTRIBUTION, VEFF (s) ####
    for (func in c("standard", "freqdep_k", "freqdep_f")) {
        network_term <- get_network_term(
            transmission_func = func,
            is_distribution = TRUE,
            num_networks = 2,
            veff_params = c("s")
        )
        expected <- switch(func,
            standard = "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime[network,id] * dot_product(A[network][id, ],Z[trial][time_step, ]);
}",
            freqdep_k = "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape[network]);
  net_effect += s_prime[network,id] * dini_transformed;
}",
            freqdep_f = "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f[network]) / (pow(active, f[network]) + pow(inactive, f[network]));
  }
  net_effect += s_prime[network,id] * frac;
}"
        )
        testthat::expect_equal(network_term, expected)
    }

    #### SINGLE NETWORK CASE ####
    for (func in c("standard", "freqdep_k", "freqdep_f")) {
        network_term <- get_network_term(
            transmission_func = func,
            is_distribution = FALSE,
            num_networks = 1,
            veff_params = c()
        )
        expected <- switch(func,
            standard = "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
}",
            freqdep_k = "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape);
  net_effect += s_prime * dini_transformed;
}",
            freqdep_f = "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f) / (pow(active, f) + pow(inactive, f));
  }
  net_effect += s_prime * frac;
}"
        )
        testthat::expect_equal(network_term, expected)
    }

    #### SINGLE NETWORK + VEFF F K CASE ####
    for (func in c("standard", "freqdep_k", "freqdep_f")) {
        network_term <- get_network_term(
            transmission_func = func,
            is_distribution = FALSE,
            num_networks = 1,
            veff_params = c("f", "k")
        )
        expected <- switch(func,
            standard = "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
}",
            freqdep_k = "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape[id]);
  net_effect += s_prime * dini_transformed;
}",
            freqdep_f = "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f[id]) / (pow(active, f[id]) + pow(inactive, f[id]));
  }
  net_effect += s_prime * frac;
}"
        )
        testthat::expect_equal(network_term, expected)
    }
})
