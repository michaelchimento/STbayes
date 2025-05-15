test_that("get_network_term generates network terms correctly", {
  #### NO DISTRIBUTION ####
  # standard one
  network_term <- STbayes::get_network_term(
    transmission_func = "standard",
    is_distribution = FALSE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime[network] * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
}"
  )

  # complex_k one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_k",
    is_distribution = FALSE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape);
  net_effect += s_prime[network] * dini_transformed;
}"
  )

  # complex_f one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_f",
    is_distribution = FALSE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f) / (pow(active, f) + pow(inactive, f));
  }
  net_effect += s_prime[network] * frac;
}"
  )
  #### IS DISTRIBUTION ####
  # standard one
  network_term <- STbayes::get_network_term(
    transmission_func = "standard",
    is_distribution = TRUE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime[network] * dot_product(A[network][id, ],Z[trial][time_step, ]);
}"
  )

  # complex_k one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_k",
    is_distribution = TRUE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape);
  net_effect += s_prime[network] * dini_transformed;
}"
  )

  # complex_f one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_f",
    is_distribution = TRUE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f) / (pow(active, f) + pow(inactive, f));
  }
  net_effect += s_prime[network] * frac;
}"
  )

  #### NO DISTRIBUTION + VEFF ####
  # standard one
  network_term <- STbayes::get_network_term(
    transmission_func = "standard",
    is_distribution = FALSE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c("s")
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime[network,id] * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
}"
  )

  # complex_k one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_k",
    is_distribution = FALSE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c("s")
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape);
  net_effect += s_prime[network,id] * dini_transformed;
}"
  )

  # complex_f one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_f",
    is_distribution = FALSE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c("s")
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f) / (pow(active, f) + pow(inactive, f));
  }
  net_effect += s_prime[network,id] * frac;
}"
  )
  #### IS DISTRIBUTION + VEFF ####
  # standard one
  network_term <- STbayes::get_network_term(
    transmission_func = "standard",
    is_distribution = TRUE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c("s")
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime[network,id] * dot_product(A[network][id, ],Z[trial][time_step, ]);
}"
  )

  # complex_k one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_k",
    is_distribution = TRUE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c("s")
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape);
  net_effect += s_prime[network,id] * dini_transformed;
}"
  )

  # complex_f one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_f",
    is_distribution = TRUE,
    separate_s = TRUE,
    num_networks = 2,
    veff_ID = c("s")
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f) / (pow(active, f) + pow(inactive, f));
  }
  net_effect += s_prime[network,id] * frac;
}"
  )

  #### NO DISTRIBUTION + w[n] ####
  # standard one
  network_term <- STbayes::get_network_term(
    transmission_func = "standard",
    is_distribution = FALSE,
    separate_s = FALSE,
    num_networks = 1,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
}"
  )

  # complex_k one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_k",
    is_distribution = FALSE,
    separate_s = FALSE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape);
  net_effect += s_prime * dini_transformed;
}"
  )

  # complex_f one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_f",
    is_distribution = FALSE,
    separate_s = FALSE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
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
  #### IS DISTRIBUTION ####
  # standard one
  network_term <- STbayes::get_network_term(
    transmission_func = "standard",
    is_distribution = TRUE,
    separate_s = FALSE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime * dot_product(A[network][id, ],Z[trial][time_step, ]);
}"
  )

  # complex_k one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_k",
    is_distribution = TRUE,
    separate_s = FALSE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape);
  net_effect += s_prime * dini_transformed;
}"
  )

  # complex_f one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_f",
    is_distribution = TRUE,
    separate_s = FALSE,
    num_networks = 2,
    veff_ID = c()
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f) / (pow(active, f) + pow(inactive, f));
  }
  net_effect += s_prime * frac;
}"
  )

  #### NO DISTRIBUTION + VEFF ####
  # standard one
  network_term <- STbayes::get_network_term(
    transmission_func = "standard",
    is_distribution = FALSE,
    separate_s = FALSE,
    num_networks = 2
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime * dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
}"
  )

  # complex_k one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_k",
    is_distribution = FALSE,
    separate_s = FALSE,
    num_networks = 2,
    veff_ID = c("k")
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network, trial, time_step][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network, trial, time_step][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape[id]);
  net_effect += s_prime * dini_transformed;
}"
  )

  # complex_f one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_f",
    is_distribution = FALSE,
    separate_s = FALSE,
    num_networks = 2,
    veff_ID = c("f")
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
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
  #### IS DISTRIBUTION + VEFF ####
  # standard one
  network_term <- STbayes::get_network_term(
    transmission_func = "standard",
    is_distribution = TRUE,
    separate_s = FALSE,
    num_networks = 2
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  net_effect += s_prime * dot_product(A[network][id, ],Z[trial][time_step, ]);
}"
  )

  # complex_k one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_k",
    is_distribution = TRUE,
    separate_s = FALSE,
    num_networks = 2,
    veff_ID = c("k")
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real numer = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real denom = numer + dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, k_shape[id]);
  net_effect += s_prime * dini_transformed;
}"
  )

  # complex_f one
  network_term <- STbayes::get_network_term(
    transmission_func = "freqdep_f",
    is_distribution = TRUE,
    separate_s = FALSE,
    num_networks = 2,
    veff_ID = c("f")
  )
  testthat::expect_equal(
    network_term,
    "real net_effect = 0;
for (network in 1:N_networks) {
  real active = dot_product(A[network][id, ],Z[trial][time_step, ]);
  real inactive = dot_product(A[network][id, ], (1 - Zn[trial][time_step, ]));
  real frac = 0;
  if ((active + inactive)>0){
    frac = pow(active, f[id]) / (pow(active, f[id]) + pow(inactive, f[id]));
  }
  net_effect += s_prime * frac;
}"
  )
})

