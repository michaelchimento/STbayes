#' Helper function to generate network terms for STAN
#'
#' @param transmission_func simple or complex transmission
#' @param is_distribution did the user supply a posterior distribution for edges
#' @param separate_s whether or not s should be constrained to be the same for all networks
#' @param veff_ID vector of varying effects
#' @param s_var string used to represent s term ("s_prime" or "s_direct")
#' @param net_var string used to represent network variable "A"
#' @param net_index string used to index networks "network"
#' @param id_var string used to index individuals. "id" for tada, need to use "id" and "j" for oada
#' @param trial_var string used to index trial "trial"
#' @param time_var string used to index time "time_step"
#' @export
#' @return string of stan code to be used in the model for calculating network effects
get_network_term <- function(transmission_func="standard", is_distribution = FALSE,
                             separate_s = FALSE, num_networks = 1, veff_ID = c(), s_var="s_prime", net_var = "A",
                             net_index = "network", id_var = "id", trial_var = "trial", time_var = "time_step", high_res=F) {
  net_effect_term <- if (id_var == "j") "net_effect_j" else "net_effect"

  s_term <- if (num_networks > 1 & separate_s) {
    if ("s" %in% veff_ID) glue::glue("{s_var}[{net_index},{id_var}]") else glue::glue("{s_var}[{net_index}]")
  } else {
    if ("s" %in% veff_ID) glue::glue("{s_var}[{id_var}]") else glue::glue("{s_var}")
  }

  # w_term <- if ("w" %in% veff_ID) glue::glue("w[network,{id_var}]") else "w[network]"

  net_expr <- if (is_distribution) {
    glue::glue("{net_var}[{net_index}][{id_var}, ]")
  } else {
    glue::glue("{net_var}[{net_index}, {trial_var}, {time_var}][{id_var}, ]")
  }

  base_term <- glue::glue("dot_product({net_expr},Z[{trial_var}][{time_var}, ])")

  if (transmission_func == "standard") {
    return(glue::glue("real {net_effect_term} = 0;
for (network in 1:N_networks) {{
  {net_effect_term} += {s_term} * {base_term};
}}"))
  }

  if (transmission_func == "freqdep_f") {
    f_term <- if ("f" %in% veff_ID) glue::glue("f[{id_var}]") else "f"
    return(glue::glue("real {net_effect_term} = 0;
for (network in 1:N_networks) {{
  real active = {base_term};
  real inactive = dot_product({net_expr}, (1 - Zn[{trial_var}][{time_var}, ]));
  real frac = 0;
  if ((active + inactive)>0){{
    frac = pow(active, {f_term}) / (pow(active, {f_term}) + pow(inactive, {f_term}));
  }}
  {net_effect_term} += {s_term} * frac;
}}"))
  }

  if (transmission_func == "freqdep_k" & !high_res) {
    k_term <- if ("k" %in% veff_ID) glue::glue("k_shape[{id_var}]") else "k_shape"
    return(glue::glue("real {net_effect_term} = 0;
for (network in 1:N_networks) {{
  real numer = {base_term};
  real denom = numer + dot_product({net_expr}, (1 - Zn[{trial_var}][{time_var}, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, {k_term});
  {net_effect_term} += {s_term} * dini_transformed;
}}"))
  }

  stop("Unknown transmission_func")
}
