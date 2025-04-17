#' Helper function to generate network terms for STAN
#'
#' @param transmission_func
#' @param is_distribution
#' @param separate_s
#' @param veff_ID
#' @param net_var
#' @param net_index
#' @param id_var
#' @param trial_var
#' @param time_var
#' @export
#' @return
#'
#' @examples
get_network_term <- function(transmission_func, is_distribution = FALSE,
                             separate_s = FALSE, num_networks=1, veff_ID = c(), net_var = "A",
                             net_index = "network", id_var = "id", trial_var = "trial", time_var = "time_step") {

  net_effect_term = if(id_var=="j") "net_effect_j" else "net_effect"

  s_prime_term <- if ("s" %in% veff_ID) glue::glue("s_prime[network,{id_var}]") else "s_prime[network]"

    w_term <- if ("w" %in% veff_ID) glue::glue("w[network,{id_var}]") else "w[network]"

    net_expr <- if (is_distribution) {
        glue::glue("{net_var}[{net_index}][{id_var}, ]")
    } else {
        glue::glue("{net_var}[{net_index}, {trial_var}, {time_var}][{id_var}, ]")
    }

    base_term <- glue::glue("sum({net_expr} .* Z[{trial_var}][{time_var}, ])")

    if (transmission_func == "standard") {
        return(glue::glue("real {net_effect_term} = 0;
for (network in 1:N_networks) {{
  {net_effect_term} += {if (separate_s & num_networks>1) {s_prime_term} else if (!separate_s & num_networks>1) {w_term} else '1.0'} * {base_term};
}}"))
    }

    if (transmission_func == "freqdep_f") {
        f_term <- if ("f" %in% veff_ID) glue::glue("f[{id_var}]") else "f"
        return(glue::glue("real {net_effect_term} = 0;
for (network in 1:N_networks) {{
  real active = {base_term};
  real inactive = sum({net_expr} .* (1 - Z[{trial_var}][{time_var}, ]));
  real frac = pow(active, {f_term}) / (pow(active, {f_term}) + pow(inactive, {f_term}));
  {net_effect_term} += {if (separate_s & num_networks>1) {s_prime_term} else if (!separate_s & num_networks>1) {w_term} else '1.0'} * frac;
}}"))
    }

    if (transmission_func == "freqdep_k") {
        k_term <- if ("k" %in% veff_ID) glue::glue("k_shape[{id_var}]") else "k_shape"
        return(glue::glue("real {net_effect_term} = 0;
for (network in 1:N_networks) {{
  real numer = {base_term};
  real denom = numer + sum({net_expr} .* (1 - Zn[{trial_var}][{time_var}, ]));
  real prop = denom > 0 ? numer / denom : 0.0;
  real dini_transformed = dini_func(prop, {k_term});
  {net_effect_term} += {if (separate_s & num_networks>1) {s_prime_term} else if (!separate_s & num_networks>1) {w_term} else '1.0'} * dini_transformed;
}}"))
    }

    stop("Unknown transmission_func")
}
