#' get_network_term()
#'
#' Helper function to generate network terms for STAN
#'
#' @param transmission_func simple or complex transmission
#' @param is_distribution did the user supply a posterior distribution for edges
#' @param num_networks integer of number of networks
#' @param veff_params vector of varying effects
#' @param s_var string used to represent s term ("s_prime" or "s_direct")
#' @param net_var string used to represent network variable "A"
#' @param net_index string used to index networks "network"
#' @param id_var string used to index individuals. "id" for tada, need to use "id" and "j" for oada
#' @param veff_type string or character vector of whether veffs are id, trial or both
#' @param veff_idx string used to index parameters for varying effects
#' @param trial_var string used to index trial "trial"
#' @param time_var string used to index time "time_step"
#' @param high_res boolean indicating if high res
#' @return string of stan code to be used in the model for calculating network effects

get_network_term <- function(transmission_func = "standard", is_distribution = FALSE, num_networks = 1, veff_params = c(), s_var = "s_prime", net_var = "A",
                             net_index = "network", id_var = "id", veff_type = c(), veff_idx = "id", trial_var = "trial", time_var = "time_step", high_res = F) {
    net_effect_term <- if (id_var == "j") "net_effect_j" else "net_effect"
    if (id_var == "j" & is.element("id", veff_type)) veff_idx <- gsub("id", "j", veff_idx)

    s_term <- if (num_networks > 1) {
        if ("s" %in% veff_params) glue::glue("{s_var}[{net_index},{veff_idx}]") else glue::glue("{s_var}[{net_index}]")
    } else {
        if ("s" %in% veff_params) glue::glue("{s_var}[{veff_idx}]") else glue::glue("{s_var}")
    }

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
        f_term <- if ("f" %in% veff_params & num_networks == 1) {
            glue::glue("f[{veff_idx}]")
        } else if ("f" %in% veff_params & num_networks > 1) {
            glue::glue("f[{net_index},{veff_idx}]")
        } else if (!is.element("f", veff_params) & num_networks > 1) {
            glue::glue("f[{net_index}]")
        } else {
            "f"
        }
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
        k_term <- if ("k" %in% veff_params & num_networks == 1) {
            glue::glue("k_shape[{veff_idx}]")
        } else if ("k" %in% veff_params & num_networks > 1) {
            glue::glue("k_shape[{net_index},{veff_idx}]")
        } else if (!is.element("k", veff_params) & num_networks > 1) {
            glue::glue("k_shape[{net_index}]")
        } else {
            "k_shape"
        }
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
