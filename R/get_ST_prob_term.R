#' Generate Stan code to compute P(e via ST) for %ST calculation
#'
#' @param transmission_func String (e.g. "standard", "freqdep_f", "freqdep_k")
#' @param is_distribution Boolean: Are edges drawn from posterior?
#' @param separate_s Boolean: Is s estimated separately per network?
#' @param veff_ID Character vector: Which params vary by ID
#' @param num_networks Integer: Number of networks
#' @param id_var String: e.g., "id"
#' @param trial_var String: e.g., "trial"
#' @param time_var String: e.g., "time_step"
#' @param net_var String: e.g., "A"
#' @return String of Stan code for GQ block to accumulate psoc and psocn
get_ST_prob_term <- function(transmission_func, is_distribution = FALSE,
                             separate_s = FALSE, veff_ID = c(),
                             num_networks = 1,
                             s_var = "s_prime",
                             id_var = "id", trial_var = "trial", time_var = "time_step",
                             net_var = "A",
                             ILVs_variable_effects = "",
                             weibull_term = "",
                             high_res=F) {
  # choose s term
  s_term <- if (num_networks > 1 && separate_s) {
    if ("s" %in% veff_ID) glue::glue("{s_var}[network, {id_var}]") else glue::glue("{s_var}[network]")
  } else {
    if ("s" %in% veff_ID) glue::glue("{s_var}[{id_var}]") else glue::glue("{s_var}")
  }

  full_s_term <- glue::glue("{s_term} * D[{trial_var}, {time_var}] {ILVs_variable_effects} {weibull_term}")

  net_expr <- if (is_distribution) {
    glue::glue("{net_var}[network][{id_var}, ]")
  } else {
    glue::glue("{net_var}[network, {trial_var}, {time_var}][{id_var}, ]")
  }

  active_expr <- glue::glue("dot_product({net_expr}, Z[{trial_var}][{time_var}, ])")
  inactive_expr <- glue::glue("dot_product({net_expr}, (1 - Zn[{trial_var}][{time_var}, ]))")

  st_lines <- switch(transmission_func,
    "standard" = glue::glue("
    for (network in 1:N_networks) {{
        real Tn = {active_expr};
        psocn_sum[network] += ({full_s_term} * Tn) / lambda;
    }}
    count_ST += 1;
    "),
    "freqdep_f" = {
      f_term <- if ("f" %in% veff_ID) glue::glue("f[{id_var}]") else "f"
      glue::glue("
      for (network in 1:N_networks) {{
          real active = {active_expr};
          real inactive = {inactive_expr};
          real frac = pow(active, {f_term}) / (pow(active, {f_term}) + pow(inactive, {f_term}));
          psocn_sum[network] += ({full_s_term} * frac) / lambda;
      }}
      count_ST += 1;
      ")
    },
    "freqdep_k" = {
      k_term <- if ("k" %in% veff_ID) glue::glue("k_shape[{id_var}]") else "k_shape"
      if (!high_res){
        glue::glue("
      for (network in 1:N_networks) {{
          real numer = {active_expr};
          real denom = numer + dot_product({net_expr}, (1 - Zn[{trial_var}][{time_var}, ]));
          real prop = denom > 0 ? numer / denom : 0.0;
          real dini = dini_func(prop, {k_term});
          psocn_sum[network] += ({full_s_term} * dini) / lambda;
      }}
      count_ST += 1;
      ")
      }
    },
    stop("Unsupported transmission_func")
  )

  return(st_lines)
}
