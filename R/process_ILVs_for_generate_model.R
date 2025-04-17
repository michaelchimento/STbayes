#' Helper function for generating ctada
#'
#' @param ilv_vars
#' @param ilv_vars_clean
#' @param veff_ID
#' @param suffix
#' @param STb_data
#' @param count_start
#'
#' @return
#' @export
#'
#' @examples
process_ILVs <- function(ilv_vars, ilv_vars_clean, veff_ID, suffix, STb_data, count_start) {
    transformed <- c()
    modified_vars <- ilv_vars
    param_lines <- c()
    prior_lines <- c()
    count <- count_start

    if (length(ilv_vars) < 1) {
        return(list(
            param = "",
            prior = "",
            term = if (suffix == "i") "1.0" else "",
            transformed = transformed,
            count = count,
            modified_vars = modified_vars
        ))
    }

    for (ilv in ilv_vars) {
        veff_tag <- paste0(ilv, "_", suffix)
        if (ilv %in% veff_ID) {
            transformed <- append(transformed, glue::glue("vector[P] {veff_tag} = beta_ILV{suffix}_{ilv} + v_ID[,{count}];"))
            modified_vars[modified_vars == ilv] <- glue::glue("{veff_tag}[id]")
            count <- count + 1
        } else {
            modified_vars[modified_vars == ilv] <- glue::glue("beta_ILV{suffix}_{ilv}")
        }
    }

    param_lines <- paste0("real beta_ILV", suffix, "_", ilv_vars_clean, ";")
    prior_lines <- paste0("beta_ILV", suffix, "_", ilv_vars_clean, " ~ normal(0, 1);")

    term <- paste0(
        if (suffix == "m") "exp(" else if (suffix == "s") "* exp(" else "exp(",
        paste0(
            modified_vars,
            " * ",
            sapply(ilv_vars_clean, function(var_clean) {
                var_expr <- paste0("ILV_", var_clean)
                if (!is.null(dim(STb_data[[var_expr]]))) {
                    paste0(var_expr, "[trial,time_step,id]")
                } else {
                    paste0(var_expr, "[id]")
                }
            }),
            collapse = " + "
        ),
        ")",
        if (suffix == "m") " *" else ""
    )

    return(list(
        param = paste0(param_lines, collapse = "\n"),
        prior = paste0(prior_lines, collapse = "\n"),
        term = term,
        transformed = transformed,
        count = count,
        modified_vars = modified_vars
    ))
}
