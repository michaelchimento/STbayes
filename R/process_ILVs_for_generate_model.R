#' process_ILVs()
#'
#' Helper function that takes ilv_vars and creates appropriately formatted stan code
#'
#' @param ilv_vars names of ilvs with prefixes
#' @param ilv_vars_clean "clean" names of ilvs without prefixes
#' @param ilv_datatypes named vector of data types "boolean", "continuous" or "categorical
#' @param ilv_n_levels named vector of number of levels for categorical data types
#' @param ilv_timevarying named vector of boolean vals for whether or not the ILV is timevarying
#' @param veff_ID vector of parameters that need veffs
#' @param suffix character "i" or "j"
#' @param STb_data user imported data
#' @param count_start integer count variable for indexing
#' @param prior_beta string defining prior eg "normal(0,1)"
#'
#' @return list containing different bits of stan declarations
process_ILVs <- function(ilv_vars, ilv_vars_clean, ilv_datatypes, ilv_n_levels, ilv_timevarying, veff_ID, suffix,
                         STb_data, count_start, prior_beta) {
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
        datatype <- ilv_datatypes[[paste0("ILV_", ilv)]]
        if (is.null(datatype) || is.na(datatype)) {
            stop(glue::glue("ILV datatype missing for variable '{ilv}'"))
        }
        veff_tag <- paste0(ilv, "_", suffix)
        is_timevarying <- ilv_timevarying[[paste0("ILV_", ilv)]]

        # parameter declaration
        if (datatype %in% c("categorical", "boolean")) {
            n_levels <- ilv_n_levels[[paste0("ILV_", ilv)]]
            n_par <- n_levels - 1
            param_lines <- c(param_lines, glue::glue("vector[{n_par}] beta_ILV{suffix}_{ilv};"))
            prior_lines <- c(prior_lines, glue::glue("beta_ILV{suffix}_{ilv} ~ {prior_beta};"))
        } else if (datatype == "continuous") {
            param_lines <- c(param_lines, glue::glue("real beta_ILV{suffix}_{ilv};"))
            prior_lines <- c(prior_lines, glue::glue("beta_ILV{suffix}_{ilv} ~ {prior_beta};"))
        } else {
            stop(glue::glue("unknown ILV datatype for '{ilv}': {datatype}"))
        }

        # transformed parameter declaration
        if (is_timevarying) { # if timevarying
            if (ilv %in% veff_ID) { # if variable is in veff_ID
                if (datatype != "continuous") {
                    transformed <- c(transformed, glue::glue(
                        "array[K, T_max] vector[P] {veff_tag};\nfor (trial in 1:K) for (timestep in 1:T_max) {veff_tag}[trial][timestep] = ILV_{ilv}[trial][timestep] * beta_ILV{suffix}_{ilv} + v_ID[,{count}];"
                    ))
                    modified_vars[modified_vars == ilv] <- glue::glue("{veff_tag}[trial,time_step,id]")
                    count <- count + 1
                } else {
                    transformed <- c(transformed, glue::glue(
                        "array[K, T_max] vector[P] {veff_tag};\nfor (trial in 1:K) for (timestep in 1:T_max) {veff_tag}[trial][timestep] = ILV_{ilv}[trial][timestep] * beta_ILV{suffix}_{ilv} + v_ID[,{count}];"
                    ))
                    modified_vars[modified_vars == ilv] <- glue::glue("{veff_tag}[trial,time_step,id]")
                    count <- count + 1
                }
            } else { # if variable is NOT veff_ID
                if (datatype != "continuous") {
                    transformed <- c(transformed, glue::glue(
                        "array[K, T_max] vector[P] {veff_tag};\nfor (trial in 1:K) for (timestep in 1:T_max) {veff_tag}[trial][timestep] = ILV_{ilv}[trial][timestep] * beta_ILV{suffix}_{ilv};"
                    ))
                    modified_vars[modified_vars == ilv] <- glue::glue("{veff_tag}[trial,time_step,id]")
                } else {
                    transformed <- c(transformed, glue::glue(
                        "array[K, T_max] vector[P] {veff_tag};\nfor (trial in 1:K) for (timestep in 1:T_max) {veff_tag}[trial][timestep] = ILV_{ilv}[trial][timestep] * beta_ILV{suffix}_{ilv};"
                    ))
                    modified_vars[modified_vars == ilv] <- glue::glue("{veff_tag}[trial,time_step,id]")
                }
            }
        } else {
            if (datatype != "continuous") {
                if (ilv %in% veff_ID) {
                    transformed <- append(transformed, glue::glue(
                        "vector[P] {veff_tag} = ILV_{ilv} * beta_ILV{suffix}_{ilv} + v_ID[,{count}];"
                    ))
                    modified_vars[modified_vars == ilv] <- glue::glue("{veff_tag}[id]")
                    count <- count + 1
                } else {
                    transformed <- append(transformed, glue::glue(
                        "vector[P] {veff_tag} = ILV_{ilv} * beta_ILV{suffix}_{ilv};"
                    ))
                    modified_vars[modified_vars == ilv] <- glue::glue("{veff_tag}[id]")
                }
            } else {
                if (ilv %in% veff_ID) {
                    transformed <- append(transformed, glue::glue(
                        "vector[P] {veff_tag} = ILV_{ilv} * beta_ILV{suffix}_{ilv} + v_ID[,{count}];"
                    ))
                    modified_vars[modified_vars == ilv] <- glue::glue("{veff_tag}[id]")
                    count <- count + 1
                } else {
                    transformed <- append(transformed, glue::glue(
                        "vector[P] {veff_tag} = ILV_{ilv} * beta_ILV{suffix}_{ilv};"
                    ))
                    modified_vars[modified_vars == ilv] <- glue::glue("{veff_tag}[id]")
                }
            }
        }
    }

    # construct multiplicative term
    term <- paste0(
        if (suffix == "m") "exp(" else if (suffix == "s") "* exp(" else "exp(",
        paste(modified_vars, collapse = " + "),
        ")",
        if (suffix == "m") " *" else ""
    )

    return(list(
        param = paste(param_lines, collapse = "\n"),
        prior = paste(prior_lines, collapse = "\n"),
        term = term,
        transformed = transformed,
        count = count,
        modified_vars = modified_vars
    ))
}
