#' process_ILVs()
#'
#' Helper function that takes ilv_vars and creates appropriately formatted stan code
#'
#' @param ilv_vars names of ilvs with prefixes
#' @param ilv_vars_clean "clean" names of ilvs without prefixes
#' @param ilv_datatypes named vector of data types "boolean", "continuous" or "categorical
#' @param ilv_n_levels named vector of number of levels for categorical data types
#' @param ilv_timevarying named vector of boolean vals for whether or not the ILV is timevarying
#' @param veff_params vector of parameters that need veffs
#' @param veff_type string for declaring vectors, either pop size "P" or trial "K"
#' @param suffix character "i" or "j"
#' @param STb_data user imported data
#' @param count_start integer count variable for indexing
#' @param prior_beta string defining prior eg "normal(0,1)"
#'
#' @return list containing different bits of stan declarations
process_ILVs <- function(ilv_vars, ilv_vars_clean, ilv_datatypes, ilv_n_levels, ilv_timevarying, veff_params, veff_type, suffix,
                         STb_data, count_start, prior_beta) {
    transformed_decl <- c()
    transformed_calc <- c()
    modified_vars <- ilv_vars
    param_lines <- c()
    prior_lines <- c()
    count <- count_start

    if (length(ilv_vars) < 1) {
        return(list(
            param = "",
            prior = "",
            term = if (suffix == "i") "1.0" else "",
            transformed_decl = transformed_decl,
            transformed_calc = transformed_calc,
            count = count,
            modified_vars = modified_vars
        ))
    }

    v_term <- c()
    if (is.element("id", veff_type)) v_term <- append(v_term, "v_id[,{count}]")
    if (is.element("trial", veff_type)) v_term <- append(v_term, "v_trial[trial,{count}]")
    v_term <- paste(v_term, collapse = " + ")

    for (ilv in ilv_vars) {
        datatype <- ilv_datatypes[[paste0("ILV_", ilv)]]
        if (is.null(datatype) || is.na(datatype)) {
            stop(glue::glue("ILV datatype missing for variable '{ilv}'"))
        }
        var_name_tag <- paste0(ilv, "_", suffix)
        is_timevarying <- ilv_timevarying[[paste0("ILV_", ilv)]]

        # parameter declaration
        if (datatype %in% c("categorical", "boolean")) {
            n_levels <- ilv_n_levels[[paste0("ILV_", ilv)]]
            n_par <- n_levels - 1
            param_lines <- c(param_lines, glue::glue("vector[{n_par}] beta_ILV{suffix}_{ilv};"))
        } else if (datatype == "continuous") {
            param_lines <- c(param_lines, glue::glue("real beta_ILV{suffix}_{ilv};"))
        } else {
            stop(glue::glue("unknown ILV datatype for '{ilv}': {datatype}"))
        }

        # set prior declaration
        prior_lines <- c(prior_lines, glue::glue("beta_ILV{suffix}_{ilv} ~ {prior_beta};"))

        # build the v_id term only if needed
        has_veff <- ilv %in% veff_params
        v_term <- if (has_veff) glue::glue(" + {glue::glue(v_term)}") else ""

        # transformed parameter declaration
        if (is_timevarying) { # if timevarying
            transformed_decl <- c(
                transformed_decl,
                glue::glue(
                    "array[K, T_max] vector[P] {var_name_tag};\n"
                )
            )

            transformed_calc <- c(
                transformed_calc,
                glue::glue(
                    "for (trial in 1:K) ",
                    "for (timestep in 1:T_max) ",
                    "{var_name_tag}[trial][timestep] = ",
                    "ILV_{ilv}[trial][timestep] * beta_ILV{suffix}_{ilv}{v_term};"
                )
            )

            # update variable names and increment veff count
            modified_vars[modified_vars == ilv] <- glue::glue("{var_name_tag}[trial,time_step,id]")
            if (has_veff) count <- count + 1
        } else {
            if (is.element("id", veff_type) & !is.element("trial", veff_type)) {
                transformed_decl <- append(
                    transformed_decl,
                    glue::glue(
                        "vector[P] {var_name_tag};\n"
                    )
                )
                transformed_calc <- append(
                    transformed_calc,
                    glue::glue(
                        "{var_name_tag} = ILV_{ilv} * beta_ILV{suffix}_{ilv}{v_term};"
                    )
                )
                modified_vars[modified_vars == ilv] <- glue::glue("{var_name_tag}[id]")
            } else if (is.element("trial", veff_type)) {
                transformed_decl <- append(
                    transformed_decl,
                    glue::glue(
                        "array[K] vector[P] {var_name_tag};\n"
                    )
                )
                transformed_calc <- append(
                    transformed_calc,
                    glue::glue(
                        "for (trial in 1:K) ",
                        "{var_name_tag}[trial] = ILV_{ilv} * beta_ILV{suffix}_{ilv}{v_term};"
                    )
                )
                modified_vars[modified_vars == ilv] <- glue::glue("{var_name_tag}[trial,id]")
            }

            # increment veff count
            if (has_veff) count <- count + 1
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
        transformed_decl = transformed_decl,
        transformed_calc = transformed_calc,
        count = count,
        modified_vars = modified_vars
    ))
}
