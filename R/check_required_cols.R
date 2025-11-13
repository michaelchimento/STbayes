#' check_required_cols()
#'
#' Helper function that makes sure user has put in all of the reqd cols
#'
#' @param df dataframe
#' @param required_cols vector of col names
#' @param df_name string for name of df
#'
#' @return nothing, just errors if something wrong
check_required_cols <- function(df, required_cols, df_name = "dataframe") {
    missing <- setdiff(required_cols, names(df))
    if (length(missing) > 0) {
        stop(glue::glue(
            "\u274C {df_name} is missing required column(s): {paste(missing, collapse = ', ')}"
        ), call. = FALSE)
    }
    invisible(TRUE)
}
