#' Check required cols
#'
#' @param df dataframe
#' @param required_cols vector of col names
#' @param df_name
#'
#' @returns
#' @export
#'
#' @examples
check_required_cols <- function(df, required_cols, df_name = "dataframe") {
    missing <- setdiff(required_cols, names(df))
    if (length(missing) > 0) {
        stop(glue::glue(
            "❌ {df_name} is missing required column(s): {paste(missing, collapse = ', ')}"
        ), call. = FALSE)
    }
    invisible(TRUE)
}
