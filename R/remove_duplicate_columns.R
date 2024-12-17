#' Utility function to remove duplicated rows
#'
#' @param df a dataframe with possibly duplicated columns
#'
#' @return dataframe
remove_duplicate_columns <- function(df) {
    duplicated_cols <- which(duplicated(as.list(df), fromLast = FALSE))
    df <- df[, -duplicated_cols, drop = FALSE]
    return(df)
}
