#' optimize_transformed_params: helper function to get rid of duplicate loops
#'
#' @param transformed_params vector of strings
#'
#' @returns string
optimize_transformed_params <- function(transformed_params) {
    # join to single string
    code <- paste(transformed_params, collapse = "\n")

    # regex to match for loops
    pattern <- "(for \\([^)]*\\)\\s*)([^;]+;)"
    matches <- gregexpr(pattern, code, perl = TRUE)
    loops <- regmatches(code, matches)[[1]]

    if (length(loops) == 0) {
        return(code)
    }

    # extract loop header and statement
    df <- data.frame(
        full = loops,
        header = sub(pattern, "\\1", loops, perl = TRUE),
        stmt = sub(pattern, "\\2", loops, perl = TRUE),
        stringsAsFactors = FALSE
    )

    # group identical headers
    grouped <- lapply(split(df, df$header), function(x) {
        header <- unique(x$header)
        body <- paste(x$stmt, collapse = "\n  ")
        paste0(header, "{\n  ", body, "\n}")
    })

    # replace all grouped loops back into code
    # remove originals, then append grouped loops
    code_no_loops <- gsub(pattern, "", code, perl = TRUE)
    new_code <- paste(code_no_loops, paste(unlist(grouped), collapse = "\n"), sep = "\n")

    # clean blank lines
    new_code <- gsub("\n{2,}", "\n", new_code)
    trimws(new_code)
}
