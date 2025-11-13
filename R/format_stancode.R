#' format_stancode
#'
#' @param stan_code string of stan code, duh
#' @param indent integer of num spaces for indent
#'
#' @return formatted string of stan code
format_stancode <- function(stan_code, indent = 4) {
    lines <- unlist(strsplit(stan_code, "\n"))
    level <- 0
    indented <- character(length(lines))

    for (i in seq_along(lines)) {
        line <- trimws(lines[i])
        # decrease indent if line starts with }
        if (grepl("^\\}", line)) level <- max(0, level - 1)
        indented[i] <- paste0(strrep(" ", indent * level), line)
        # increase indent if line contains {
        if (grepl("\\{\\s*$", line)) level <- level + 1
    }
    paste(indented, collapse = "\n")
}
