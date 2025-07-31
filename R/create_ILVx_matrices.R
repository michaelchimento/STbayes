#' generate_X_matrix()
#'
#' Helper function to create model matrices for categorical ILVs
#'
#' @param ilv_vector c() of ILV values
#' @param ilv_name name to use for column
#' @param n_levels number of levels in factor
#'
#' @returns matrix of size [N,L-1] where L=num levels
#' @export
#'
generate_X_matrix <- function(ilv_vector, ilv_name, n_levels) {
    # if there's only one level, return 0-col matrix
    if (n_levels <= 1) {
        message(sprintf("\u26A0\uFE0F Categorical ILV '%s' has only one level; no contrasts created", ilv_name))
        return(matrix(ncol = 0, nrow = length(ilv_vector)))
    }

    # drop reference level (1)
    mm <- matrix(0, nrow = length(ilv_vector), ncol = n_levels - 1)

    for (lvl in 2:n_levels) {
        mm[, lvl - 1] <- as.integer(ilv_vector == lvl)
    }

    colnames(mm) <- paste0(ilv_name, "_", 2:n_levels)

    return(mm)
}
