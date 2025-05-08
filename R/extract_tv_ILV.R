#' extract_tv_ILV: private function to extract the tv ILV matrices and transform them into ILV_tv dataframe for further processing
#' This is a really silly way of doing it:
#' 1. the original TV matrix is impossible to reconstruct from NBDA objs, since NBDA only saves values for when individuals are naive.
#' This shouldn't matter for basic functionality, but if you have a use-case where you wanted to do counterfactuals etc, one needs that data.
#' 2. this basically converts vector > matrix > vector, which will then be converted again into a matrix by the import_NBDA_STb function..
#' However it works, will come back later to make it more efficient.
#'
#' @param nbda_object NBDAdata object
#' @param ILV_type character indicating types of ILVs: @asocILVdata = asocial, @intILVdata = social, and @multiILVdata = multiplicative (affects both social and asocial equally)
#'
#' @return List of vectors for each ILV
extract_tv_ILV <- function(nbda_object, ILV_type) {
  long_format_list <- list()
  ILV_matrix <- methods::slot(nbda_object, ILV_type) # Access the slot via its name
  if (sum(ILV_matrix) == 0) message("Warning, ILV values are all 0, is this correct?")
  nID <- length(methods::slot(nbda_object, "idname"))

  for (col in 1:ncol(ILV_matrix)) {
    flattened_data <- ILV_matrix[, col]
    # Get acquisition events and naive individuals
    orderAcq <- methods::slot(nbda_object, "orderAcq")
    nAcq <- length(orderAcq) # number of acquisition events

    # initialize an empty matrix for reconstruction
    reconstructed_matrix <- matrix(0, nrow = nID, ncol = nAcq)

    # track row index in flattened data
    row_index <- 1

    # iterate over acquisition events
    for (event in 1:nAcq) {
      nonlearners <- which(nbda_object@statusMatrix[, event] == 0) # naive individuals before event
      num_nonlearners <- length(nonlearners)

      if (num_nonlearners > 0) {
        # fill corresponding rows in reconstructed matrix
        reconstructed_matrix[nonlearners, event] <- flattened_data[row_index:(row_index + num_nonlearners - 1)]

        # Move to the next set of rows in flattened_data
        row_index <- row_index + num_nonlearners
      }
    }

    values <- as.vector(reconstructed_matrix) # Flatten matrix to vector

    # Create data frame
    long_data <- data.frame(value = values)
    long_data <- stats::setNames(data.frame(values), colnames(ILV_matrix)[col])

    # Store in list
    long_format_list[[col]] <- long_data
  }

  # Combine all ILVs into a single data frame
  ILV_tv <- do.call(cbind, long_format_list)
  # ILV_tv$id <- rep(1:nID, times = nAcq)  # Repeat each id for all time points
  # ILV_tv$time <- rep(1:nAcq, each = nID)  # Repeat each time across all ids

  return(ILV_tv)
}
