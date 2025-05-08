#' extract_ILV: private function to extract the ILV matrix and back transform into vectors for STb format
#'
#' @param nbda_object NBDAdata object
#' @param ILV_type character indicating types of ILVs: @asocILVdata = asocial, @intILVdata = social, and @multiILVdata = multiplicative (affects both social and asocial equally)
#'
#' @return List of vectors for each ILV
extract_ILV <- function(nbda_object, ILV_type) {
  ILV_matrix <- methods::slot(nbda_object, ILV_type) # Access the slot via its name

  if (sum(ILV_matrix) == 0) message("Warning, ILV values are all 0, is this correct?")
  # transform each col of the ILV matrix into a vector length N
  ILV_vectors <- lapply(1:ncol(ILV_matrix), function(row_idx) {
    # match ILVs to IDs using the order of acquisition
    ILV_vector <- ILV_matrix[order(nbda_object@idname), row_idx]
    names(ILV_vector) <- nbda_object@idname # Add ID names for clarity
    return(ILV_vector)
  })
  # combine ILV vectors into a data frame
  ILV_df <- as.data.frame(ILV_vectors)
  colnames(ILV_df) <- colnames(ILV_matrix)

  return(ILV_df)
}
