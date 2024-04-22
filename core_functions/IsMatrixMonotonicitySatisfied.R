#' This function checks if the toxicity matrix is monotone
#'
#' @param toxicity_matrix: Matrix with the probabilities a toxicity
#' @param return TRUE if the toxicity matrix is monotone

IsMatrixMonotonicitySatisfied <- function(toxicity_matrix){
  dim <- nrow(toxicity_matrix)*ncol(toxicity_matrix)
  if(sum(t(apply(toxicity_matrix, 1, sort) ) == toxicity_matrix) != dim){
    return(FALSE)
  }
  if(sum(apply(toxicity_matrix,2, sort, decreasing = T) == toxicity_matrix) != dim){
    return(FALSE)
  }
  return(TRUE)
}
