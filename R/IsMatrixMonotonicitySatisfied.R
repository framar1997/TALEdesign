# check if the toxicity matrix is monotone. If yes, return TRUE

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
