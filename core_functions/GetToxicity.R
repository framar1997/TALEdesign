#' This function extracts the probability of toxicity of DL (i,j) from the toxicity matrix
#'
#' @param i: level of the first drug
#' @param j: level of the second drug drug
#' @param toxicity_matrix: matrix of toxicities
#' @param return the probability of toxicity at DL (i,j)
#'

GetToxicity = function(i, j, toxicity_matrix)
{
  x = which(colnames(toxicity_matrix) == toString(i) )
  y = which(rownames(toxicity_matrix) == toString(j) )
  return(toxicity_matrix[y,x])
}
