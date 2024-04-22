#' This function updates the DLT rate of a DL and the related reccomendation
#'
#' @param row_explore: data frame containing updated information on the latest explored DLs.
#' @param lambda_e: escalation parameter
#' @param lambda_d: de-escalation parameter
#' @param lambda_r: de-escalate and remove parameter
#' @param return recommendations for the latest explored DLs


GetReccomendation <- function(row_explore, lambda_e, lambda_d, lambda_r){
  total_enrolled <- row_explore[8]
  total_DLT <- row_explore[9]
  observed_tox_rate <- total_DLT/total_enrolled
  if(observed_tox_rate <= lambda_e ){
    reccomendation <- "E"
  } else if(observed_tox_rate > lambda_e & observed_tox_rate <= lambda_d ){
    reccomendation <- "S"
  } else if(observed_tox_rate > lambda_d & observed_tox_rate <= lambda_r){
    reccomendation <- "D+R"
  } else if(observed_tox_rate > lambda_r){
    reccomendation <- "D"
  }
  
}
