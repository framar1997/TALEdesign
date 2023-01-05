#' This function applies the pruning procedure to the list of explorable DLs
#'
#' @param trial_data: data frame summarizing the event of the trials.
#' @param explorable_DL: a list of explorable DLs to be pruned
#' @param return a pruned list of explorable DLs

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
