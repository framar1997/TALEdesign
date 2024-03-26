#' This function explores the explorable DLs and it updates the trial information
#' @param explorable_DL: the list of the DLs explorable at the current stage
#' @param toxicity_matrix: data frame summarizing the event of the trials.
#' @param trial_data: data frame summarizing the event of the trials.
#' @param n_cohort:
#' @param trial_stage:
#' @param select_lambda:
#' @param return an updated version of trial_data


SimulateToxicity <- function(explorable_DL, toxicity_matrix, trial_data, n_cohort, trial_stage, select_lambda,
                       lambda_e, lambda_d, lambda_r){
  lambdas = c(lambda_e, lambda_d, lambda_r)
  #get  toxicity of explorable DLs
  tox <- apply(explorable_DL, 1, function(x) GetToxicity(as.integer(x[1]), as.integer(x[2]), toxicity_matrix) )

  # Start build "explore"
  explore <- cbind(explorable_DL, tox) #DLs to explore and respective toxicities
  explore$DLT <- unlist(lapply(explore$tox, rbinom, n = 1, size = n_cohort )) # Observed DLTs
  explore$enrolledPatient <- n_cohort
  #  Obtain the past information over the DLs in "explore" (enrolled patients & DLT)
  trial_dataSummary <- as.data.frame(trial_data[!is.na(trial_data$enrolledPatient),] %>%
                                      group_by(indA,indB) %>%
                                      summarise(pastEnrolledPatients = sum(as.integer(enrolledPatient) ), pastDLT= sum(DLT)) )
  explore <- left_join(explore, trial_dataSummary, by = c("indA", "indB"))
  explore$totalEnrolledPatients <- apply(explore[,c(5,6)], 1, sum, na.rm=T)
  explore$totalDLT <- apply(explore[,c(4,7)], 1, sum, na.rm=T)

  # obtain the new suggestions
  if(select_lambda == TRUE){
    explore$reccomendation = apply(explore, 1, function(y) GetReccomendation(y, lambda_e, lambda_d, lambda_r))
  }else if(select_lambda == FALSE){
    explore$reccomendation = apply(explore, 1, function(y) DoseState(y))
  }

  explore$stage <- trial_stage


  # update trial_data with the new information
  trial_data <- rbind(trial_data[trial_data$stage != trial_stage,], explore[,c(11, 1,2,5,4,10)] )

  #out <- list()
  #out$explore <- explore
  #out$trial_data <- trial_data
  return(trial_data)

}


