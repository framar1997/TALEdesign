#' This function explores the explorable DLs and it updates the trial information
#'
#' @param explorable_DL: the list of the DLs explorable at the current stage
#' @param N: Enrolled patients (must be multiple of n_cohort)
#' @param n_cohort: Cohort size
#' @param trial_data: data frame summarizing the event of the trials
#' @param return a list of DLs that can be explored in the next stage of the trial


AssignCohorts <- function(explorable_DL, N, n_cohort, trial_data){
  if(nrow(explorable_DL)> (N - sum(trial_data$enrolledPatient, na.rm=T))/n_cohort){
    cohorts_left <- (N - sum(trial_data$enrolledPatient, na.rm=T))/n_cohort
    # assign the cohorts left giving priority to the DLs with a lower amount of enrolled patients
    # when the amount of enrolled patients is equal, the DLs which are mantained in the list of explorable DLs
    # are randomly choosen

    u <- left_join(explorable_DL, trial_data[,c(2,3,4)], by = c("indA", "indB"))
    u[which(is.na(u$enrolledPatient)),3] = 0
    u = as.data.frame(u %>%
                        group_by(indA,indB) %>%
                        summarise(TOT_enrolledPatients = sum(as.integer(enrolledPatient))))
    # Randomize the DLs in the list
    u <- u[sample(nrow(u)), ]

    # Order the list by the total number of enrolled patients
    u <-u[order(u$TOT_enrolledPatients),]

    # Maintain in the list only the first cohorts_left rows of the dataset
    explorable_DL <- u[1:cohorts_left, c(1,2)]
  } else{
    # nothing
  }
  return(explorable_DL)
}
