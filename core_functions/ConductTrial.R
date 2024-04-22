#' This function simulate a single trial
#'
#' @param N: Enrolled patients (must be multiple of n_cohort)
#' @param n_cohort: Cohort size
#' @param toxicity_matrix: Matrix with the probabilities a toxicity
#' @param n_max: Maximum number of patients to assign at DL (i,j), for all i and for all j
#' @param lambda_e: escalation parameter
#' @param lambda_d: de-escalation parameter
#' @param lambda_r: de-escalate and remove parameter
#' @param return a list of objects containing information on the simulated trial.

ConductTrial <- function(N, n_cohort, toxicity_matrix, n_max, select_lambda = TRUE, lambda_e = 0.17, lambda_d = 0.25, lambda_r = 0.45){
  options(dplyr.summarise.inform = FALSE)
  colnames(toxicity_matrix) = c(1:ncol(toxicity_matrix))
  rownames(toxicity_matrix) = rev(c(1:nrow(toxicity_matrix)))
  # Return a warning message if N is not multiple of n_cohort
  if( N %% n_cohort != 0 ) warning('N must be a multiple of n_cohort')
  #if( n_max %% n_cohort != 0 ) warning('n_max must be a multiple of n_cohort')
  if( n_max < n_cohort ) warning('n_max must be > n_cohort')
  if(select_lambda == FALSE){
    lambda_e = NULL; lambda_d = NULL; lambda_r = NULL
  }
  if(!IsMatrixMonotonicitySatisfied(toxicity_matrix)) warning('Toxicity probability of each drug while holding the other one fixed should be monotonically increasing with dose level.')

  ## initialization: start from dose (1,1)
  trial_data = data.frame(stage =1,
                          indA=1,
                          indB =1,
                          enrolledPatient = NA,
                          DLT =NA,
                          reccomendation = NA)
  Stop = NULL

  ## list of explorable DLs at stage 1
  explorable_DL =data.frame(indA = 1, indB = 1)
  ncol <- ncol(toxicity_matrix)
  nrow <- nrow(toxicity_matrix)
  # All doses:

  all_DLs <- expand.grid(1:ncol, 1:nrow)
  all_DLs$explore <- FALSE

  while(TRUE)
  {

    trial_stage <- trial_data[nrow(trial_data),1]

    ########################################################################################################
    ##### 1: Given the current list of explorable DLs, explore these DLs and update the trialData information
    ########################################################################################################

    trial_data <- SimulateToxicity(explorable_DL, toxicity_matrix, trial_data, n_cohort, trial_stage, select_lambda,
                                   lambda_e, lambda_d, lambda_r)


    ########################################################################################################
    ########## Stop if the maximum sample size is reached (we don't need to know which DL is explorable next)
    ########################################################################################################

    if( sum(trial_data$enrolledPatient) >= N)
    {
      Stop <- "Maximum Sample Size Reached"
      break()

    }

    ########################################################################################################
    ## 2: Start build the list of explorable DLs for next the next stage of the trial
    ########################################################################################################

    # check for escalation, de-escalation, stay and revisit

    all_DLs$explore <- apply(all_DLs[,c(1,2)], 1, function(y) IsExplorable(y, trial_data, all_DLs, N, n_cohort, trial_stage, nrow, ncol))

    ####### Stop if the list of explorable DLs is empty
    if(sum(all_DLs$explore)==0)
    {
      Stop <- "List of explorable DLs empty before pruning"
      break() # break if the list of explorable DLs for the next stage is empty

    } else{ # If it is not empty, remove duplicates
      explorable_DL <- all_DLs[all_DLs$explore == TRUE, c(1,2)]
      colnames(explorable_DL) <- c("indA", "indB")
    }


    ########################################################################################################
    ## 3: Prune the list of explorable DLs
    ########################################################################################################

    explorable_DL <- PruneExplorableDLs(trial_data, explorable_DL, n_max)



    ####### Stop if the list of explorable DLs is empty after the pruning procedure

    if(nrow(explorable_DL) ==0)
    {
      Stop <- "List of explorable DLs empty after pruning"
      break()
    }

    ########################################################################################################
    ## 4: Check if there are enough cohorts for each DL in the list of explorable DLS
    ########################################################################################################

    # If there are not enough cohorts...

    explorable_DL <- AssignCohorts(explorable_DL, N, n_cohort, trial_data)

    # Final update of trial_data prior to go to the next stage
    trial_stage                = trial_stage   +1
    u <- cbind(rep(trial_stage, nrow(explorable_DL)), explorable_DL,rep(NA, nrow(explorable_DL)),rep(NA, nrow(explorable_DL)),rep(NA, nrow(explorable_DL)) )
    names(u) = names(trial_data)


    trial_data <- rbind(trial_data,u)
  }

  out = list()
  out$trial_data <- trial_data
  all_DLs <- data.frame(NA,all_DLs[,c(1,2)],NA, NA, NA)
  colnames(all_DLs) <- colnames(trial_data)
  trial_data_summary <- as.data.frame(rbind(trial_data, all_DLs) %>%
                                        group_by(indA,indB) %>%
                                        summarise(TOT_enrolledPatients = sum(as.integer(enrolledPatient), na.rm = T ),
                                                  TOT_DLT = sum(as.integer(DLT), na.rm = T)))
  #trial_data_summary$DLT <- ifelse(trial_data_summary$TOT_enrolledPatients == 0, NA, trial_data_summary$DLT)
  trial_data_summary$TOT_DLT[trial_data_summary$TOT_enrolledPatients == 0] <- NA
  out$trial_data_summary <- trial_data_summary

  out$stop <- Stop
  return(out)
}
