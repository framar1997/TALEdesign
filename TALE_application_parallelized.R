library(shiny)
library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(shinythemes)
library(gtable)
library(knitr)
library(kableExtra)
library(rmdformats)
library(shinyMatrix)
library(DT)
library(foreach)
library(doParallel)
library(npreg)
library(parallel)


#filesSources = paste0(".../RSim/",list.files(".../RSim/"))
#sapply(filesSources, source)




#


#' This function executes a complete trial simulation
#'
#' @param N: Enrolled patients (must be multiple of n_cohort)
#' @param n_cohort: Cohort size
#' @param toxicity_matrix: Matrix with the probabilities a toxicity
#' @param n_max: Maximum number of patients to assign at DL (i,j), for all i and for all j

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


ConductTrial <- function(N, n_cohort, toxicity_matrix, n_max, select_lambda = TRUE, lambda_e = 0.17, lambda_d = 0.25, lambda_r = 0.33){
  options(dplyr.summarise.inform = FALSE)
  # Return a warning message if N is not multiple of n_cohort
  if( N %% n_cohort != 0 ) warning('N must be a multiple of n_cohort')
  #if( n_max %% n_cohort != 0 ) warning('n_max must be a multiple of n_cohort')
  if( n_max < n_cohort ) warning('n_max must be > n_cohort')
  if(select_lambda == FALSE){
    lambda_e = NULL; lambda_d = NULL; lambda_r = NULL
  }

  ## initialization: start from dose (1,1)
  trial_data = data.frame(step =1,
                          indA=1,
                          indB =1,
                          enrolledPatient = NA,
                          DLT =NA,
                          reccomendation = NA)
  Stop = NULL

  ## list of explorable DLs at step 1
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

    trial_data <- UpdateTrialData(explorable_DL, toxicity_matrix, trial_data, n_cohort, trial_stage, select_lambda,
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
    ## 2: Start build the list of explorable DLs for next the next step of the trial
    ########################################################################################################

    # check for escalation, de-escalation, stay and revisit

    all_DLs$explore <- apply(all_DLs[,c(1,2)], 1, function(y) IsExplorable(y, trial_data, all_DLs, N, n_cohort, trial_stage, nrow, ncol))

    ####### Stop if the list of explorable DLs is empty
    if(sum(all_DLs$explore)==0)
    {
      Stop <- "List of explorable DLs empty before pruning"
      break() # break if the list of explorable DLs for the next step is empty

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
    }

    # Final update of trial_data prior to go to the next step
    trial_stage                = trial_stage   +1
    u <- cbind(rep(trial_stage, nrow(explorable_DL)), explorable_DL,rep(NA, nrow(explorable_DL)),rep(NA, nrow(explorable_DL)),rep(NA, nrow(explorable_DL)) )
    names(u) = names(trial_data)


    trial_data <- rbind(trial_data,u)
  }

  out = list()
  #out$trial_data <- trial_data
  trial_data_summary <- as.data.frame(trial_data %>%
                                        group_by(indA,indB) %>%
                                        summarise(TOT_enrolledPatients = sum(as.integer(enrolledPatient)),
                                                  TOT_DLT = sum(as.integer(DLT))))
  trial_data_summary$TOT_DLT[trial_data_summary$TOT_enrolledPatients == 0] <- NA
  out$trial_data_summary <- trial_data_summary
  g <- matrix(NA, nrow, ncol)
  colnames(g) <- 1:ncol
  rownames(g) <- rev(1:nrow)
  for( i in 1:nrow(trial_data_summary)){
    ind1 <- which(colnames(g) == trial_data_summary[i,1])
    ind2 <- which(rownames(g) == trial_data_summary[i,2])
    g[ind2, ind1] <- trial_data_summary[i,4]/trial_data_summary[i,3]
  }
  out$rates <- round(g,2)
  out$stop <- Stop

  return(out)
}

#

# dose(i,j) allowable states
# E = ESCALATION ALLOWED
# S = STAY AT DOSE
# D+R = DEESCALATE WITH REVISTING
# D = DEESCALATE WITHOUT REVISTING

# Table 2 functions to determine allowable states
# n=3

n3 <- function(DLT) {
  #if(is.integer(DLT)==FALSE) stop("noninteger value")
  ifelse(DLT==0, 'E',
         ifelse(DLT==1, 'D+R',
                ifelse(DLT>1, 'D', NA)))
}
# n=4
n4 <- function(DLT) {
  #if(is.integer(DLT)==FALSE) stop("noninteger value")
  ifelse(DLT == 0, 'E',
         ifelse(DLT == 1, 'S',
                ifelse(DLT > 1, 'D', NA)))
}

# n=5
n5 <- function(DLT) {
  #if(is.integer(DLT)==FALSE) stop("noninteger value")
  ifelse(DLT == 0, 'E',
         ifelse(DLT == 1, 'S',
                ifelse(DLT == 2, 'D+R',
                       ifelse(DLT > 2, 'D', NA))))
}

# n=6
n6 <- function(DLT) {
  #if(is.integer(DLT)==FALSE) stop("noninteger value")
  ifelse(DLT < 1, 'E',
         ifelse(DLT==1, 'S',
                ifelse(DLT == 2, 'D+R',
                       ifelse(DLT > 2, 'D', NA))))
}

# n=7
n7 <- function(DLT) {
  #if(is.integer(DLT)==FALSE) stop("noninteger value")
  ifelse(DLT < 2, 'E',
         ifelse(DLT == 2, 'S',
                ifelse(DLT == 3, 'D+R',
                       ifelse(DLT > 3, 'D', NA))))
}

# n=7
n8 <- function(DLT) {
  #if(is.integer(DLT)==FALSE) stop("noninteger value")
  ifelse(DLT < 2, 'E',
         ifelse(DLT == 2, 'S',
                ifelse(DLT == 3, 'D+R',
                       ifelse(DLT > 3, 'D', NA))))
}

n9 <- function(DLT) {
  #if(is.integer(DLT)==FALSE) stop("noninteger value")
  ifelse(DLT < 3, 'E',
         ifelse(DLT >= 3, 'D', NA))
}


#' function for the permissiable action allowed in dose(ii,jj)
#'
#' @param row_explore: a row of the dataset "explore"
#' Inside the function:
#' @param nn: total enrolled patients at DL (i,j)
#' @param ii: current i dose
#' @param jj: current j dose
#' @param DLT:  DLT at (ii,jj)
#'
#' @export
DoseState <- function(row_explore) {
  nn <- row_explore[8]
  DLT <- row_explore[9]#$totalDLT
  tmp.st <- ifelse(nn==3, n3(DLT),
                   ifelse(nn==4, n4(DLT),
                          ifelse(nn==5, n5(DLT),
                                 ifelse(nn==6, n6(DLT),
                                        ifelse(nn==7, n7(DLT),
                                               ifelse(nn==8, n8(DLT),
                                                      ifelse(nn>=9, n9(DLT))))))))
  out <- toString(tmp.st)
  return(out)
}
#

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

#
HasOneHigherRedOrange<- function(DL, trial_data){
  ind1 <- DL[1]
  ind2 <- DL[2]
  #DL (i,j) is unexplored. Check if we can escalate to DL
  l1 <- DL + c(1,0)
  l2 <- DL + c(0,1)

  redorange1 <- FALSE
  redorange2 <- FALSE
  if(nrow(trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2],]) != 0){
    last_step_l1 <- max(trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2],]$step)
    reccomendation_l1 <- trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2] & trial_data$step == last_step_l1,]$reccomendation
    if(reccomendation_l1 == "D" |reccomendation_l1 == "D+R"){
      redorange1 <- TRUE
    }
  }

  if(nrow(trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2],]) != 0){
    last_step_l2 <- max(trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2],]$step)
    reccomendation_l2 <- trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2] & trial_data$step == last_step_l2,]$reccomendation
    if(reccomendation_l2 == "D" |reccomendation_l2 == "D+R"){
      redorange2 <- TRUE
    }
  }

  if(redorange1 | redorange2){
    return(TRUE)
  } else{
    return(FALSE)
  }

}
#

HasTwoLowerGreen<- function(DL, trial_data){
  ind1 <- DL[1]
  ind2 <- DL[2]
  #DL (i,j) is unexplored. Check if we can escalate to DL
  l1 <- DL - c(1,0)
  l2 <- DL - c(0,1)

  green1 <- FALSE
  green2 <- FALSE

  if(l1[1] < 1 | l1[2] < 1){ #this DL doesn't exists
    green1 <- TRUE
  }else if(nrow(trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2],]) != 0){
    last_step_l1 <- max(trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2],]$step)
    reccomendation_l1 <- trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2] & trial_data$step == last_step_l1,]$reccomendation
    if(reccomendation_l1 == "E"){
      green1 <- TRUE
    }
  }

  if(l2[1] < 1 | l2[2] < 1){
    green2 <- TRUE
  }else if(nrow(trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2],]) != 0){
    last_step_l2 <- max(trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2],]$step)
    reccomendation_l2 <- trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2] & trial_data$step == last_step_l2,]$reccomendation
    if(reccomendation_l2 == "E"){
      green2 <- TRUE
    }
  }

  if(green1 & green2){
    return(TRUE)
  } else{
    return(FALSE)
  }

}
#

IsExplorable<- function(DL, trial_data, all_DLs, N,n_cohort, trial_stage, nrow, ncol){
  ind1 <- DL[1]
  ind2 <- DL[2]
  explore <- FALSE
  reccomendation <- NULL

  # Check if (i,j) has already been explored. If yes, obtain its latest reccomendation
  if(nrow(trial_data[trial_data$indA == ind1 & trial_data$indB == ind2,]) != 0){
    last_step <- max(trial_data[trial_data$indA == ind1 & trial_data$indB == ind2,]$step)
    reccomendation <- trial_data[trial_data$indA == ind1 & trial_data$indB == ind2 & trial_data$step == last_step,]$reccomendation
  }

  # Check if
  two_green <- HasTwoLowerGreen(DL, trial_data)  # TRUE if (i,j) has two lower green neighbors
  one_orangered <- HasOneHigherRedOrange(DL, trial_data) # TRUE if (i,j) has at least one orange or red higher neighbor

  #ESCALATION TO DL (i,j) IS ALLOWED IF:

  if(is.null(reccomendation)){ # (i,j) is unexplored
    if(two_green){ # (i,j) has two green lower adjacent neighbours
      explore <- TRUE
    }
  }


  # STAY at DL (i,j) IF ITS LATEST STATE IS S:

  if(!is.null(reccomendation)){
    if(reccomendation == "S"){
      explore <- TRUE
    }

  }

  # DE-ESCALATE to DL (i,j) IF:

  if(one_orangered){ # (i,j) has at least one orange or red higher adjacent neighbour
    if(two_green){ # (i,j) has two green lower adjacent neighbour
      explore <- TRUE
    }
  }
  # RECONSIDER DL (i,j) IF:

  if(!is.null(reccomendation)){
    if(reccomendation == "D+R"){ # (i,j) is orange
      if(two_green){ # (i,j) has two green lower adjacent neighbours
        if(all_DLs[all_DLs[,1]==ind1 & all_DLs[,2] == ind2,3]==FALSE){ # (i,j) is not in the current list of explorable DLs
          explore <- TRUE
        }
      }
    }

  }

  return(explore)

}

#

# Function which performs isotonic regression over the output of a trial

PerformIsotonic <- function(outputdata, nrow, ncol){
  DLT_mat = matrix(0.5, nrow, ncol) #input$ncol,nrow
  N_mat = matrix(0.1, nrow, ncol)
  tds <- na.omit(outputdata)
  #tds$prop <- tds$TOT_DLT/tds$TOT_enrolledPatients
  for(i in 1:NROW(tds)){
    ii = which( c(1:nrow) == tds[i,2])
    jj = which( c(1:ncol) == tds[i,1])
    DLT_mat[ii,jj] = tds$TOT_DLT[i]/tds$TOT_enrolledPatients[i]
    N_mat[ii,jj] = tds$TOT_enrolledPatients[i]
  }

  isoReg = suppressWarnings(Iso::biviso(DLT_mat, w = N_mat))
  d1 <- expand.grid(x = 1:nrow(isoReg), y = 1:ncol(isoReg))
  out <- transform(d1, z = isoReg[as.matrix(d1)])
  summary_iso <- out[order(out$x), c(2,1,3) ]
  colnames(summary_iso) <- c("indA", "indB", "rates")
  return(summary_iso)
}

#
#' This function applies the pruning procedure to the list of explorable DLs
#'
#' @param trial_data: data frame summarizing the event of the trials.
#' @param explorable_DL: a list of explorable DLs to be pruned
#' @param return a pruned list of explorable DLs
#'
PruneExplorableDLs <- function(trial_data, explorable_DL, n_max){
  # a) Remove overly-toxic DLs
  # Look for DLs with reccomendation = "D": remove these DLs and higher DLs from the list of explorable_DLs
  indeces_to_remove <- trial_data[which(trial_data$reccomendation=="D"), c(2,3)]
  if(nrow(indeces_to_remove >0)){
    for(i in 1:nrow(indeces_to_remove)){
      ind1 <- indeces_to_remove[i,1]
      ind2 <- indeces_to_remove[i,2]
      to_remove <- which(explorable_DL$indA >= ind1 & explorable_DL$indB >= ind2)
      if(sum(to_remove) != 0){
        explorable_DL <- explorable_DL[-to_remove, ]
      }
    }
  }
  # b)  a DL (i,j) is removed from the list if n_max patients have already been treated at that DL
  u = as.data.frame(trial_data %>%
                      group_by(indA,indB) %>%
                      summarise(TOT_enrolledPatients = sum(as.integer(enrolledPatient) ), TOT_DLT= sum(DLT)) )
  indeces_to_remove2 <- u[which(u$TOT_enrolledPatients >= n_max), c(1,2)]
  if(nrow(indeces_to_remove2) >0){
    for(i in 1:nrow(indeces_to_remove2)){
      ind1 <- indeces_to_remove2[i,1]
      ind2 <- indeces_to_remove2[i,2]
      to_remove <- which(explorable_DL$indA == ind1 & explorable_DL$indB == ind2)
      if(sum(to_remove) != 0){
        explorable_DL <- explorable_DL[-to_remove, ]
      }
    }
  }
  # c) a DL (i,j) is removed from the list if a higher DL is in the list
  indexToRemove = c()
  for(i in 1:NROW(explorable_DL))
  {
    if(any(explorable_DL[-i,1] >= explorable_DL[i,1] &  explorable_DL[-i,2] >= explorable_DL[i,2]))
    {
      indexToRemove = c(indexToRemove,i)
    }
  }
  if(is.null(indexToRemove) == F)
  {
    explorable_DL <- explorable_DL[setdiff(1:NROW(explorable_DL), indexToRemove ),]
  }
  return(explorable_DL)
}
#
# Function which takes as imput the data_summary (final output) of a trial and returns the MTD

SelectMTD <- function(data_summary, MTD){
  #data_summary$rates <- data_summary$TOT_DLT/data_summary$TOT_enrolledPatients
  data_summary$MTDdist <- (data_summary$rates - MTD)^2
  minimum <- min(data_summary$MTDdist)
  candidates <- data_summary[data_summary$MTDdist == minimum, ]
  if(nrow(candidates) == 1){
    selected_mtd <- candidates
  } else{

    # if there are multiple MTD candidates, we need to select the lowest dose

    # se sommo indA e indB ottengo dei "livelli di bassezza". Cioé (1,1) somma a 2 ed é il piú basso DL,
    # (1,2) (2,2) (2,1) sommano a 3 e sono i secondi piú bassi DL. Quindi sfrutto ció.

    candidates$indexSum <- candidates$indA + candidates$indB
    minimum2 <- min(candidates$indexSum)
    new_candidates <- candidates[candidates$indexSum == minimum2,]
    if(nrow(new_candidates) == 1){
      selected_mtd <- new_candidates
    } else{
      #if the candidates MTD have the same "level of lowness" then chose one randomly
      random <- sample(1:nrow(new_candidates), 1, replace=F)
      selected_mtd <- new_candidates[random,]
    }
  }
  return(selected_mtd)
}

#
#' This function explores the explorable DLs and it updates the trial information
#' @param explorable_DL: the list of the DLs explorable at the current step
#' @param toxicity_matrix: data frame summarizing the event of the trials.
#' @param trial_data: data frame summarizing the event of the trials.
#' @param n_cohort:
#' @param trial_stage:
#' @param select_lambda:
#' @param return an updated version of trial_data


UpdateTrialData <- function(explorable_DL, toxicity_matrix, trial_data, n_cohort, trial_stage, select_lambda,
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

  explore$step <- trial_stage


  # update trial_data with the new information
  trial_data <- rbind(trial_data[trial_data$step != trial_stage,], explore[,c(11, 1,2,5,4,10)] )

  #out <- list()
  #out$explore <- explore
  #out$trial_data <- trial_data
  return(trial_data)

}

#' GetToxicity: this function extracts the probability of toxicity of DL (i,j) from the toxicity matrix
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


#
SimulateTrial<- function(input,output){

  observeEvent(input$simulate, {

    #set.seed(input$set_seed)

    tox_table <- input$toxicity_matrix
    colnames(tox_table) = c(1:ncol(tox_table))
    rownames(tox_table) = rev(c(1:nrow(tox_table)))
    n_sim = input$n_sim
    sim = n_sim
    lambda_e = input$lambda_e
    lambda_d = input$lambda_d
    lambda_r = input$lambda_r
    n_max = input$n_max
    n_cohort = input$n_cohort
    nrow = input$nrow
    ncol = input$ncol
    inputMTD = input$MTD
    number_of_cohorts =  input$number_of_cohorts
    l = list()
    n_dlt<- c()
    n_enrolled <- c()
    stop <- c()
    l2 <- matrix(0, nrow, ncol)
    # Perform input$n_sim independent trials


    out = list()
    iso = list()
    all_mtd = list()

    #max_workers <- 1  # Set a reasonable maximum limit

    # Determine the number of available cores dynamically
    available_cores <- detectCores()
    #showNotification(paste(available_cores))
    # Calculate the number of workers based on available resources

    cl <- makeCluster(available_cores)


    showNotification(paste("Connecting to ", available_cores, " clusters . . ."))

    clusterEvalQ(cl, {
      library(shiny)
      library(tidyr)
      library(dplyr)
      library(ggplot2)
      library(grid)
      library(gridExtra)
      library(shinythemes)
      library(gtable)
      library(knitr)
      library(kableExtra)
      library(rmdformats)
      library(shinyMatrix)
      library(DT)
      library(foreach)
      library(doParallel)
      library(npreg)
      library(parallel)
    })

    clusterExport(cl, c("tox_table", "n_sim", "lambda_e", "lambda_d", "lambda_r", "n_max", "n_cohort",
                        "nrow", "ncol", "inputMTD", "number_of_cohorts", "l", "n_dlt", "n_enrolled",
                        "stop", "l2","out", "iso","all_mtd", "sim"), envir = environment())

    function_names <- c("IsMatrixMonotonicitySatisfied", "ConductTrial", "PerformIsotonic", "SelectMTD", "PruneExplorableDLs",
                        "GetReccomendation", "HasOneHigherRedOrange", "HasTwoLowerGreen", "IsExplorable", "DoseState",
                        "n3", "n4", "n5", "n6", "n7", "n8", "n9", "UpdateTrialData", "GetToxicity")

    # Export all functions to the clusters
    clusterExport(cl, function_names, envir = environment())


    # Show error if monotonicity assumption is not observed
    if(!IsMatrixMonotonicitySatisfied(input$toxicity_matrix)){
      output$error <- renderText({ "Toxicity probability of each drug while holding
        the other one fixed should be monotonically increasing with dose level."})
      output$result1 = NULL
      output$toxmat = NULL
      output$enrolled = NULL
      output$dlt = NULL
      output$toxrate = NULL
      output$MTDperc = NULL
    }else{


      showNotification(paste( sim, " replications started . . .") )

      #######################################################################################################################################
      ############### Simulation execution
      #######################################################################################################################################

      showNotification("Processing...", duration = NULL, closeButton = FALSE, id = "loading_notification")

      parallel_out <- function(i) {
        set.seed(i)
        # if(i%%50 == 0){
        #   h <- showNotification(paste("Processed more than", i, "repetitions"))
        # }
        out <- ConductTrial(N = n_cohort * number_of_cohorts,
                              n_cohort = n_cohort,
                              toxicity_matrix = tox_table,
                              n_max = n_max, select_lambda = TRUE,
                              lambda_e = lambda_e,
                              lambda_d = lambda_d,
                              lambda_r = lambda_r)

        iso = PerformIsotonic(out$trial_data_summary, nrow, ncol)
        out_iso = list(out, iso)
        return(out_iso)
      }
      clusterExport(cl, "parallel_out", envir = environment())
      out_iso <- parLapply(cl, 1:n_sim, parallel_out)

      removeNotification("loading_notification")

      showNotification(paste( sim, " replications executed . . .") )
      get_out <- function(i){
        out = out_iso[[2]][1]
        return(out[[1]])
      }
      get_iso <- function(i){
        iso = out_iso[[i]][2]
        return(iso[[1]])
      }
      clusterExport(cl, "out_iso", envir = environment())
      clusterExport(cl, "get_out", envir = environment())
      clusterExport(cl, "get_iso", envir = environment())
      out <- parLapply(cl, 1:n_sim, get_out)
      iso <- parLapply(cl, 1:n_sim, get_iso)

      showNotification("Preparing results . . .")

      clusterExport(cl, "out", envir = environment())
      # parallel_iso <- function(i){
      #   iso = PerformIsotonic(out[[i]]$trial_data_summary, nrow, ncol)
      #   return(iso)
      # }
      #
      # clusterExport(cl, "parallel_iso", envir = environment())
      # iso <- parLapply(cl, 1:n_sim, parallel_iso)
      clusterExport(cl, "iso", envir = environment())

      parallel_mtd <- function(i){
        mtd = SelectMTD(iso[[i]], inputMTD)
        return(mtd)
      }
      clusterExport(cl, "parallel_mtd", envir = environment())
      all_mtd <- parLapply(cl, 1:n_sim, parallel_mtd)

      stopCluster(cl)


      for( sim in 1:n_sim){
        mtd <- all_mtd[[sim]]

        # collect results in a common output:
        MTD <- matrix(0, nrow, ncol)
        rownames(MTD) <- rev(1: nrow)
        MTD[which(as.integer(rownames(MTD)) == mtd$indB), mtd$indA] <- 1
        l2 <- l2 + MTD

        summary <- out[[sim]]$trial_data_summary
        n_dlt <- c(n_dlt,sum(summary$TOT_DLT))
        n_enrolled <-c(n_enrolled, sum(summary$TOT_enrolledPatients) )
        stop <- c(stop, out[[sim]]$stop)
        l[[sim]] <- summary
        #if(sim%%50 == 0){
        #  showNotification(paste("Processed", sim, "computations over",input$n_sim))
        #}
      }


      #showNotification( "pure qui")

      l2 <- l2/n_sim*100

      n_dlt <- c(n_dlt,sum(summary$TOT_DLT))
      n_enrolled <-c(n_enrolled, sum(summary$TOT_enrolledPatients) )
      #stop <- c(stop, data$stop)

      # Close the cluster when finished

      # Extract individual elements from results

      #######################################################################################################################################
      ############### OUTPUT:
      #######################################################################################################################################
      # BUILD OUTPUT
      a <- do.call("rbind", l)
      b <- as.data.frame(a %>%group_by(indA,indB) %>%
                           summarise(average_enrolledPatients = mean(as.integer(TOT_enrolledPatients))))
      a <- as.data.frame(a %>%group_by(indA,indB) %>%
                           summarise(total_enrolledPatients = sum(as.integer(TOT_enrolledPatients)), ###########
                                     average_DLT = mean(as.integer(TOT_DLT))))

      a <- cbind( a$indA, a$indB, round(a$total_enrolledPatients,2),
                  round(a$average_DLT,2),
                  round(a$average_DLT/b$average_enrolledPatients,2)) #a$total_enrolledPatients
      colnames(a) <- c("indA", "indB", "Total Enrolled Patients", "Average DLT", "Average Toxicity Rate") ###########
      average_enrolled <- matrix(NA, nrow = nrow, ncol = ncol)

      average_dlt <- matrix(NA, nrow = nrow, ncol = ncol)
      average_toxicity_rate <- matrix(NA, nrow = nrow, ncol = ncol)
      rownames(average_enrolled) <- rev(1: nrow)
      rownames(average_dlt) <- rev(1: nrow)
      rownames(average_toxicity_rate) <- rev(1: nrow)

      for(i in 1:nrow(a)){
        ind1 <- a[i,1]
        ind2 <- a[i,2]
        average_enrolled[which(as.integer(rownames(average_enrolled)) == ind2), ind1] <- a[i,3] # TOTAL ENROLLED PATIENTS
        average_dlt[which(as.integer(rownames(average_dlt)) == ind2), ind1] <- a[i,4]
        average_toxicity_rate[which(as.integer(rownames(average_toxicity_rate)) == ind2), ind1] <- a[i,5]
      }

      ###########
      average_enrolled[is.na(average_enrolled)] = 0
      average_enrolled <- round(average_enrolled/(n_sim*n_cohort*number_of_cohorts)*100,2) #percentage of patients enrolled ###########
      #average_enrolled[is.na(average_enrolled)] = 0
      #average_dlt[is.na(average_dlt)] = 0
      #average_toxicity_rate[is.na(average_toxicity_rate)] = 0

      ########### OUTPUT 1: SIMULATION SETTING

      mtd_tox_table = which(tox_table == inputMTD)
      tot_perc_mtd = 0
      for(i in mtd_tox_table){
        tot_perc_mtd = tot_perc_mtd + l2[i]
      }

      average_enrolled_mtd = 0 ###########
      for(i in mtd_tox_table){###########
        average_enrolled_mtd =  average_enrolled_mtd + as.matrix(average_enrolled)[i]###########
      }###########

      paste("B",rev(1: 5))
      rownames(average_enrolled) <- paste("B",rev(1: nrow), sep = "")
      colnames(average_enrolled) <- paste("A",(1: ncol), sep = "")
      rownames(average_dlt) <- rownames(average_enrolled)
      colnames(average_dlt) <- colnames(average_enrolled)
      rownames(average_toxicity_rate) <- rownames(average_enrolled)
      colnames(average_toxicity_rate) <- colnames(average_enrolled)
      rownames(tox_table) = rownames(average_enrolled)
      colnames(tox_table) = colnames(average_enrolled)

      rownames(l2) = rownames(average_enrolled)
      colnames(l2) = colnames(average_enrolled)



      max_sample_size <- length(stop[stop == "Maximum Sample Size Reached"])
      empty_expl_dl_after_pruning <- length(stop[stop == "List of explorable DLs empty after pruning"])
      empty_expl_dl_before_pruning <- length(stop[stop == "List of explorable DLs empty before pruning"])
      res1 <- matrix(round(c(mean(n_dlt),mean(n_enrolled), tot_perc_mtd,
                             average_enrolled_mtd, ###########
                             max_sample_size/length(stop)*100,
                             empty_expl_dl_after_pruning/length(stop)*100,
                             empty_expl_dl_before_pruning/length(stop)*100),2),7 ,1)
      colnames(res1) <- " "
      rownames(res1) <- c("Average number of total DLT", "Average number of total enrolled patients",
                          "Selection percent of MTD", "% of patients treated at MTD",
                          "% of times in which the maximum sample size was reached",
                          "% of times in which the list of explorable DLs was empty after pruning",
                          "% of times in which the list of explorable DLs was empty before pruning"
      )




      output$result1 = DT::renderDataTable(res1, server = FALSE, class = 'cell-border stripe',
                                           options = list(pageLength = 5, dom = 'tip'),
                                           caption = 'Table 1: Operating characteristics summary')
      output$toxmat = DT::renderDataTable(tox_table, server = FALSE, class = 'cell-border stripe',
                                          options = list(pageLength = 5, dom = 'tip'),
                                          caption = 'Table 2: True toxicity probabilities')
      output$toxrate = DT::renderDataTable(average_toxicity_rate, server = FALSE, class = 'cell-border stripe',
                                           options = list(pageLength = 5, dom = 'tip'),
                                           caption = 'Table 3: Average toxicity rates')
      output$dlt = DT::renderDataTable(average_dlt, server = FALSE, class = 'cell-border stripe',
                                       options = list(pageLength = 5, dom = 'tip'),
                                       caption = 'Table 4: Average dose-limiting toxicities')
      output$enrolled = DT::renderDataTable(average_enrolled, server = FALSE, class = 'cell-border stripe',
                                            options = list(pageLength = 5, dom = 'tip'),
                                            caption = 'Table 5: % of enrolled patients') ###########



      output$MTDperc = DT::renderDataTable(round(l2,2), server = FALSE, class = 'cell-border stripe',
                                           options = list(pageLength = 5, dom = 'tip'),
                                           caption = 'Table 6: % of MTD selection. In each trial, the selectd MTD is the DL
                                         whose simulated toxicity rate is closest to the MTD target in terms of Euclidean distance; if
                                         there are multiple MTD candidates, the MTD is selected by giving priority to the lowest DL when
                                         it is possible; otherwise, the MTD is randomly selected between the candidates.
                                         ')

    } #else
  })


}




#################################

#ui
ui <- fluidPage( theme = shinytheme("lumen"),
                 navbarPage("TALE design",

                            tabPanel("Simulation",
                                     sidebarPanel(tags$head(
                                       tags$style(
                                         HTML(".shiny-notification { position:fixed;
                             top: calc(50%);
                             left: calc(50%);
                             }"
                                         )
                                       )
                                     ),
                                     fluidRow(
                                       column(width = 5,
                                              numericInput("ncol", label = "Drug A doses:", 4, min=1, max = 7,
                                                           step = 1)
                                       ),
                                       column(width = 5,
                                              numericInput("nrow", label = "Drug B doses:", 3, min=1, max = 7,
                                                           step = 1)
                                       )
                                     ),
                                     hr(),
                                     fluidRow(
                                       column(width = 5,
                                              numericInput("number_of_cohorts", label = "Number of cohorts:", 10,
                                                           step = 1)
                                       ),
                                       column(width = 5,
                                              numericInput("n_cohort", label =  "Cohort size:", 3, min=1, max = 9,
                                                           step = 1)
                                       )
                                     ),
                                     textOutput("text"),
                                     hr(),
                                     uiOutput("n_max"),
                                     hr(),
                                     fluidRow(
                                       column(width = 4,
                                              numericInput("lambda_e", label = HTML("&lambda;e:"), 0.16, min = 0, max = 1,
                                                           step = 0.1)
                                       ),
                                       column(width = 4,
                                              numericInput("lambda_d", label =  HTML("&lambda;d:"), 0.24, min = 0, max = 1,
                                                           step = 0.1)
                                       ),
                                       column(width = 4,
                                              numericInput("lambda_r", label =  HTML("&lambda;r:"), 0.45, min = 0, max = 1,
                                                           step = 0.1)
                                       )
                                     ),
                                     hr(),
                                     fluidRow(
                                       column(width = 6,
                                              numericInput("MTD", label = "target MTD:", 0.2, min = 0, max =1, step = 0.01)
                                       ),
                                       column(width = 6,
                                              numericInput("n_sim", label = "Number of replications:", 1000, min= 1, max = 15000, step = 20)
                                       ),
                                       # column(width = 5,
                                       #        numericInput("set_seed", label =  "Set seed:", 3,
                                       #                     step = 1)
                                       # )
                                     ),
                                     hr(),

                                     uiOutput("toxicity_matrix"),

                                     hr(),
                                     actionButton("simulate", "Start simulation"), align = "center"
                                     ),#siderbarPanel



                                     mainPanel(
                                       textOutput("error"),
                                       DT::dataTableOutput('result1'),
                                       hr(),
                                       fluidRow(
                                         column(6, DT::dataTableOutput('toxmat')),
                                         column(6, DT::dataTableOutput('toxrate'))
                                       ),
                                       hr(),
                                       fluidRow(
                                         column(6, DT::dataTableOutput('enrolled')),
                                         column(6, DT::dataTableOutput('dlt'))
                                       ),
                                       DT::dataTableOutput('MTDperc')

                                     )
                            ) #tabPanel
                 ) #navbarPage
)

server <- shinyServer(function(input, output) {

  # render INPUTS
  output$text <- renderText({ paste("The trial sample size is:", as.numeric(input$n_cohort)*as.numeric(input$number_of_cohorts),
                                    sep = " ")})

  output$n_max <- renderUI({
    sliderInput("n_max", "Maximum number of patients enrolled for a DL:", 3*as.numeric(input$n_cohort), min = as.numeric(input$n_cohort),
                max = 9*as.numeric(input$n_cohort), step = as.numeric(input$n_cohort))
  })
  output$toxicity_matrix <- renderUI({
    if(input$nrow == 3 & input$ncol == 4){
      tox_mat <- matrix(rbind(c(0.15,	0.33,	0.37,	0.50),
                              c(0.10,	0.20,	0.33,	0.45),
                              c(0   ,	0.15,	0.20,	0.37)), input$nrow, input$ncol, byrow = FALSE )
      colnames(tox_mat) = paste("A", 1:input$ncol, sep = "")
      rownames(tox_mat) = paste("B", rev(1:input$nrow),  sep = "")
    }else{
      tox_mat <- matrix(0.15, input$nrow, input$ncol)
      colnames(tox_mat) = paste("A", 1:input$ncol, sep = "")
      rownames(tox_mat) = paste("B", rev(1:input$nrow),  sep = "")
    }

    matrixInput("toxicity_matrix", "True toxicities matrix:",
                value = tox_mat,
                class = "numeric" )
  })



  # Simulate trial:
  SimulateTrial(input, output)

}) #server

# Run the application
shinyApp(ui = ui, server = server)



