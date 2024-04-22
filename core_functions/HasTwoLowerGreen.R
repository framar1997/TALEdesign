#' This function checks if a DL has two lower green neighbor.
#'
#' @param DL: a DL
#' @param trial_data: data frame summarizing the event of the trials
#' @param return TRUE if DL has two lower green neighbor
#'

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
      last_stage_l1 <- max(trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2],]$stage)
      reccomendation_l1 <- trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2] & trial_data$stage == last_stage_l1,]$reccomendation
      if(reccomendation_l1 == "E"){
        green1 <- TRUE
      }
    }
    
    if(l2[1] < 1 | l2[2] < 1){
      green2 <- TRUE
    }else if(nrow(trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2],]) != 0){
      last_stage_l2 <- max(trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2],]$stage)
      reccomendation_l2 <- trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2] & trial_data$stage == last_stage_l2,]$reccomendation
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

