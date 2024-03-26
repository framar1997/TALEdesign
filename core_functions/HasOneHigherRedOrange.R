HasOneHigherRedOrange<- function(DL, trial_data){
  ind1 <- DL[1]
  ind2 <- DL[2]
  #DL (i,j) is unexplored. Check if we can escalate to DL
  l1 <- DL + c(1,0)
  l2 <- DL + c(0,1)
 
  redorange1 <- FALSE
  redorange2 <- FALSE
  if(nrow(trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2],]) != 0){
    last_stage_l1 <- max(trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2],]$stage)
    reccomendation_l1 <- trial_data[trial_data$indA == l1[1] & trial_data$indB == l1[2] & trial_data$stage == last_stage_l1,]$reccomendation
    if(reccomendation_l1 == "D" |reccomendation_l1 == "D+R"){
      redorange1 <- TRUE
    }
  }
  
  if(nrow(trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2],]) != 0){
    last_stage_l2 <- max(trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2],]$stage)
    reccomendation_l2 <- trial_data[trial_data$indA == l2[1] & trial_data$indB == l2[2] & trial_data$stage == last_stage_l2,]$reccomendation
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
