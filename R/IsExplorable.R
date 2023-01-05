IsExplorable<- function(DL, trial_data, all_DLs, N,n_cohort, trial_stage, nrow, ncol){
  ind1 <- DL[1]
  ind2 <- DL[2]
  explore <- FALSE
  reccomendation <- NULL
  
  # Check if (i,j) has already been explored. If yes, obtain its latest reccomendation
  if(nrow(trial_data[trial_data$indA == ind1 & trial_data$indB == ind2,]) != 0){
    last_stage <- max(trial_data[trial_data$indA == ind1 & trial_data$indB == ind2,]$stage)
    reccomendation <- trial_data[trial_data$indA == ind1 & trial_data$indB == ind2 & trial_data$stage == last_stage,]$reccomendation
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

