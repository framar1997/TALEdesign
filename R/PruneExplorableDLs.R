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

