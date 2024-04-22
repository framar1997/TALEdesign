# This function selects the MTD at the end of the trial

#' @param data_summary: dataset containing the results of a simulated trial
#' @param MTD: target toxicity rate
#' @param return the selected MTD

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
