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
