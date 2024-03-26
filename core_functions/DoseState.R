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

