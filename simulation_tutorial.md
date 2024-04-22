# Simulation tutorial
Francesco Mariani

## Description

In this tutorial we provide guidelines to simulate a trial using the
provided R code. In what it follows, we consider scenarios extrapolated
from the section “simulation studies”.

We start uploading useful libraries,

``` r
require(tidyr)
require(dplyr)
require(Iso)
require(gridExtra)
require(grid)
require(parallel)
require(stringr)
```

we set seed to ensure reproducibility of results,

``` r
set.seed(3)
```

and finally we upload the R functions we are going to use to simulate a
trial,

``` r
sapply( list.files("core_functions", full.names=TRUE), source )
```

## Simulate a single trial:

In order to simulate a single trial, we need to set the TALE design
parameters, which are:

-   `N`: maximum number of patients;

-   `n_cohort`: size of the cohort;

-   `n_max`: maximum number of patients which can be enrolled at a
    certain dose level;

-   `toxicity_matrix`: the matrix of the true toxicity probabilities.
    Columns and rows represent doses of drug A and drug B, respectively;

-   `ncol`: number of doses of drug A

-   `nrow`: number of doses of drug B

-   `lambda_e`: escalation parameter;

-   `lambda_d`: de-escalation parameter;

-   `lambda_r`: de-escalation and remove parameter;

We set these parameters following **Scenario 1** of Section “numerical
study”. We want to perform a simulated trial where 30 patients are
divided in cohorts of 3 patients. The dose levels are 25 combinations of
5 doses of drug A and 5 doses of drug B. We do not want to assign more
than 9 patients to a certain dose level.

The true toxicity probabilities of the dose levels are assumed to be:

<table>
<thead>
<tr class="header">
<th></th>
<th>A 1</th>
<th>A 2</th>
<th>A 3</th>
<th>A 4</th>
<th>A 5</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><strong>B 5</strong></td>
<td>0.32</td>
<td>0.37</td>
<td>0.42</td>
<td>0.43</td>
<td>0.46</td>
</tr>
<tr class="even">
<td><strong>B 4</strong></td>
<td>0.27</td>
<td>0.33</td>
<td>0.37</td>
<td>0.40</td>
<td>0.43</td>
</tr>
<tr class="odd">
<td><strong>B 3</strong></td>
<td>0.24</td>
<td>0.29</td>
<td>0.34</td>
<td>0.37</td>
<td>0.42</td>
</tr>
<tr class="even">
<td><strong>B 2</strong></td>
<td>0.20</td>
<td>0.23</td>
<td>0.29</td>
<td>0.33</td>
<td>0.37</td>
</tr>
<tr class="odd">
<td><strong>B 1</strong></td>
<td>0.11</td>
<td>0.20</td>
<td>0.24</td>
<td>0.28</td>
<td>0.32</td>
</tr>
</tbody>
</table>

``` r
N = 42
n_cohort = 3
n_max = 9
toxicity_matrix = matrix(rbind(c(0.32,  0.37,   0.42,   0.43,   0.46),
                               c(0.27,  0.33,   0.37,   0.40,   0.43),
                               c(0.24,  0.29,   0.34,   0.37,   0.42),
                               c(0.20,  0.23,   0.29,   0.33,   0.37),
                               c(0.11,  0.20,   0.24,   0.28,   0.32)),
                               nrow = 5, ncol = 5, byrow = FALSE )
lambda_e = 0.16
lambda_d = 0.24
lambda_r = 0.45
nrow = 5
ncol = 5
```

At this point, we can simulate a single trial by using the function
`ConductTrial()`; we save results in an object called `trial_results.`

``` r
trial_results <- ConductTrial(N = N,
                              n_cohort = n_cohort,
                              toxicity_matrix = toxicity_matrix,
                              n_max = n_max,
                              lambda_e = 0.16, lambda_d = 0.24, lambda_r = 0.45)
```

`trial_results` is a list of three objects:

-   `trial_results$trial_data`: contains step-by-step details about the
    simulated trial;

-   `trial_results$trial_data_summary`: report total enrolled patients
    and total observed dose limiting toxicities for each dose level;

-   `trial_results$stop`: contains the reasons why the trial stopped.

The function `PerformIsotonic()` is used to perform the Isotonic
regression at the end of the trial. The function requires as input the
objects `trial_results$trial_data_summary`, `nrow` and `ncol`.

``` r
isotonic_results = PerformIsotonic(trial_results$trial_data_summary,
                                   nrow, ncol)
```

The output `isotonic_results` contains the estimates of the Isotonic
regression based on the accrued trial data.

Finally, let us assume that the MTD has 0.20 of true probability of
toxicity. We set the value `trueMTD` at 0.20, and we use the function
SelectMTD() to obtain the DL that the simulated trial identifies as MTD.

``` r
trueMTD = 0.20

SelectMTD(isotonic_results, trueMTD)
```

      indA indB     rates      MTDdist
    6    2    1 0.2222222 0.0004938272

## Replicate the trial multiple times to conduct a simulation study:

Suppose now we want to conduct a simulation study, thus we want to
replicate multiple times the trial described in the previous section.
The aim of the simulation study is to evaluate the operative
characteristics of the TALE design.

First, we set the number of replications, `n_sim`, we want to conduct:

``` r
n_sim = 10
```

To accelerate computations, we execute the independent trials in
parallel. We start specifying the number of clusters we are going to
use, and we send to the clusters the necessary libraries and core
functions.

``` r
# Connect to the clusters
available_cores <- detectCores() # detect available clusters
cl <- makeCluster(available_cores) # turn on the clusters

# Upload needed packages on the clusters
clusterEvalQ(cl, {
      library(tidyr)
      library(dplyr)
      library(Iso)
      library(grid)
      library(gridExtra)
      library(parallel)
    })


# Export needed functions and objects to the clusters
function_names = str_sub(list.files("core_functions", full.names=FALSE), end = -3)

clusterExport(cl, as.vector(function_names), envir = environment()) 
clusterExport(cl, c("N", "n_cohort", "n_max", "toxicity_matrix",
                    "trueMTD", "nrow", "ncol"),
              envir = environment())
```

The following piece of code executes the trials and save results.

``` r
parallel_out <- function(i) {
        set.seed(i)
        out <- ConductTrial(N = N,
                              n_cohort = n_cohort,
                              toxicity_matrix = toxicity_matrix,
                              n_max = n_max,
                              lambda_e = 0.16, lambda_d = 0.24, lambda_r = 0.45)
        
        iso = PerformIsotonic(out$trial_data_summary, nrow, ncol)
        out_iso = list(out, iso)
        return(out_iso)
      }
clusterExport(cl, "parallel_out", envir = environment())
out_iso <- parLapply(cl, 1:n_sim, parallel_out)

get_out <- function(i){
        out = out_iso[[2]][1]
        return(out[[1]])
      }
get_iso <- function(i){
        iso = out_iso[[i]][2]
        return(iso[[1]])
      }
clusterExport(cl, "out_iso", envir = environment())
clusterExport(cl, "get_iso", envir = environment())
clusterExport(cl, "get_out", envir = environment())
out <- parLapply(cl, 1:n_sim, get_out)
iso <- parLapply(cl, 1:n_sim, get_iso)
      
clusterExport(cl, "out", envir = environment())
clusterExport(cl, "iso", envir = environment())
parallel_mtd <- function(i){
  mtd = SelectMTD(iso[[i]], trueMTD)
  return(mtd)
  }
clusterExport(cl, "parallel_mtd", envir = environment())
all_mtd <- parLapply(cl, 1:n_sim, parallel_mtd)


# Close the clusters
stopCluster(cl)
```

The following piece of code is used to manipulate results in order to
let them be easy to interpret.

``` r
l = list()
n_dlt<- c()
n_enrolled <- c()
stop <- c()
l2 <- matrix(0, nrow, ncol)
for( sim in 1:n_sim){
        mtd <- all_mtd[[sim]]
        
        # collect results in a common output:
        MTD <- matrix(0, nrow, ncol)
        rownames(MTD) <- rev(1: nrow)
        MTD[which(as.integer(rownames(MTD)) == mtd$indB), mtd$indA] <- 1
        l2 <- l2 + MTD
        
        summary <- out[[sim]]$trial_data_summary
        n_dlt <- c(n_dlt,sum(summary$TOT_DLT[is.na(summary$TOT_DLT) == F]))
        n_enrolled <-c(n_enrolled, sum(summary$TOT_enrolledPatients) )
        stop <- c(stop, out[[sim]]$stop)
        l[[sim]] <- summary
}

l2 <- l2/n_sim*100
n_dlt <- c(n_dlt,sum( na.omit(summary$TOT_DLT)) )
n_enrolled <-c(n_enrolled, sum(summary$TOT_enrolledPatients) )

a <- do.call("rbind", l)
b <- as.data.frame(a %>%group_by(indA,indB) %>%
                   summarise(average_enrolledPatients = 
                               mean(as.integer(TOT_enrolledPatients))))
a <- as.data.frame(a %>%group_by(indA,indB) %>%
                   summarise(total_enrolledPatients = 
                               sum(as.integer(TOT_enrolledPatients)), 
                                     average_DLT = mean(as.integer(TOT_DLT))))
      
a <- cbind( a$indA, a$indB, round(a$total_enrolledPatients,2),
                  round(a$average_DLT,2),
                  round(a$average_DLT/b$average_enrolledPatients,2)) 
colnames(a) <- c("indA", "indB", "Total Enrolled Patients",
                 "Average DLT", "Average Toxicity Rate") 

average_enrolled <- matrix(NA, nrow = nrow, ncol = ncol)
average_dlt <- matrix(NA, nrow = nrow, ncol = ncol)
average_toxicity_rate <- matrix(NA, nrow = nrow, ncol = ncol)
rownames(average_enrolled) <- rev(1: nrow)
rownames(average_dlt) <- rev(1: nrow)
rownames(average_toxicity_rate) <- rev(1: nrow)
      
for(i in 1:nrow(a)){
        ind1 <- a[i,1]
        ind2 <- a[i,2]
        average_enrolled[which(as.integer(rownames(average_enrolled)) ==
                                 ind2),ind1] <- a[i,3] 
        average_dlt[which(as.integer(rownames(average_dlt)) ==
                            ind2), ind1] <- a[i,4]
        average_toxicity_rate[which(as.integer(rownames(average_toxicity_rate))
                                    == ind2), ind1] <- a[i,5]
      }
      
average_enrolled[is.na(average_enrolled)] = 0
average_enrolled <- round(average_enrolled/(n_sim*N)*100,2) 
      
mtd_tox_table = which(toxicity_matrix == trueMTD)
tot_perc_mtd = 0
for(i in mtd_tox_table){
  tot_perc_mtd = tot_perc_mtd + l2[i]
  }
      
average_enrolled_mtd = 0
for(i in mtd_tox_table){
  average_enrolled_mtd =  average_enrolled_mtd + as.matrix(average_enrolled)[i]
  }

rownames(average_enrolled) <- paste("B",rev(1: nrow), sep = "")
colnames(average_enrolled) <- paste("A",(1: ncol), sep = "")
rownames(average_dlt) <- rownames(average_enrolled)
colnames(average_dlt) <- colnames(average_enrolled)
rownames(average_toxicity_rate) <- rownames(average_enrolled)
colnames(average_toxicity_rate) <- colnames(average_enrolled)
rownames(toxicity_matrix) = rownames(average_enrolled)
colnames(toxicity_matrix) = colnames(average_enrolled)
rownames(l2) = rownames(average_enrolled)
colnames(l2) = colnames(average_enrolled)
      
      
max_sample_size <- length(stop[stop == "Maximum Sample Size Reached"])
empty_expl_dl_after_pruning <- length(stop[stop == "List of explorable DLs empty after pruning"])
empty_expl_dl_before_pruning <- length(stop[stop == "List of explorable DLs empty before pruning"])

sim_summary <- matrix(round(c(mean(n_dlt),mean(n_enrolled), tot_perc_mtd,
                             average_enrolled_mtd, ###########
                             max_sample_size/length(stop)*100,
                             empty_expl_dl_after_pruning/length(stop)*100,
                             empty_expl_dl_before_pruning/length(stop)*100),2),7 ,1)

colnames(sim_summary) <- " "
rownames(sim_summary) <- c("Average number of total DLT", 
                    "Average number of total enrolled patients",
                    "Selection percent of MTD", "% of patients treated at MTD",
                    "% of times in which the maximum sample size was reached",
                    "% of times in which the list of explorable DLs was empty after pruning",
                    "% of times in which the list of explorable DLs was empty before pruning")   

percentage_mtd_selection = round(l2,2)
```

The outputs of the simulation study are the following:

-   `sim_summary`: summaries across the simulated trials, such as
    average number of total dose limiting toxicities, enrolled patients,
    etc;

-   `average_toxicity_rate`: average estimated toxicity rates for each
    dose level across simulations;

-   `average_dlt`: average dose limiting toxicities for each dose level
    across simulations;

-   `percentage_mtd_selection`: percentage of times a dose level is
    selected as MTD across simulations

``` r
sim_summary
```

                                                                                  
    Average number of total DLT                                               6.00
    Average number of total enrolled patients                                21.00
    Selection percent of MTD                                                 60.00
    % of patients treated at MTD                                             28.58
    % of times in which the maximum sample size was reached                   0.00
    % of times in which the list of explorable DLs was empty after pruning  100.00
    % of times in which the list of explorable DLs was empty before pruning   0.00

``` r
average_toxicity_rate
```

        A1  A2 A3 A4 A5
    B5  NA  NA NA NA NA
    B4  NA  NA NA NA NA
    B3  NA  NA NA NA NA
    B2 0.5  NA NA NA NA
    B1 0.0 0.5 NA NA NA

``` r
average_dlt
```

       A1 A2 A3 A4 A5
    B5 NA NA NA NA NA
    B4 NA NA NA NA NA
    B3 NA NA NA NA NA
    B2  3 NA NA NA NA
    B1  0  3 NA NA NA

``` r
average_enrolled
```

          A1    A2 A3 A4 A5
    B5  0.00  0.00  0  0  0
    B4  0.00  0.00  0  0  0
    B3  0.00  0.00  0  0  0
    B2 14.29  0.00  0  0  0
    B1 21.43 14.29  0  0  0

``` r
percentage_mtd_selection 
```

       A1 A2 A3 A4 A5
    B5  0  0  0  0  0
    B4  0  0  0  0  0
    B3 10 10  0  0  0
    B2 20  0  0  0  0
    B1 10 40 10  0  0
