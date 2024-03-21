
sapply( list.files("R", full.names=TRUE), source )
N = 30
n_cohort = 3
toxicity_matrix = matrix(rbind(c(0.10,	0.12,	0.18,	0.33,	0.40),
                               c(0.10,	0.12,	0.15,	0.25,	0.27),
                               c(   0,  0.10,	0.10,	0.10,	0.20)),
                         3, 5, byrow = FALSE )

n_max = 6

trial_results <- ConductTrial(N = N,
                              n_cohort = n_cohort,
                              toxicity_matrix = toxicity_matrix,
                              n_max = n_max,
                              select_lambda = TRUE,
                              lambda_e = 0.17, lambda_d = 0.25, lambda_r = 0.33)

# Data accrued during the trial stage by stage
trial_results$trial_data

# Overall data accrued during the trial
trial_results$trial_data_summary

# Reason why the trial stopped
trial_results$stop
