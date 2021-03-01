# Code for "Improving inferences about private land conservation by accounting for incomplete reporting"

This repository contains the code for the simulations and case study contained in Williamson et. al 2021 "Improving inferences about private land conservation by accounting for incomplete reporting". Note that the shell scripts are HPC dependent and may require some adjustment to work with your local setup.

## Simulation Study
### Shell scripts
`simulate_shell.sh`and `simulate_shell_noCAR.sh` are shell scripts for running multiple simulations in parallel on a High Performance Computing Cluster and may need to be modified to fit your local context. The scripts require 2 arguments - the first is the simulation run (in our case a number between 1 and 3000), the second is an identifier that indicates which model to use (i.e., naive logistic regression, occupancy with conditional auto-regressive (CAR) term on both occurrence and detection, occupancy with CAR on occurrence only, occupancy with CAR on detection only, and occupancy with no CAR term).

## R scripts
`~./simrun.R` and `~./simrun_noCAR.R` simulate data according to the parameters passed from the shell scripts and then fit the appropriate Stan model to the subsequent data. `gen_spec_precMatx.R` is a utility function for generating the precision matrix necessary for imposing the desired level of spatial autocorrelation. `stan_utility.R` is a function created by [Michael Betancourt](https://github.com/betanalpha/knitr_case_studies/tree/master/principled_bayesian_workflow) for extracting diagnostics from a Stan run. 

