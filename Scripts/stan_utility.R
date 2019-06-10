##from Michael Betancourt's Robust Statistical Workflow https://github.com/betanalpha/knitr_case_studies/blob/master/qr_regression/stan_utility.R
##modified to return actual counts rather than returning msgs

count_divs <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
}

# Check transitions that ended prematurely due to maximum tree depth limit
count_td <- function(fit, max_depth = 10) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
}

# Checks the energy Bayesian fraction of missing information (E-BFMI)
count_energy <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  n_energies <- 0
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      n_energies <- n_energies + 1
    } else {
      n_energies <- n_energies
    }
  }
return(n_energies)}

# Checks the effective sample size per iteration
count_n_eff <- function(fit) {
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  iter <- dim(extract(fit)[[1]])[[1]]
  
  n_warning <- 0
  for (n in 1:N) {
    ratio <- fit_summary[,5][n] / iter
    if (ratio < 0.001) {
      n_warning <- n_warning + 1
    } else {
      n_warning <- n_warning
    }
  }
return (n_warning)}

# Checks the potential scale reduction factors
count_rhat <- function(fit) {
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  n_warning <- 0
  for (n in 1:N) {
    rhat <- fit_summary[,6][n]
    if (rhat > 1.15 || is.infinite(rhat) || is.nan(rhat)) {
      n_warning <- n_warning + 1
    } else {
      n_warning <- n_warning
    }
  }
return(n_warning)}

