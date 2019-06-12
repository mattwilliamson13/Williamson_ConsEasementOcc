#!/usr/bin/Rscript
library(lhs)
library(tidyverse)

# Get environment params from slurm ---------------------------------------
args <- commandArgs(trailingOnly = TRUE)
print(args)
run <- as.numeric(args[1])
print(paste0("run = ", run))
stan.idx <- as.numeric(args[2])
print(paste0("stan.idx = ", stan.idx))

inputn <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Set seed and simulaiton parms -------------------------------------------


set.seed(082980 * run)
h <- 300
latin.hyp <- maximinLHS(h, 6)
# Parameters that affect estimation ---------------------------------------

occ.rho.min <- 0.5
occ.rho.max <- 0.999
p.rho.min <- 0.5 #per M. Hooten anything below 0.5 is not likely to matter
p.rho.max <- 0.999
sub.occ.min <- 0.2 #this is the propbability that at subunit is available
sub.occ.max <- 0.8
occ.min <- 0.2
occ.max <- 0.8 
p.min <- 0.3
p.max <- 0.98
tau.min <- 0.1
tau.max <- 1


#sdev <- seq(1, 3, length.out = 1000)
#tau <- 1/(sdev^2)
#plot(tau ~ sdev)
# Map params to hypercube samples -----------------------------------------
params.set <- cbind(occRHO = latin.hyp[,1] * (occ.rho.max - occ.rho.min) + occ.rho.min,
                    pRHO = latin.hyp[,2] * (p.rho.max - p.rho.min) + p.rho.min,
                    subOcc = latin.hyp[,3] * (sub.occ.max - sub.occ.min) + sub.occ.min,
                    occ = latin.hyp[,4] * (occ.max - occ.min) + occ.min,
                    p = latin.hyp[,5] * (p.max - p.min) + p.min,
                    tau = latin.hyp[,6] * (tau.max - tau.min) + tau.min)


#mean.psi <- seq(from=0.2, to=0.8, by=0.2)

# Set design parameters ---------------------------------------------------
#m.cty <- 2 #Numver of cty-level params on occurrence
m.tct <- 4 #number of tract-level params on occurence
m.bg <- 3 #Number of block-group level params on reporting

# get a dummy geography ---------------------------------------------------
library(tigris)
library(sf)
options(tigris_use_cache = TRUE)
prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
st <- "IA"
geog.tct <- tracts(state=st) %>%
  as(., "sf") %>% st_transform(., crs= prj)
geog.bg <- block_groups(state=st) %>% 
  as(., "sf") %>% st_transform(., crs=prj)


# Generate spatially correlated observation values ------------------------
source("./Scripts/gen_sp_precMatx.R")
obs.prec.matx.tct <- gen_sp_precMatx(geog = geog.tct, rho=0.3, tau=1) #moderate level of autocorrelation in covariate observations
obs.cov.matx.tct <- solve(obs.prec.matx.tct[[2]]) #necessary for rmvn
obs.prec.matx.bg <- gen_sp_precMatx(geog=geog.bg, rho=0.3, tau=1)
obs.cov.matx.bg <- solve(obs.prec.matx.bg[[2]]) #necessary for rmvn

#X.cty <- matrix(rnorm((m.cty)*length(unique(geog.tct$COUNTYFP)), mean=0, sd=1.25), ncol=m.cty)
#X.cty <- apply(X.cty, 2, scale)
#X.tct <- scale(t(mvnfast::rmvn(m.tct, rep(0,NROW(obs.cov.matx.tct)), obs.cov.matx.tct, ncores=10)))

X.tct <- matrix(c(rep(1, nrow(geog.tct)), mvnfast::rmvn(m.tct-1, rep(0,NROW(obs.cov.matx.tct)),
                                                       obs.cov.matx.tct, ncores=10)), ncol=m.tct)
X.tct[,2:m.tct] <- apply(X.tct[,2:m.tct], 2, scale)

X.bg <- matrix(c(rep(1, nrow(geog.bg)), mvnfast::rmvn(m.bg-1, rep(0,NROW(obs.cov.matx.bg)),
                                                              obs.cov.matx.bg, ncores=10)), ncol=m.bg)
X.bg[,2:m.bg] <- apply(X.bg[,2:m.bg], 2, scale)


#X.bg <- scale(t(mvnfast::rmvn(m.bg, rep(0,NROW(obs.cov.matx.bg)), obs.cov.matx.bg, ncores=10)))
#tct.in.cty <- rethinking::coerce_index(geog.tct[order(geog.tct$GEOID),]$COUNTYFP) #because the covariance matrices reorder the geog file, have to reorder here to get proper indices

# Generate occurrence data ------------------------------------------------
#inputn <- 1

#b.cty <- runif(m.cty, -1,1)
b.tct <-c(unname(boot::logit(params.set[inputn, 4])), runif(m.tct-1, -3,3))

err.prec.matx.tct <- gen_sp_precMatx(geog = geog.tct, rho=params.set[inputn,1], tau=params.set[inputn,6])
err.cov.matx.tct <- solve(err.prec.matx.tct[[2]])
phi.occ <- mvnfast::rmvn(1, rep(0, nrow(err.cov.matx.tct)), err.cov.matx.tct, ncores = 10)
#mu.cty <- qlogis(mean.psi[inputn]) + X.cty %*% b.cty
logit.psi <- X.tct %*% b.tct + t(phi.occ)
z <- rbinom(nrow(geog.tct), size = 1, prob = boot::inv.logit(logit.psi))


# Generate detection data -------------------------------------------------
b.bg <- c(unname(boot::logit(params.set[inputn, 5])), runif(m.bg-1, -3,3))

survey_summary <- geog.bg %>% 
  group_by(., STATEFP, COUNTYFP, TRACTCE) %>% 
  summarise(count = n()) #NOTE: This orders the data; fine if CAR function also orders, but could screw things up

n_survey <- as.vector(survey_summary$count)
total_surveys <- nrow(geog.bg)
#X.bg <- matrix(c(rep(1, total_surveys), mvnfast::rmvn(m.bg-1, rep(0,NROW(obs.cov.matx.bg)),
#                                                              obs.cov.matx.bg, ncores=10)), ncol=m.bg)
#X.bg[,2:3] <- apply(X.bg[,2:3], 2, scale)

err.prec.matx.bg <- gen_sp_precMatx(geog = geog.bg, rho=params.set[inputn,2], tau=params.set[inputn,6])
err.cov.matx.bg <- solve(err.prec.matx.bg[[2]])
phi.p <- mvnfast::rmvn(1, rep(0, nrow(err.cov.matx.bg)), err.cov.matx.bg, ncores = 10)
logit.p <- X.bg %*% b.bg + t(phi.p)
a <- unlist(lapply(n_survey, function(x) rbinom(x,1, params.set[inputn, 3])))
p <- boot::inv.logit(logit.p)

survey.df <- tibble(site = rep(1:nrow(geog.tct), n_survey),
                    siteID = rep(geog.tct$GEOID, n_survey)) %>%
  mutate(y = rbinom(n = total_surveys, size = 1, prob = z[site] * a * p))
# get start and end indices to extract slices of y for each site
start_idx <- rep(0, nrow(geog.tct))
end_idx <- rep(0, nrow(geog.tct))
for (i in 1:nrow(geog.tct)) {
  if (n_survey[i] > 0) {
    site_indices <- which(survey.df$site == i)
    start_idx[i] <- site_indices[1]
    end_idx[i] <- site_indices[n_survey[i]]
  }
}

any_seen <- rep(0, nrow(geog.tct))
for (i in 1:nrow(geog.tct)) {
  if (n_survey[i] > 0) {
    any_seen[i] <- max(survey.df$y[start_idx[i]:end_idx[i]])
  }
}

# Prepare for stan --------------------------------------------------------
####MAKE SURE SEED IS WORKING TO ENSURE DATA IS THE SAME ACROSS RUNS!!
stan_data_mbj <- list(
  n_site = nrow(geog.tct),
  m_psi = m.tct,
  X_tct = X.tct,
  total_surveys = nrow(geog.bg),
  m_p = m.bg,
  X_bg = X.bg, 
  site = survey.df$site,
  y = survey.df$y,
  start_idx = start_idx,
  end_idx = end_idx,
  any_seen = any_seen,
  n_survey = n_survey,
  W_tct = err.prec.matx.tct[[1]],
  W_n_tct = sum(err.prec.matx.tct[[1]]/2),
  W_bg = err.prec.matx.bg[[1]],
  W_n_bg = sum(err.prec.matx.bg[[1]]/2)
)



# Run Stan ----------------------------------------------------------------
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

stan.model <- list.files("./Scripts/", pattern = "*.stan")

model <- stan_model(paste0("./Scripts/", stan.model[stan.idx]))

n_iter <- 3700
n_wmp <- 3200
n_c <- 4

start.stan <- Sys.time()
fit <- sampling(model,
                data = stan_data_mbj,
                iter = n_iter, warmup=n_wmp, chains=n_c, seed = 8029,
                control=list(adapt_delta = 0.98, max_treedepth = 16))

end.stan <- Sys.time()
stan.time <- round(end.stan - start.stan, digits = 2)

source("./Scripts/stan_utility.R")
n.divs <- count_divs(fit)
n.trdpth <- count_td(fit)
n.energy <- count_energy(fit)
n.rhat <- count_rhat(fit)

diagnostics.df <- data.frame(mod = stan.model[stan.idx],
                             runID = run,
                             occRHO = params.set[inputn, 1],
                             pRHO = params.set[inputn, 2],
                             avail = params.set[inputn, 3],
                             occ = params.set[inputn, 4],
                             p = params.set[inputn, 5],
                             tau = params.set[inputn, 6],
                             smplrtime = stan.time,
                             nDivs = n.divs,
                             TD = n.trdpth,
                             bfmi = n.energy,
                             rHats = n.rhat)
# Extract samples from model ----------------------------------------------
library(tidybayes)

fit.rec <- fit %>%  recover_types(stan_data_mbj)
reg.draws <- gather_draws(fit.rec, `beta_.*`, regex=TRUE) %>%  
  to_broom_names() 
regexp <- "[[:digit:]]+"

occ.bias <- reg.draws %>% 
  filter(., grepl("*psi.*", term)) %>% 
  mutate(truth = b.tct[as.numeric(str_extract(term, regexp))],
         relBias = (estimate - truth)/abs(truth))

p.bias <- reg.draws %>% 
  filter(., !grepl("*psi.*", term)) %>% 
  mutate(truth = b.bg[as.numeric(str_extract(term, regexp))],
         relBias = (estimate - truth)/abs(truth))
phi.draws <- NA
#phi.draws <- gather_draws(fit.rec, `phi_.*`, regex=TRUE) %>%  
#  to_broom_names()
phi.occ.bias <- NA
#phi.occ.bias <- phi.draws %>% 
#  filter(., grepl("*occ.*", term)) %>% 
#  mutate(truth = phi.occ[as.numeric(str_extract(term, regexp))],
#         relBias = (estimate - truth)/abs(truth))
phi.det.bias
#phi.det.bias <- phi.draws %>% 
#  filter(., grepl("*det.*", term)) %>% 
#  mutate(truth = phi.p[as.numeric(str_extract(term, regexp))],
#         relBias = (estimate - truth)/abs(truth))

# Save output -------------------------------------------------------------

outputs <- list(
  diagnostics = diagnostics.df,
  post_draws = reg.draws,
  occ_bias = occ.bias,
  p_bias = p.bias,
  phi_occ_bias = phi.occ.bias,
  phi_det_bias = phi.det.bias
)

saveRDS(outputs, file=paste0(stan.model[stan.idx],"LHS", inputn, "run", run, ".rds"), compress = "xz")


