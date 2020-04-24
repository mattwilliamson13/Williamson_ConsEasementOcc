library(rstan)
library(tidybayes)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

stan_data <- readRDS("D:/Data/IDMTDataArchive/IDMTstanData.rds")

stan.model <- list.files("./IDMT_CaseStudy/", pattern = "*.stan")

model <- stan_model(paste0("./IDMT_CaseStudy/", stan.model[1]))


n_iter <- 4700
n_wmp <- 3700
n_c <- 5

binom.fit <- sampling(model,
                data = stan_data,
                iter = n_iter, warmup=n_wmp, chains=n_c, seed = 8029,
                control=list(adapt_delta = 0.98, max_treedepth = 16))
saveRDS(binom.fit, "D:/Data/IDMTDataArchive/IDMT_binom.rds")

model <- stan_model(paste0("./IDMT_CaseStudy/", stan.model[2]))
psiCARdetCAR.fit <- sampling(model,
                      data = stan_data,
                      iter = n_iter, warmup=n_wmp, chains=n_c, seed = 8029,
                      control=list(adapt_delta = 0.98, max_treedepth = 16))
saveRDS(psiCARdetCAR.fit, "D:/Data/IDMTDataArchive/IDMT_psiCARdetCAR.rds")


model <- stan_model(paste0("./IDMT_CaseStudy/", stan.model[3]))
psiCARdetSTD.fit <- sampling(model,
                             data = stan_data,
                             iter = n_iter, warmup=n_wmp, chains=n_c, seed = 8029,
                             control=list(adapt_delta = 0.98, max_treedepth = 16))
saveRDS(psiCARdetSTD.fit, "D:/Data/IDMTDataArchive/IDMT_psiCARdetSTD.rds")

model <- stan_model(paste0("./IDMT_CaseStudy/", stan.model[4]))
psiSTDdetCAR.fit <- sampling(model,
                             data = stan_data,
                             iter = n_iter, warmup=n_wmp, chains=n_c, seed = 8029,
                             control=list(adapt_delta = 0.98, max_treedepth = 16))
saveRDS(psiSTDdetCAR.fit, "D:/Data/IDMTDataArchive/IDMT_psiSTDdetCAR.rds")

model <- stan_model(paste0("./IDMT_CaseStudy/", stan.model[5]))
psiSTDdetSTD.fit <- sampling(model,
                             data = stan_data,
                             iter = n_iter, warmup=n_wmp, chains=n_c, seed = 8029,
                             control=list(adapt_delta = 0.98, max_treedepth = 16))
saveRDS(psiSTDdetSTD.fit, "D:/Data/IDMTDataArchive/IDMT_psiSTDdetSTD.rds")
