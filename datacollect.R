library(tidyverse)
library(magrittr)
library(tidybayes)
library(scales)
cluster.output.folder <- "D:/OccRuns/"

model.names <- c("binom", "psiCARdetCAR", "psiCARdetSTD" ,"psiSTDdetCAR", "psiSTDdetSTD")
task.id <- as.character(seq(from = 1, to = 300, by = 1))

get_file_list <- function(loc, nruns){
  pattern <- paste0("run", nruns, "\\b")
  model.list <- list.files(loc, pattern) 
}
n.runs <- 10
model.list <- map(1:n.runs, function(x) get_file_list(loc = cluster.output.folder, nruns = x))
model.tally <- map(seq_along(model.list), function (y) 
  map(model.names, function(x) sum(str_count(model.list[[y]], pattern = x))))



pat <- "(\\d)+"
successful.runs <- map(seq_along(model.list), function (y) 
  map(model.names, function(x) str_subset(model.list[[y]], pattern = x, 
                                          negate = FALSE) %>% 
        str_extract(., pattern = pat)))
missing <- map(seq_along(successful.runs), 
               function (x) map(seq_along(successful.runs[[x]]), function (y) 
                 subset(task.id, !(task.id %in% successful.runs[[x]][[y]]))))


# load all saved files ----------------------------------------------------
file.list.master <- unlist(model.list)

#file.list.sample <- unlist(map(model.names, function(x) str_subset(file.list.master, 
#                                                                   pattern = x, 
#                                                                   negate = FALSE) %>% 
#                                 sample(., size = 150)))


read_model_runs <- function(folder, filename, iter){
  input.file <- readRDS(paste0(folder, filename))
  occ.frame <- magrittr::extract2(input.file, "occ_bias")
  det.frame <- magrittr::extract2(input.file, "p_bias")
  diag.frame <- magrittr::extract2(input.file, "diagnostics")
  df.occ <- occ.frame %>% add_column(!!! diag.frame) %>% mutate(iterID = iter)
  df.det <- det.frame %>% add_column(!!! diag.frame) %>% mutate(iterID = iter)
  return(list(df.occ, df.det))
}



start.read <- Sys.time()  
psi.det.bias.list <- map(seq_along(file.list.master), function (x) 
  read_model_runs(folder = cluster.output.folder, filename = file.list.master[x], iter = x))
end.read <- Sys.time()

saveRDS(psi.det.bias.list, here::here("Outputs", "bias_list_070519.rds"))

###UPDATE THESE AFTER YOU READ IN


psi.reg.bias <-  map(seq_along(psi.det.bias.list), 
                     function(x) magrittr::extract2(psi.det.bias.list[[x]],1)) %>% 
  bind_rows() %>% 
  filter(term != "beta_psi[1]")

mean.p.bias <-  map(seq_along(psi.det.bias.list), 
                    function(x) magrittr::extract2(psi.det.bias.list[[x]],2)) %>% 
  bind_rows() %>% 
  filter(term == "beta_p[1]")

p.reg.bias <-  map(seq_along(psi.det.bias.list), 
                   function(x) magrittr::extract2(psi.det.bias.list[[x]],2)) %>% 
  bind_rows() %>% 
  filter(term != "beta_p[1]")

# Plot functions -------------------------------------------------------------------
mean.psi.bias <- map(seq_along(psi.det.bias.list), 
                     function(x) magrittr::extract2(psi.det.bias.list[[x]],1)) %>% 
  bind_rows() %>% 
  filter(term == "beta_psi[1]")



# Create custom log-style y axis transformer (0,1,3,10,...)
custom_log_y_trans <- function()
  trans_new("custom_log_y",
            transform = function (x) ( sign(x)*log10(abs(x)+1) ),
            inverse = function (y) ( sign(y)*( 10^(abs(y))-1) ),
            domain = c(-Inf,Inf))

inverse_logit_trans <- function()
  trans_new("inverse_logit",
            transform = function (x) boot::inv.logit(x) ,
            inverse = function(y) identity(y),
            domain = c(0,1))

# Custom log y breaker (0,1,3,10,...)
custom_y_breaks <- function(x)
{ 
  range <- max(abs(x), na.rm=TRUE)
  
  return (sort( c(0,
                  sapply(0:log10(range), function(z) (10^z) ),
                  sapply(0:log10(range/5), function(z) (5*10^z) ),
                  sapply(0:log10(range), function(z) (-10^z) ),
                  sapply(0:log10(range/5), function(z) (-5*10^z) )
  )))
}



gen_summary_plots <- function(base.data.frame, mod.name, plot.group){
  ##subset data and round based on the plotting variable of interest
  data.subset <- base.data.frame %>% 
    #filter(., rHats == 0) %>%
    filter(., mod == mod.name)  %>% 
    mutate(rounded = round(!!rlang::ensym(plot.group)/0.05) *0.05)
  
  ##get the median of the posterior of each individual run
  post.med.fit <- data.subset %>% 
    group_by(iterID, truth,p, avail, tau, occ) %>% 
    summarise(medbias = median(relBias))
  
  ##get the median for all models at a given rounded value for plotting
  post.med.mod <- data.subset %>%
    group_by(rounded) %>% 
    summarise(medbias = median(relBias))
  
  ##get the median for a rounded value within a run
  post.med.run <- data.subset %>%
    group_by(runID, rounded) %>% 
    summarise(medbias = median(relBias))
  
  ggplot(data = post.med.fit, mapping = aes(x=!!rlang::ensym(plot.group), y = medbias)) + 
    geom_point(color = "gray80", alpha=0.2) +
    geom_line(data = post.med.run, mapping = aes(x = rounded, y = medbias, group=runID), color = "gray60") +
    geom_line(data = post.med.mod, mapping = aes(x = rounded, y = medbias), color = "blue") +
    geom_point(data = subset(post.med.fit, p <= quantile(p, 0.1) & avail <= quantile(p, 0.1)), 
               mapping = aes(x=!!rlang::ensym(plot.group), y = medbias), color = "red") +
    scale_y_continuous(trans = 'custom_log_y', breaks=custom_y_breaks(post.med.fit$medbias)) +
    theme_minimal() +
    ylab("Relative bias")
  
}




mean.occ.plot.list.OCC <- map(seq_along(model.names), function(x) gen_summary_plots(base.data.frame = mean.psi.bias, 
                                                                                    mod.name = paste0(model.names[x],".stan"), plot.group = "occ"))

mean.occ.plot.list.DET <- map(seq_along(model.names), function(x) gen_summary_plots(base.data.frame = mean.psi.bias, 
                                                                                    mod.name = paste0(model.names[x],".stan"), plot.group = "p"))

mean.occ.plot.list.AVAIL <- map(seq_along(model.names), function(x) gen_summary_plots(base.data.frame = mean.psi.bias, 
                                                                                      mod.name = paste0(model.names[x],".stan"), plot.group = "avail"))

mean.occ.plot.list.TRUTH <- map(seq_along(model.names), function(x) gen_summary_plots(base.data.frame = mean.psi.bias, 
                                                                                      mod.name = paste0(model.names[x],".stan"), plot.group = "truth"))

saveRDS(mean.occ.plot.list.OCC,"D:/Data/occ_intercept_plotlist_occ.rds")
saveRDS(mean.occ.plot.list.DET,"D:/Data/occ_intercept_plotlist_det.rds")
saveRDS(mean.occ.plot.list.AVAIL,"D:/Data/occ_intercept_plotlist_avail.rds")
saveRDS(mean.occ.plot.list.TRUTH,"D:/Data/occ_intercept_plotlist_truth.rds")



occ.reg.plot.list.OCC <- map(seq_along(model.names), function(x) gen_summary_plots(base.data.frame = psi.reg.bias, 
                                                                                   mod.name = paste0(model.names[x],".stan"), plot.group = "occ"))


occ.reg.plot.list.DET <- map(seq_along(model.names), function(x) gen_summary_plots(base.data.frame = psi.reg.bias, 
                                                                                   mod.name = paste0(model.names[x],".stan"), plot.group = "p"))

occ.reg.plot.list.AVAIL <- map(seq_along(model.names), function(x) gen_summary_plots(base.data.frame = psi.reg.bias, 
                                                                                     mod.name = paste0(model.names[x],".stan"), plot.group = "avail"))

occ.reg.plot.list.TRUTH <- map(seq_along(model.names), function(x) gen_summary_plots(base.data.frame = psi.reg.bias, 
                                                                                     mod.name = paste0(model.names[x],".stan"), plot.group = "truth"))

saveRDS(occ.reg.plot.list.OCC,"D:/Data/occ_reg_plotlist_occ.rds")
saveRDS(occ.reg.plot.list.DET,"D:/Data/occ_reg_plotlist_det.rds")
saveRDS(occ.reg.plot.list.AVAIL,"D:/Data/occ_reg_plotlist_avail.rds")
saveRDS(occ.reg.plot.list.TRUTH,"D:/Data/occ_reg_plotlist_truth.rds")

det.models <- model.names[-1]
mean.p.plot.list.OCC <- map(seq_along(det.models), function(x) gen_summary_plots(base.data.frame = mean.p.bias, 
                                                                                 mod.name = paste0(det.models[x],".stan"), plot.group = "occ"))


mean.p.plot.list.DET <- map(seq_along(det.models), function(x) gen_summary_plots(base.data.frame = mean.p.bias, 
                                                                                 mod.name = paste0(det.models[x],".stan"), plot.group = "p"))

mean.p.plot.list.AVAIL <- map(seq_along(det.models), function(x) gen_summary_plots(base.data.frame = mean.p.bias, 
                                                                                   mod.name = paste0(det.models[x],".stan"), plot.group = "avail"))

mean.p.plot.list.TRUTH <- map(seq_along(det.models), function(x) gen_summary_plots(base.data.frame = mean.p.bias, 
                                                                                   mod.name = paste0(det.models[x],".stan"), plot.group = "truth"))
saveRDS(mean.p.plot.list.OCC,"D:/Data/p_intercept_plotlist_occ.rds")
saveRDS(mean.p.plot.list.DET,"D:/Data/p_intercept_plotlist_det.rds")
saveRDS(mean.p.plot.list.AVAIL,"D:/Data/p_intercept_plotlist_avail.rds")
saveRDS(mean.p.plot.list.TRUTH,"D:/Data/p_intercept_plotlist_truth.rds")

p.reg.plot.list.OCC <- map(seq_along(det.models), function(x) gen_summary_plots(base.data.frame = p.reg.bias, 
                                                                                mod.name = paste0(det.models[x],".stan"), plot.group = "occ"))


p.reg.plot.list.DET <- map(seq_along(det.models), function(x) gen_summary_plots(base.data.frame = p.reg.bias, 
                                                                                mod.name = paste0(det.models[x],".stan"), plot.group = "p"))

p.reg.plot.list.AVAIL <- map(seq_along(det.models), function(x) gen_summary_plots(base.data.frame = p.reg.bias, 
                                                                                  mod.name = paste0(det.models[x],".stan"), plot.group = "avail"))

p.reg.plot.list.TRUTH <- map(seq_along(det.models), function(x) gen_summary_plots(base.data.frame = p.reg.bias, 
                                                                                  mod.name = paste0(det.models[x],".stan"), plot.group = "truth"))

saveRDS(p.reg.plot.list.OCC,"D:/Data/p_reg_plotlist_occ.rds")
saveRDS(p.reg.plot.list.DET,"D:/Data/p_reg_plotlist_det.rds")
saveRDS(p.reg.plot.list.AVAIL,"D:/Data/p_reg_plotlist_avail.rds")
saveRDS(p.reg.plot.list.TRUTH,"D:/Data/p_reg_plotlist_truth.rds")



pdf("./Outputs/meanoccplots.pdf")
mean.occ.plot.list
dev.off()


p <- ggplot(data = occ.intercept.bias.run, mapping = aes(x = occRd, y = medbias, group = factor(runID)))
beta0.psi.plot <- p + geom_point(data = subset(occ.intercept.bias.iter, 
                                               mod = "binom.stan" & medbias >= quantile(medbias, 0.1) & medbias <= quantile(medbias, 0.9)), 
                                 mapping = aes(x = occ, y = medbias), color = "gray80", alpha = 0.2, inherit.aes = FALSE) +
  geom_point(data = subset(occ.intercept.bias.iter,  mod = "binom.stan" & p < quantile(p, 0.1) & avail < quantile(avail, 0.1)), 
             mapping = aes(x = occ, y = medbias), color = "red", inherit.aes = FALSE) +
  geom_line(color = "gray40") +
  geom_line(data = occ.intercept.bias.mod, mapping = aes(x = occRd, y = medbias), color = "blue", inherit.aes = FALSE) +
  scale_y_continuous(trans = 'custom_log_y', breaks=custom_y_breaks(occ.intercept.bias.iter$medbias)) 


+
  geom_point(data = subset(occ.intercept.bias.iter, p < quantile(p, 0.1) & avail < quantile(avail, 0.1)), 
             mapping = aes(x = occ, y = medbias), color = "red", inherit.aes = FALSE) +
  geom_line(color = "gray40") +
  geom_line(data = occ.intercept.bias.mod, mapping = aes(x = occRd, y = medbias), color = "blue", inherit.aes = FALSE) +
  scale_y_continuous(trans = 'asinh') +
  theme_bw() +
  xlab("Mean occupancy probability") +
  ylab("Relative Bias") +
  facet_wrap(. ~ mod, scales = "free_y")

ggsave(here::here("Outputs", "beta0_psi_plot.png"), plot = beta0.psi.plot)

occ.reg.bias.run <- psi.reg.bias %>% 
  #filter(., rHats == 0) %>% 
  mutate(occRd = plyr::round_any(occ, 0.05)) %>% 
  group_by(mod, occRd, runID) %>% 
  summarise(mnbias = mean(relBias),
            medbias = median(relBias))

occ.reg.bias.mod <-  psi.reg.bias %>% 
  #filter(., rHats == 0) %>% 
  mutate(occRd = plyr::round_any(occ, 0.05)) %>% 
  group_by(mod, occRd) %>% 
  summarise(mnbias = mean(relBias),
            medbias = median(relBias))

occ.reg.bias.iter <- psi.reg.bias %>% 
  #filter(., rHats == 0) %>% 
  group_by(iterID,mod, p, avail, occ) %>% 
  summarise(mnbias = mean(relBias),
            medbias = median(relBias))


p2 <- ggplot(data = occ.reg.bias.run, mapping = aes(x = occRd, y = medbias, group = factor(runID)))
betaAll.psi.plot <- p2 + geom_point(data = occ.reg.bias.iter, mapping = aes(x = occ, y = medbias), color = "gray80", alpha = 0.2, inherit.aes = FALSE) +
  geom_point(data = subset(occ.reg.bias.iter, p < quantile(p, 0.1) & avail < quantile(avail, 0.1)), 
             mapping = aes(x = occ, y = medbias), color = "red", inherit.aes = FALSE) +
  geom_line(color = "gray40") +
  geom_line(data = occ.reg.bias.mod, mapping = aes(x = occRd, y = medbias), color = "blue", inherit.aes = FALSE) +
  theme_bw() +
  xlab("Mean occupancy probability") +
  ylab("Relative Bias") +
  facet_wrap(. ~ mod, scales = "free_y")
ggsave(here::here("Outputs", "betaAll_psi_plot.png"), plot = betaAll.psi.plot)




p.intercept.bias.run <- mean.p.bias %>% 
  #filter(., rHats == 0) %>% 
  mutate(detRd = plyr::round_any(p, 0.05)) %>% 
  group_by(mod, detRd, runID) %>% 
  summarise(mnbias = mean(relBias),
            medbias = median(relBias))

p.intercept.bias.mod <-  mean.p.bias %>% 
  #filter(., rHats == 0) %>% 
  mutate(detRd = plyr::round_any(p, 0.05)) %>% 
  group_by(mod, detRd) %>% 
  summarise(mnbias = mean(relBias),
            medbias = median(relBias))

p.intercept.bias.iter <- mean.p.bias %>% 
  #filter(., rHats == 0) %>% 
  group_by(iterID,mod, p, avail, tau, occ) %>% 
  summarise(mnbias = mean(relBias),
            medbias = median(relBias))


p3 <- ggplot(data = p.intercept.bias.run, mapping = aes(x = detRd, y = medbias, group = factor(runID)))
beta0.p.plot <- p3 + geom_point(data = subset(p.intercept.bias.iter, medbias >= quantile(medbias, 0.1) & medbias <= quantile(medbias, 0.9)), 
                                mapping = aes(x = p, y = medbias), color = "gray80", alpha = 0.2, inherit.aes = FALSE) +
  geom_point(data = subset(p.intercept.bias.iter, occ < quantile(p, 0.1) & avail < quantile(avail, 0.1)), 
             mapping = aes(x = p, y = medbias), color = "red", inherit.aes = FALSE) +
  geom_line(color = "gray40") +
  geom_line(data = p.intercept.bias.mod, mapping = aes(x = detRd, y = medbias), color = "blue", inherit.aes = FALSE) +
  theme_bw() +
  xlab("Mean detection probability") +
  ylab("Relative Bias") +
  facet_wrap(. ~ mod, scales = "free_y")


p.reg.bias.run <- p.reg.bias %>% 
  #filter(., rHats == 0) %>% 
  mutate(detRd = plyr::round_any(p, 0.05)) %>% 
  group_by(mod, detRd, runID) %>% 
  summarise(mnbias = mean(relBias),
            medbias = median(relBias))

p.reg.bias.mod <-  p.reg.bias %>% 
  #filter(., rHats == 0) %>% 
  mutate(detRd = plyr::round_any(p, 0.05)) %>% 
  group_by(mod, detRd) %>% 
  summarise(mnbias = mean(relBias),
            medbias = median(relBias))

p.reg.bias.iter <- p.reg.bias %>% 
  #filter(., rHats == 0) %>% 
  group_by(iterID,mod, p, avail, occ) %>% 
  summarise(mnbias = mean(relBias),
            medbias = median(relBias))

p4 <- ggplot(data = p.reg.bias.run, mapping = aes(x = detRd, y = medbias, group = factor(runID)))
betaAll.p.plot <- p4 + geom_point(data = p.reg.bias.iter, mapping = aes(x = p, y = medbias), color = "gray80", alpha = 0.2, inherit.aes = FALSE) +
  geom_point(data = subset(p.reg.bias.iter, occ < quantile(occ, 0.1) & avail < quantile(occ, 0.1)), 
             mapping = aes(x = p, y = medbias), color = "red", inherit.aes = FALSE) +
  geom_line(color = "gray40") +
  geom_line(data = p.reg.bias.mod, mapping = aes(x = detRd, y = medbias), color = "blue", inherit.aes = FALSE) +
  theme_bw() +
  xlab("Mean detection probability") +
  ylab("Relative Bias") +
  facet_wrap(. ~ mod, scales = "free_y")

# Code Graveyard ----------------------------------------------------------


d <- mean.psi.bias %>% 
  filter(., rHats == 0) %>% 
  mutate(occRd = round(occ, digits =1)) %>% 
  group_by(mod, occRd) %>% 
  median_qi(relBias, .width = c(.50, .80, .95))  

d$.lower <- if_else(d$.lower < -100, -100, d$.lower)

weird <- scales::trans_new("signed_log",
                           transform=function(x) sign(x)*log10(abs(x)),
                           inverse=function(x) sign(x)*10^(abs(x)))
###NEED TO FIGURE OUT HOW TO GET AXES BETTER PLOTTED
p <- ggplot(data = d)
p + geom_interval(aes(x=occRd, y = relBias))+
  #scale_y_continuous(trans = scales::atanh_trans()) +
  scale_color_brewer() +
  geom_line(data = mean.bias.occ, mapping = aes(x = occRd, y = mnbias)) +
  facet_wrap(.~mod, scales = "free_y")

geom_smooth(se=FALSE) + 
  geom_interval(aes(ymin = .lower, ymax = .upper, color = fct_rev(ordered(.width)))) +
  scale_color_brewer() +
  geom_smooth(data = mean.psi.bias, mapping(x = occRd, y = mnbias)) +
  facet_wrap(.~mod, scales = "free_y")

occ + geom_line() + facet_wrap(.~mod, scales = 'free_y')
occ + geom_jitter(data =  mean.psi.bias %>% 
                    group_by(mod,runID,occRHO,pRHO,avail,occ,p,tau) %>% sample_n(50), 
                  mapping = aes(x = occ, y = relBias), col = "blue", alpha = 0.3) +
  geom_line(data = psi.bias.summary, mapping = aes(x = occ, y = bias)) +
  geom_ribbon(mapping = aes(ymin = lo, ymax = hi), fill = "red") + 
  facet_wrap(.~mod, scales = 'free_y')

det <- ggplot(data = psi.bias.summary, mapping = aes(x = p, y = bias))

det + geom_jitter(data =  mean.psi.bias %>% 
                    group_by(mod,runID,occRHO,pRHO,avail,occ,p,tau) %>% sample_n(50), 
                  mapping = aes(x = p, y = relBias, color = avail), alpha = 0.3) +
  scale_fill_viridis_c() +
  geom_line(data = psi.bias.summary, mapping = aes(x = p, y = bias)) +
  geom_ribbon(mapping = aes(ymin = lo, ymax = hi)) + 
  facet_wrap(.~mod, scales = 'free_y')


avail <- ggplot(data = psi.bias.summary, mapping = aes(x = avail, y = bias))

avail + geom_jitter(data =  mean.psi.bias %>% 
                      group_by(mod,runID,occRHO,pRHO,avail,occ,p,tau) %>% sample_n(50), 
                    mapping = aes(x = avail, y = relBias, color = p), alpha = 0.3) +
  scale_fill_viridis_c() +
  geom_line(data = psi.bias.summary, mapping = aes(x = avail, y = bias)) +
  geom_ribbon(mapping = aes(ymin = lo, ymax = hi)) + 
  facet_wrap(.~mod, scales = 'free_y')


psi.bias.heatmap <- mean.psi.bias %>% 
  mutate(occGroup = round(occ, digits = 1),
         detGroup = round(p, digits = 1),
         availGroup = round(avail, digits = 1)) %>% 
  group_by(mod,runID,occGroup, detGroup, availGroup) %>% 
  summarise_at("relBias", mean)

p <- ggplot(data = psi.bias.heatmap, mapping = aes(x = occGroup, y =availGroup))
p + geom_tile(aes(fill = log(abs(relBias))), interpolate = TRUE) +
  scale_fill_viridis_c(option = "magma")+ 
  facet_wrap(.~factor(mod))




# Regression Params -------------------------------------------------------






psi.bias.heatmap.occ.p <- psi.bias.heatmap %>% 
  group_by(occGroup, detGroup, mod) %>% 
  summarise(bias = mean(relBias),
            lo = quantile(relBias, 0.1),
            hi = quantile(relBias, 0.9))


