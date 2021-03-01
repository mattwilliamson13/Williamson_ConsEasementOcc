library(tidyverse)
library(magrittr)
library(tidybayes)
library(scales)
library(patchwork)
cluster.output.folder <- "D:/OccRuns/" #location where you transfered all successful model runs

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

#this is used to make sure that the correct number of simulations actually finished
successful.runs <- map(seq_along(model.list), function (y) 
  map(model.names, function(x) str_subset(model.list[[y]], pattern = x, 
                                          negate = FALSE) %>% 
        str_extract(., pattern = pat)))
missing <- map(seq_along(successful.runs), 
               function (x) map(seq_along(successful.runs[[x]]), function (y) 
                 subset(task.id, !(task.id %in% successful.runs[[x]][[y]]))))


# load all saved files ----------------------------------------------------
file.list.master <- unlist(model.list)


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

###UPDATE THESE AFTER YOU READ IN (prevent memory overruns)

psi.det.bias.list <- readRDS(here::here("DataArchive", "bias_list_070519.rds"))
mean.psi.bias <- map(seq_along(psi.det.bias.list), 
                     function(x) magrittr::extract2(psi.det.bias.list[[x]],1)) %>% 
  bind_rows() %>% 
  filter(grepl("beta_psi*", term)) 

psi.bias <- mean.psi.bias %>%
  mutate(param = if_else(term == "beta_psi[1]", "Intercept", "Coefficient"))

mean.p.bias <-  map(seq_along(psi.det.bias.list), 
                    function(x) magrittr::extract2(psi.det.bias.list[[x]],2)) %>% 
  bind_rows() %>% 
  filter(term == "beta_p[1]")

p.reg.bias <-  map(seq_along(psi.det.bias.list), 
                   function(x) magrittr::extract2(psi.det.bias.list[[x]],2)) %>% 
  bind_rows() %>% 
  filter(term != "beta_p[1]")

# Plot functions -------------------------------------------------------------------



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
  pg <- plot.group
  data.subset.05 <- base.data.frame %>% 
    #filter(., rHats == 0) %>%
    filter(., mod == mod.name)  %>% 
    mutate(rounded = round(!!rlang::ensym(pg)/0.05) *0.05) %>% 
    group_by(rounded) %>% 
    median_qi(relBias, .width = c(0.5, 0.75, 0.9)) 
 
   data.subset.1 <-base.data.frame %>% 
    #filter(., rHats == 0) %>%
    filter(., mod == mod.name)  %>% 
    mutate(rounded = round(!!rlang::ensym(pg)/0.1) *0.1) %>% 
    group_by(rounded) %>% 
    median_qi(relBias, .width = c(0.5, 0.75, 0.9))
  
  p <- ggplot() + 
    geom_pointinterval(data = data.subset.1, mapping = aes(x = rounded, y = relBias)) +
    geom_hline(yintercept = 0, color = "gray60", linetype = 2, size = 1.1) +
    geom_line(data = data.subset.05, mapping = aes(x = rounded, 
                                                   y = relBias),color = "blue", size = 1.5) +
    ylab("Relative bias") +
    xlab(plot.group)+
    theme_bw()
  }

psi.bias.smple <- psi.bias %>% filter(param == "Coefficient")
psi.bias.smple.int <- psi.bias %>% filter(param == "Intercept")
plt.grps <- c("p","avail", "occRHO")
binom.plots <- map(seq_along(plt.grps), function (x)
                   gen_summary_plots(base.data.frame = psi.bias.smple, mod.name = "binom.stan", plot.group = as.character(plt.grps[x])))  

binom.plot.group <- binom.plots[[1]] + xlab("Reporting probability") + labs(title = "Naive Logistic Regression", subtitle = "Reporting") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  binom.plots[[2]] + xlab("Availability") + ylab("") + labs(subtitle = "Availability") + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) +
  binom.plots[[3]] + xlab(expression(rho[occurrence])) + ylab("") + labs(subtitle = "Spatial \nAutocorrelation") + theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) + plot_layout(ncol = 3)  
  
psi.CAR.det.CAR.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = psi.bias.smple, mod.name = "psiCARdetCAR.stan", plot.group = as.character(plt.grps[x])))  

psi.CAR.det.CAR.plot.group <- psi.CAR.det.CAR.plots[[1]] + xlab("Reporting probability") + labs(title = "Occupancy Model with CAR on both components \n(OCC-CAR1)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.CAR.det.CAR.plots[[2]] + xlab("Availability") + ylab("") + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.CAR.det.CAR.plots[[3]] + xlab(expression(rho[occurrence])) + ylab("") + theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) + plot_layout(ncol = 3)  


psi.CAR.det.STD.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = psi.bias.smple, mod.name = "psiCARdetSTD.stan", plot.group = as.character(plt.grps[x])))  

psi.CAR.det.STD.plot.group <- psi.CAR.det.STD.plots[[1]] + xlab("Reporting probability") + labs(title = "Occupancy Model with CAR on occurrence only \n(OCC-CAR2)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.CAR.det.STD.plots[[2]] + xlab("Availability") + ylab("")  + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.CAR.det.STD.plots[[3]] + xlab(expression(rho[occurrence])) + ylab("") + theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) + plot_layout(ncol = 3)  

psi.STD.det.CAR.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = psi.bias.smple, mod.name = "psiSTDdetCAR.stan", plot.group = as.character(plt.grps[x])))  

psi.STD.det.CAR.plot.group <- psi.STD.det.CAR.plots[[1]] + xlab("Reporting probability") + labs(title = "Occupancy Model with CAR on reporting only \n(OCC-CAR3)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.STD.det.CAR.plots[[2]] + xlab("Availability") + ylab("")  + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.STD.det.CAR.plots[[3]] + xlab(expression(rho[occurrence])) + ylab("") + theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) + plot_layout(ncol = 3)  


psi.STD.det.STD.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = psi.bias.smple, mod.name = "psiSTDdetSTD.stan", plot.group = as.character(plt.grps[x])))  

psi.STD.det.STD.plot.group <- psi.STD.det.STD.plots[[1]] + xlab("Reporting probability") + labs(title = "Occupancy Model with no CAR \n(OCC)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.STD.det.STD.plots[[2]] + xlab("Availability") + ylab("")  + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.STD.det.STD.plots[[3]] + xlab(expression(rho[occurrence])) + ylab("") +  theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) + plot_layout(ncol = 3)  


all.occ.reg.plots <-  binom.plots[[1]] + xlab("") + labs(subtitle = "Naive Logistic Regression (NLR)", title = "Reporting") + theme(panel.background = element_rect(fill = "lightgray"), 
                                                                                                                              plot.background = element_rect(fill= "lightgray"),
                                                                                                                              plot.title = element_text(size=14,face = "bold",hjust = 0.5)) +
  binom.plots[[2]] + xlab("") + ylab("") + labs(title = "Availability") + theme(panel.background = element_rect(fill = "lightcyan"), 
                                                                                plot.background = element_rect(fill= "lightcyan"),
                                                                                plot.title = element_text(hjust=0.5, face="bold", size=14)) +
  binom.plots[[3]] + xlab("") + ylab("") + labs(title = "Spatial \nAutocorrelation") + theme(panel.background = element_rect(fill = "cornsilk"), 
                                                                                             plot.background = element_rect(fill= "cornsilk"),
                                                                                             plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +   
  psi.CAR.det.CAR.plots[[1]] + xlab("") + labs(subtitle = "Occupancy Model with CAR on both components (OCC-CAR1)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.CAR.det.CAR.plots[[2]] + xlab("") + ylab("") + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.CAR.det.CAR.plots[[3]] + xlab("") + ylab("") + theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) +
  psi.CAR.det.STD.plots[[1]] + xlab("") + labs(subtitle = "Occupancy Model with CAR on occurrence only (OCC-CAR2)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.CAR.det.STD.plots[[2]] + xlab("") + ylab("")  + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.CAR.det.STD.plots[[3]] + xlab("") + ylab("") + theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) + 
  psi.STD.det.CAR.plots[[1]] + xlab("") + labs(subtitle = "Occupancy Model with CAR on reporting only (OCC-CAR3)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.STD.det.CAR.plots[[2]] + xlab("") + ylab("")  + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.STD.det.CAR.plots[[3]] + xlab("") + ylab("") + theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) +
  psi.STD.det.STD.plots[[1]] + xlab("Reporting probability") + labs(subtitle = "Occupancy Model with no CAR (OCC)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.STD.det.STD.plots[[2]] + xlab("Availability") + ylab("")  + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.STD.det.STD.plots[[3]] + xlab(expression(rho[occurrence])) + ylab("") +  theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) + plot_layout(ncol = 3)  

ggsave(here::here("Outputs", "occregressionplots.png"), all.occ.reg.plots, height = 9, width = 7.5, units = "in")


# occupancy intercept plots -----------------------------------------------
binom.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = psi.bias.smple.int, mod.name = "binom.stan", plot.group = as.character(plt.grps[x])))  

psi.CAR.det.CAR.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = psi.bias.smple.int, mod.name = "psiCARdetCAR.stan", plot.group = as.character(plt.grps[x])))  

psi.CAR.det.STD.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = psi.bias.smple.int, mod.name = "psiCARdetSTD.stan", plot.group = as.character(plt.grps[x])))  

psi.STD.det.CAR.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = psi.bias.smple.int, mod.name = "psiSTDdetCAR.stan", plot.group = as.character(plt.grps[x])))  

psi.STD.det.STD.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = psi.bias.smple.int, mod.name = "psiSTDdetSTD.stan", plot.group = as.character(plt.grps[x])))  

all.occ.int.plots <-  binom.plots[[1]] + xlab("") + labs(subtitle = "Naive Logistic Regression (NLR)", title = "Reporting") + theme(panel.background = element_rect(fill = "lightgray"), 
                                                                                                                                    plot.background = element_rect(fill= "lightgray"),
                                                                                                                                    plot.title = element_text(size=14,face = "bold",hjust = 0.5)) +
  binom.plots[[2]] + xlab("") + ylab("") + labs(title = "Availability") + theme(panel.background = element_rect(fill = "lightcyan"), 
                                                                                plot.background = element_rect(fill= "lightcyan"),
                                                                                plot.title = element_text(hjust=0.5, face="bold", size=14)) +
  binom.plots[[3]] + xlab("") + ylab("") + labs(title = "Spatial \nAutocorrelation") + theme(panel.background = element_rect(fill = "cornsilk"), 
                                                                                             plot.background = element_rect(fill= "cornsilk"),
                                                                                             plot.title = element_text(hjust = 0.5, face = "bold", size =14)) +   
  psi.CAR.det.CAR.plots[[1]] + xlab("") + labs(subtitle = "Occupancy Model with CAR on both components (OCC-CAR1)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.CAR.det.CAR.plots[[2]] + xlab("") + ylab("") + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.CAR.det.CAR.plots[[3]] + xlab("") + ylab("") + theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) +
  psi.CAR.det.STD.plots[[1]] + xlab("") + labs(subtitle = "Occupancy Model with CAR on occurrence only (OCC-CAR2)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.CAR.det.STD.plots[[2]] + xlab("") + ylab("")  + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.CAR.det.STD.plots[[3]] + xlab("") + ylab("") + theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) + 
  psi.STD.det.CAR.plots[[1]] + xlab("") + labs(subtitle = "Occupancy Model with CAR on reporting only (OCC-CAR3)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.STD.det.CAR.plots[[2]] + xlab("") + ylab("")  + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.STD.det.CAR.plots[[3]] + xlab("") + ylab("") + theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) +
  psi.STD.det.STD.plots[[1]] + xlab("Reporting probability") + labs(subtitle = "Occupancy Model with no CAR (OCC)") + theme(panel.background = element_rect(fill = "lightgray"), plot.background = element_rect(fill= "lightgray")) +
  psi.STD.det.STD.plots[[2]] + xlab("Availability") + ylab("")  + theme(panel.background = element_rect(fill = "lightcyan"), plot.background = element_rect(fill= "lightcyan")) + 
  psi.STD.det.STD.plots[[3]] + xlab(expression(rho[occurrence])) + ylab("") +  theme(panel.background = element_rect(fill = "cornsilk"), plot.background = element_rect(fill= "cornsilk")) + plot_layout(ncol = 3)  

ggsave(here::here("Outputs", "occintplots.png"), all.occ.int.plots, height = 9, width = 7.5, units = "in")


plt.grps <- c("occ", "p", "avail", "truth")

binom.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = mean.psi.bias, mod.name = "binom.stan", plot.group = as.character(plt.grps[x])))  

psi.CAR.det.CAR.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = mean.psi.bias, mod.name = "psiCARdetCAR.stan", plot.group = as.character(plt.grps[x])))  

psi.CAR.det.STD.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = mean.psi.bias, mod.name = "psiCARdetSTD.stan", plot.group = as.character(plt.grps[x])))  

psi.STD.det.CAR.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = mean.psi.bias, mod.name = "psiSTDdetCAR.stan", plot.group = as.character(plt.grps[x])))  

psi.STD.det.STD.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = mean.psi.bias, mod.name = "psiSTDdetSTD.stan", plot.group = as.character(plt.grps[x])))  

all.occ.int.plots <- binom.plots[[1]] + xlab("") + labs(subtitle = "Logistic (CAR on occurrence)") +
  binom.plots[[2]] + xlab("") + ylab("") + 
  binom.plots[[3]] + xlab("") + ylab("") + 
  binom.plots[[4]] + xlab("") + ylab("") +
  
  psi.CAR.det.CAR.plots[[1]] + xlab("") + labs(subtitle = "Occupancy (CAR on both components)") +
  psi.CAR.det.CAR.plots[[2]] + xlab("") + ylab("") + 
  psi.CAR.det.CAR.plots[[3]] + xlab("") + ylab("") + 
  psi.CAR.det.CAR.plots[[4]] + xlab("") + ylab("") +
  psi.CAR.det.STD.plots[[1]] + xlab("") + labs(subtitle = "Occupancy (CAR on occupancy only)") +
  psi.CAR.det.STD.plots[[2]] + xlab("") + ylab("") + 
  psi.CAR.det.STD.plots[[3]] + xlab("") + ylab("") + 
  psi.CAR.det.STD.plots[[4]] + xlab("") + ylab("") + 
  psi.STD.det.CAR.plots[[1]] + xlab("") + labs(subtitle = "Occupancy (CAR on detection only)") +
  psi.STD.det.CAR.plots[[2]] + xlab("") + ylab("") + 
  psi.STD.det.CAR.plots[[3]] + xlab("") + ylab("") + 
  psi.STD.det.CAR.plots[[4]] + xlab("") + ylab("") +
  psi.STD.det.STD.plots[[1]] + xlab("Occupancy probability") + labs(subtitle = "Occupancy (no CAR)") +
  psi.STD.det.STD.plots[[2]] + xlab("Detection probability") + ylab("") + 
  psi.STD.det.STD.plots[[3]] + xlab("Availability") + ylab("") + 
  psi.STD.det.STD.plots[[4]] + xlab("True value") + ylab("") + plot_layout(nrow = 5)

ggsave(here::here("Outputs", "occinterceptplots.png"), all.occ.int.plots, height = 8.25, width = 10.75, units = "in")


# Regression coeffs for detection -----------------------------------------
psi.CAR.det.CAR.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = p.reg.bias, mod.name = "psiCARdetCAR.stan", plot.group = as.character(plt.grps[x])))  

psi.CAR.det.STD.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = p.reg.bias, mod.name = "psiCARdetSTD.stan", plot.group = as.character(plt.grps[x])))  

psi.STD.det.CAR.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = p.reg.bias, mod.name = "psiSTDdetCAR.stan", plot.group = as.character(plt.grps[x])))  

psi.STD.det.STD.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = p.reg.bias, mod.name = "psiSTDdetSTD.stan", plot.group = as.character(plt.grps[x])))  

all.p.reg.plots <- psi.CAR.det.CAR.plots[[1]] + xlab("") + labs(subtitle = "Occupancy (CAR on both components)") +
  psi.CAR.det.CAR.plots[[2]] + xlab("") + ylab("") + 
  psi.CAR.det.CAR.plots[[3]] + xlab("") + ylab("") + 
  psi.CAR.det.CAR.plots[[4]] + xlab("") + ylab("") +
  psi.CAR.det.STD.plots[[1]] + xlab("") + labs(subtitle = "Occupancy (CAR on occupancy only)") +
  psi.CAR.det.STD.plots[[2]] + xlab("") + ylab("") + 
  psi.CAR.det.STD.plots[[3]] + xlab("") + ylab("") + 
  psi.CAR.det.STD.plots[[4]] + xlab("") + ylab("") + 
  psi.STD.det.CAR.plots[[1]] + xlab("") + labs(subtitle = "Occupancy (CAR on detection only)") +
  psi.STD.det.CAR.plots[[2]] + xlab("") + ylab("") + 
  psi.STD.det.CAR.plots[[3]] + xlab("") + ylab("") + 
  psi.STD.det.CAR.plots[[4]] + xlab("") + ylab("") +
  psi.STD.det.STD.plots[[1]] + xlab("Occupancy probability") + labs(subtitle = "Occupancy (no CAR)") +
  psi.STD.det.STD.plots[[2]] + xlab("Detection probability") + ylab("") + 
  psi.STD.det.STD.plots[[3]] + xlab("Availability") + ylab("") + 
  psi.STD.det.STD.plots[[4]] + xlab("True value") + ylab("") + plot_layout(nrow = 4)

ggsave(here::here("Outputs", "detregplots.png"), all.p.reg.plots, height = 8.25, width = 10.75, units = "in")


# Det intercept plots -----------------------------------------------------
psi.CAR.det.CAR.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = mean.p.bias, mod.name = "psiCARdetCAR.stan", plot.group = as.character(plt.grps[x])))  

psi.CAR.det.STD.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = mean.p.bias, mod.name = "psiCARdetSTD.stan", plot.group = as.character(plt.grps[x])))  

psi.STD.det.CAR.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = mean.p.bias, mod.name = "psiSTDdetCAR.stan", plot.group = as.character(plt.grps[x])))  

psi.STD.det.STD.plots <- map(seq_along(plt.grps), function (x)
  gen_summary_plots(base.data.frame = mean.p.bias, mod.name = "psiSTDdetSTD.stan", plot.group = as.character(plt.grps[x])))  

all.p.int.plots <- psi.CAR.det.CAR.plots[[1]] + xlab("") + labs(subtitle = "Occupancy (CAR on both components)") +
  psi.CAR.det.CAR.plots[[2]] + xlab("") + ylab("") + 
  psi.CAR.det.CAR.plots[[3]] + xlab("") + ylab("") + 
  psi.CAR.det.CAR.plots[[4]] + xlab("") + ylab("") +
  psi.CAR.det.STD.plots[[1]] + xlab("") + labs(subtitle = "Occupancy (CAR on occupancy only)") +
  psi.CAR.det.STD.plots[[2]] + xlab("") + ylab("") + 
  psi.CAR.det.STD.plots[[3]] + xlab("") + ylab("") + 
  psi.CAR.det.STD.plots[[4]] + xlab("") + ylab("") + 
  psi.STD.det.CAR.plots[[1]] + xlab("") + labs(subtitle = "Occupancy (CAR on detection only)") +
  psi.STD.det.CAR.plots[[2]] + xlab("") + ylab("") + 
  psi.STD.det.CAR.plots[[3]] + xlab("") + ylab("") + 
  psi.STD.det.CAR.plots[[4]] + xlab("") + ylab("") +
  psi.STD.det.STD.plots[[1]] + xlab("Occupancy probability") + labs(subtitle = "Occupancy (no CAR)") +
  psi.STD.det.STD.plots[[2]] + xlab("Detection probability") + ylab("") + 
  psi.STD.det.STD.plots[[3]] + xlab("Availability") + ylab("") + 
  psi.STD.det.STD.plots[[4]] + xlab("True value") + ylab("") + plot_layout(nrow = 4)

ggsave(here::here("Outputs", "detintplots.png"), all.p.int.plots, height = 8.25, width = 10.75, units = "in")


