library(tidybayes)
library(ggridges)
library(tidyverse)
library(viridis)
library(ggstance)
library(extrafont)
windowsFonts(myFont = windowsFont("TT Times New Rome"))

stan.data <- read_rds("D:/Data/IDMTDataArchive/IDMTstanData.rds")
design.df <- read_rds("D:/Data/IDMTDataArchive/designdf.rds")
df.tct <- read_rds("D:/Data/IDMTDataArchive/dftct.rds")
binom <- read_rds("D:/Data/IDMTDataArchive/IDMT_binom.rds") %>% 
  recover_types(stan.data) %>% 
  gather_draws(., `beta_psi.*`, regex=TRUE) %>%  
  to_broom_names() %>%
  median_hdi(.width = c(0.5, 0.8, 0.99)) %>% 
  mutate(modname = "Logistic")

psiCARdetCAR <- read_rds("D:/Data/IDMTDataArchive/IDMT_psiCARdetCAR.rds") %>% 
  recover_types(stan.data) %>% 
  gather_draws(., `beta_psi.*`, regex=TRUE) %>%  
  to_broom_names() %>% 
  median_hdi(.width = c(0.5, 0.8, 0.99)) %>% 
  mutate(modname = "Occupancy (CAR on both components)")

psiCARdetSTD <- read_rds("D:/Data/IDMTDataArchive/IDMT_psiCARdetSTD.rds") %>% 
  recover_types(stan.data) %>% 
  gather_draws(., `beta_psi.*`, regex=TRUE) %>%  
  to_broom_names() %>%
  median_hdi(.width = c(0.5, 0.8, 0.99)) %>% 
  mutate(modname = "Occupancy (CAR on occupancy only)")

psiSTDdetCAR <- read_rds("D:/Data/IDMTDataArchive/IDMT_psiSTDdetCAR.rds") %>% 
  recover_types(stan.data) %>% 
  gather_draws(., `beta_psi.*`, regex=TRUE) %>%  
  to_broom_names() %>%
  median_hdi(.width = c(0.5, 0.8, 0.99)) %>% 
  mutate(modname = "Occupancy (CAR on detection only)")

psiSTDdetSTD <- read_rds("D:/Data/IDMTDataArchive/IDMT_psiSTDdetSTD.rds") %>% 
  recover_types(stan.data) %>% 
  gather_draws(., `beta_psi.*`, regex=TRUE) %>%  
  to_broom_names() %>% 
  median_hdi(.width = c(0.5, 0.8, 0.99)) %>% 
  mutate(modname = "Occupancy (no CAR)")

all.models <- rbind(binom, psiCARdetCAR, psiCARdetSTD, psiSTDdetCAR, psiSTDdetSTD)
p <- ggplot(data=all.models, mapping = aes(x = estimate, y = term, color = modname, fill = modname)) +
  scale_fill_viridis_d(option = "D", alpha = 0.1, guide = guide_legend(reverse = TRUE)) +
  scale_color_viridis_d(option = "D", guide = guide_legend(reverse = TRUE)) +
  geom_pointintervalh(position = position_dodgev(height = .3)) +
  scale_y_discrete(name = NULL, labels = c("beta_psi[1]" = "Median income", "beta_psi[2]" = "Percent degree",
                                           "beta_psi[3]" = "Land use diversity", "beta_psi[4]" = "Rarity-weighted richness (max)",
                                           "beta_psi[5]" = "Wildness variance")) +
  theme_bw() +
  xlab("Estimate") +
  theme(legend.title = element_blank(),
        legend.position="none",  
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(size = 20, family = "myFont"))
  
ggsave("G:/My Drive/Williamson_ConsEasementOcc/Outputs/IDMT_paramestforPres.png", p, width = 10, height = 7, units = "in")


# Fitting lines -----------------------------------------------------------
state.id <- unlist(lapply(1:44, function(x) paste0("alpha_cty[", x, "]"))) #the first 44 counties are in ID

binom <- read_rds("D:/Data/IDMTDataArchive/IDMT_binom.rds") 

bin.posterior.mu.int <- as.matrix(binom, pars="mu_alpha")
bin.posterior.ints <- as.matrix(binom, pars="alpha_cty")
bin.posterior.E <- as.matrix(binom, pars = "beta_psi")[,4]

bin.int.medians <- apply(bin.posterior.ints, 2, median)

binom_fit <- data.frame(cty = names(bin.int.medians),
                        medInt = bin.int.medians,
                        X = rep(seq(from = min(design.df$mxRWR), to = max(design.df$mxRWR), length.out = 100), 
                                each = length(bin.int.medians))) %>% 
  mutate(fit_val = median(bin.posterior.mu.int) + medInt + (median(bin.posterior.E) * X),
         st = if_else(cty %in% state.id, "Idaho", "Montana"),
         modtype = "Logistic")

ggplot(binom_fit, mapping = aes(x = X, y = plogis(fit_val), group = cty, color = st))+
  geom_line() + facet_wrap(~st)

state.id.occ <- unlist(lapply(1:44, function(x) paste0("alpha_cty_occ[", x, "]"))) #the first 44 counties are in ID

psiCARdetCAR <- read_rds("D:/Data/IDMTDataArchive/IDMT_psiCARdetCAR.rds")

occ.posterior.mu.int <- as.matrix(psiCARdetCAR, pars="mu_alpha_occ")
occ.posterior.ints <- as.matrix(psiCARdetCAR, pars="alpha_cty_occ")
occ.posterior.E <- as.matrix(psiCARdetCAR, pars = "beta_psi")[,4]

occ.int.medians <- apply(occ.posterior.ints, 2, median)

occ_fit <- data.frame(cty = names(occ.int.medians),
                        medInt = occ.int.medians,
                        X = rep(seq(from = min(design.df$mxRWR), to = max(design.df$mxRWR), length.out = 100), 
                                each = length(occ.int.medians))) %>% 
  mutate(fit_val = median(occ.posterior.mu.int) + medInt + (median(occ.posterior.E) * X),
         st = if_else(cty %in% state.id.occ, "Idaho", "Montana"),
         modtype = "Occupancy \n(CAR on both components)")


both_fit <- rbind(binom_fit, occ_fit) %>% 
  mutate(rwrConv = (X * sd(df.tct$mxRWR)) + mean(df.tct$mxRWR))

mxRWR.plot <- ggplot(both_fit, mapping = aes(x = rwrConv, y = plogis(fit_val), group = cty, color = modtype)) +
  geom_line() +
  scale_color_viridis_d(option = "D", begin = 0, end = 0.4, guide = guide_legend(reverse = TRUE)) +
  facet_wrap(~st) +
  theme_bw() +
  xlim(c(0,.5))+
  labs(y = "Probability of easement", x="Maximum rarity-weighted richness") +
  theme(legend.title = element_blank(),
        legend.text.align = 0,
        legend.justification=c(1,0), 
        legend.position= "none",  
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        text = element_text(size = 20, family = "myFont"),
        panel.spacing.x = unit(1,"lines") )



bin.posterior.mu.int <- as.matrix(binom, pars="mu_alpha")
bin.posterior.ints <- as.matrix(binom, pars="alpha_cty")
bin.posterior.E <- as.matrix(binom, pars = "beta_psi")[,2]

bin.int.medians <- apply(bin.posterior.ints, 2, median)

binom_fit2 <- data.frame(cty = names(bin.int.medians),
                        medInt = bin.int.medians,
                        X = rep(seq(from = min(design.df$percDeg10tct), to = max(design.df$percDeg10tct), length.out = 100), 
                                each = length(bin.int.medians))) %>% 
  mutate(fit_val = median(bin.posterior.mu.int) + medInt + (median(bin.posterior.E) * X),
         st = if_else(cty %in% state.id, "Idaho", "Montana"),
         modtype = "Logistic")


occ.posterior.mu.int <- as.matrix(psiCARdetCAR, pars="mu_alpha_occ")
occ.posterior.ints <- as.matrix(psiCARdetCAR, pars="alpha_cty_occ")
occ.posterior.E <- as.matrix(psiCARdetCAR, pars = "beta_psi")[,2]

occ.int.medians <- apply(occ.posterior.ints, 2, median)

occ_fit2 <- data.frame(cty = names(occ.int.medians),
                      medInt = occ.int.medians,
                      X = rep(seq(from = min(design.df$percDeg10tct), to = max(design.df$percDeg10tct), length.out = 100), 
                              each = length(occ.int.medians))) %>% 
  mutate(fit_val = median(occ.posterior.mu.int) + medInt + (median(occ.posterior.E) * X),
         st = if_else(cty %in% state.id.occ, "Idaho", "Montana"),
         modtype = "Occupancy \n(CAR on both components)")


both_fit2 <- rbind(binom_fit2, occ_fit2) %>% 
  mutate(degConv = (X * sd(df.tct$percDeg10tct)) + mean(df.tct$percDeg10tct))

degree.plot <- ggplot(both_fit2, mapping = aes(x = degConv, y = plogis(fit_val), group = cty, color = modtype)) +
  geom_line() +
  scale_color_viridis_d(option = "D", begin = 0, end = 0.4, guide = guide_legend(reverse = TRUE)) +
  facet_wrap(~st) +
  theme_bw() +
  labs(y = "Probability of easement", x="Percent of population with a college degree") +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text = element_blank(),
        text = element_text(size = 20, family = "myFont"),
        panel.spacing.x = unit(1,"lines"))

case.study.marg.effect <- cowplot::plot_grid(mxRWR.plot,
                                             degree.plot,
                                             label_fontfamily = "myFont",
                                             label_x = 0.95,
                                             label_y = c(1, 1.05),
                                             scale = 0.9,
                                             nrow = 2, align = "v", axis = "l")

fig6 <- cowplot::plot_grid(case.study.marg.effect, CARBinomPlot,
                           ncol = 2, rel_widths = c(1.3, 1),
                           labels = c("", "C"),
                           label_size = 12,
                           label_x = 0.9,
                           label_y = 0.75,
                           label_fontfamily = "myFont")

ggsave("G:/My Drive/Williamson_ConsEasementOcc/Outputs/IDMT_margeffect.png", fig6, width = 8, height = 6.75, units = "in")
