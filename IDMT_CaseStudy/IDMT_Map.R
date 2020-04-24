library(tidybayes)
library(tidyverse)
library(sf)
library(tigris)
library(viridis)

stan.data <- read_rds("D:/Data/IDMTDataArchive/IDMTstanData.rds")

binom <- read_rds("D:/Data/IDMTDataArchive/IDMT_binom.rds") %>% 
  recover_types(stan.data) %>% 
  spread_draws(., `logit_psi.*`, regex=TRUE) %>%  
  to_broom_names() 
  
psiCARdetCAR <- read_rds("D:/Data/IDMTDataArchive/IDMT_psiCARdetCAR.rds")  %>% 
  recover_types(stan.data) %>% 
  spread_draws(., `logit_psi.*`, regex=TRUE) %>%  
  to_broom_names() 

psiSTDdetSTD <- read_rds("D:/Data/IDMTDataArchive/IDMT_psiSTDdetSTD.rds") %>% 
  recover_types(stan.data) %>% 
  spread_draws(., `logit_psi.*`, regex=TRUE) %>%  
  to_broom_names() 

prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

st <- c("ID", "MT")
tct <- reduce(
  map(st, function (x) tracts(x) %>% 
        as(. ,"sf") %>% 
        st_transform(., prj)),
  rbind)
IDMT <- states() %>%  as(. ,"sf") %>% 
  st_transform(., prj) %>% 
  filter(STUSPS %in% st)

tct$occSTD <- plogis(apply(psiSTDdetSTD[,4:572], 2, median))
tct$occCAR <- plogis(apply(psiCARdetCAR[,4:572], 2, median))
tct$binom <- plogis(apply(binom[,4:572], 2, median))
tct$diffCARvSTD <- tct$occCAR - tct$occSTD  
tct$diffCARvsBIN <- tct$occCAR - tct$binom

# Make maps ---------------------------------------------------------------
#the for mapping
theme_map <- function(...) {
  theme_minimal() +
    theme(
      #text = element_text(family = "Ubuntu Regular", color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # panel.grid.minor = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.major = element_line(color = "black", size = 0.002),
      panel.grid.minor = element_line(color = "black", size = 0.002),
      plot.background = element_rect(fill = "black", color = NA), 
      panel.background = element_rect(fill = "black", color = NA), 
      legend.background = element_rect(fill = "black", color = NA),
      panel.border = element_blank(),
      ...
    )
}
CAR.plot <- ggplot(data = tct, mapping = aes(fill = occCAR))+
  geom_sf() +
  scale_fill_viridis() +
  theme_map()

OCC.plot <- ggplot(data = tct, mapping = aes(fill = occSTD))+
  geom_sf() +
  scale_fill_viridis() +
  theme_map()

Bin.plot <- ggplot(data = tct, mapping = aes(fill = binom))+
  geom_sf() +
  scale_fill_viridis() +
  theme_map()


CARSTDplot <- ggplot(data = tct, mapping = aes(fill = diffCARvSTD))+
  geom_sf() +
  geom_sf(data = IDMT, colour= "white", fill = NA, size = 1.5) +
  scale_fill_scico("Difference in predicted probability", palette="cork") +
  theme_map() +
  guides(fill = guide_colorbar(nbin = 10, direction = "horizontal", title.position = "top", label.position = "bottom", barwidth = 20)) +
  theme(legend.position = "bottom", legend.justification = "center", legend.margin =margin(0.1,0.1,0.1,0.1),
        text=element_text(size = 20, family = "Helvetica", colour = "white", face = "bold"))

CARBinomPlot <- ggplot()+
  geom_sf(data = tct, mapping = aes(fill = diffCARvsBIN), size = 0.1) +
  geom_sf(data = IDMT, colour= "white", fill = NA, size = 1.5) +
  scale_fill_viridis("Difference in predicted probability") +
  theme_map() +
  guides(fill = guide_colorbar(nbin = 10, direction = "horizontal", title.position = "top", label.position = "bottom", barwidth = 20)) +
  theme(legend.position = "bottom", legend.justification = "center", legend.margin =margin(0.1,0.1,0.1,0.1),
        text=element_text(size = 20, family = "Helvetica", colour = "white", face = "bold"))

ggsave("G:/My Drive/Williamson_ConsEasementOcc/Outputs/CARbinomDif.png", CARBinomPlot, width = 10,  units = "in")

