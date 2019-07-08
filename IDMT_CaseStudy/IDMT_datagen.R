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
library(sf)
library(tigris)
library(tidyverse)
library(spdep)
library(rstan)
library(tidycensus)
options(mc.cores = parallel::detectCores())
options(tigris_use_cache = TRUE)
census_api_key('baaaffb5ed3accd8dfa53c6f827659d43fcdfa21') #get this from the census api webpage see help(census_api_key) for details
orig.data.folder <- "D:/Data/NCED_Data/OriginalData/"
gen.data.folder <- "D:/Data/NCED_Data/GeneratedData/"

prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# Get geography -----------------------------------------------------------
#load spatial data
st <- c("ID", "MT")

cty <- reduce(
  map(st, function (x) counties(x) %>% 
        as(. ,"sf") %>% 
        st_transform(., prj)),
  rbind)

tct <- reduce(
  map(st, function (x) tracts(x) %>% 
        as(. ,"sf") %>% 
        st_transform(., prj)),
    rbind)

bg <- reduce(
  map(st, function (x) block_groups(x) %>% 
        as(. ,"sf") %>% 
        st_transform(., prj)),
  rbind)

##Download the most recent NCED file
download.file("https://www.conservationeasement.us/data_downloads/NCED_09062018.zip", here::here("DataArchive","NCED09062018.zip"))
unzip(here::here("DataArchive","NCED09062018.zip"),exdir=here::here("DataArchive"))
ease <- read_sf(here::here("DataArchive/NCED_09062018/", "NCED_Polygons.shp")) %>% 
  subset(., state %in% st) %>% 
  st_transform(., prj)


# Get neighbors for tracts and blockgroups --------------------------------
gen_neighb_Matx <- function(geog){
  geog.order <- geog[order(geog$GEOID),] %>% as(., "Spatial");
  coords.geog <- coordinates(geog.order);
  geog.nearnb <- knn2nb(knearneigh(coords.geog, k = 1), row.names = geog.order$GEOID, sym=TRUE); #estimate distance to first neareset neighbor
  dis.nb.list <- dnearneigh(coords.geog, 0,  max( unlist(nbdists(geog.nearnb, coords.geog))),row.names = geog.order$GEOID); #neighbors within the minimum distance to ensure at least 1 neighbor
  geog.W.matx <- nb2mat(dis.nb.list, style="B", zero.policy = TRUE); #converts list of neighbors to matrix with 0s on diagonal (i.e., no self-neighbros)
}

tct.neighbors <- gen_neighb_Matx(tct)
bg.neighbors <- gen_neighb_Matx(bg)

# Tally easements in block group ------------------------------------------
int <- rowSums(st_contains(bg, ease, sparse = FALSE))
int.sf <- setNames(int, bg$GEOID)
block.group.easements <- data.frame(
  GEOID = names(int.sf),
  NumEase = int.sf,
  row.names = NULL
)
block.group.easements$EasePres <- ifelse(block.group.easements$NumEase > 0, 1, 0)
block.group.easements$GEOID <- as.character(block.group.easements$GEOID)
#reorder to match neighbor matrix
bg.easements <- block.group.easements[order(block.group.easements$GEOID),]

#get tract ids from blockgroup
bg.in.tct <- rethinking::coerce_index(str_sub(bg.easements$GEOID, 1, 11))
#Note that group_by reorders output; that's okay because the neighbor matrix code does that too

tct.ease <- block.group.easements %>% 
  group_by(tctid = str_sub(GEOID, 1,11)) %>% 
  summarise(tot_ease = sum(NumEase),
            tot_bg_w_ease = sum(EasePres)) %>% 
  mutate(EasePres = if_else(tot_ease > 0, 1, 0))
                       

# Get tract level census info ---------------------------------------------

sfips <- unique(str_sub(bg.easements$GEOID,1,2)) #create a lookup for tidycensus based on state sfips

med.inc.tct <- map_df(sfips, function(x) {
  get_acs(geography = "tract", variables = "B20002_001E", 
          state = x, year = 2015, output = "wide")})
colnames(med.inc.tct)[3] <- "medinctct"

ed.lvl.tct <- map_df(sfips, function(x) {
  get_acs(geography = "tract", variables = c("B15003_021E", "B15003_022E","B15003_023E", "B15003_024E", "B15003_025E", "B15003_001E"), 
          state = x, year = 2015, output = "wide") %>% 
    mutate(percDeg10tct = ((B15003_021E+B15003_022E+B15003_023E+B15003_024E+B15003_025E)/B15003_001E)*100)
})

med.age.tct <- map_df(sfips, function(x) {
  get_acs(geography = "tract", variables = "B01002_001E", 
          state = x, year = 2015, output = "wide")})
colnames(med.age.tct)[3] <- "medagetct"


# Get tract level institutional data --------------------------------------
#FROM THEOBALD 2014 PLOS One (had to do this manually because the files weren't downloading or unzipping properly from R)
download.file("http://csp-inc.org/public/NLUD2010_20140326.zip", here::here("DataArchive","NLUD2010.zip"))
unzip(here::here("DataArchive","NLUD2010.zip"),exdir=here::here("DataArchive"))
unzip(here::here("DataArchive","NLUD2010_W_20140326.zip"),exdir=here::here("DataArchive"))

