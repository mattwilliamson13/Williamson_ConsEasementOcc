library(sf)
library(tigris)
library(tidyverse)
library(spdep)
library(rstan)
library(tidycensus)
library(raster)
library(tidyverse)
library(magrittr)
library(stringi)

options(mc.cores = parallel::detectCores())
options(tigris_use_cache = TRUE)
census_api_key('YOUR_API_KEY_HERE') #get this from the census api webpage see help(census_api_key) for details
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
download.file("https://www.conservationeasement.us/data_downloads/NCED_09062018.zip", "D:/Data/IDMTDataArchive/NCED09062018.zip")
unzip("D:/Data/IDMTDataArchive/NCED09062018.zip",exdir="D:/Data/IDMTDataArchive")
ease <- read_sf("D:/Data/IDMTDataArchive/NCED_09062018/NCED_Polygons.shp") %>% 
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

#Note that group_by reorders output; that's okay because the neighbor matrix code does that too

tct.ease <- block.group.easements %>% 
  group_by(tctid = str_sub(GEOID, 1,11)) %>% 
  summarise(tot_ease = sum(NumEase),
            tot_bg_w_ease = sum(EasePres)) %>% 
  mutate(EasePres = if_else(tot_ease > 0, 1, 0))


# Get tract level census info ---------------------------------------------

sfips <- unique(str_sub(bg$GEOID,1,2)) #create a lookup for tidycensus based on state sfips

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


# Get maximum rarity weighted richness ------------------------------------

rwr.tract <- readRDS("D:/Data/IDMTDataArchive/rwr_tract.rds") #Users will have to contact NatureServe for this data
max.rwr <- map(seq_along(rwr.tract), function(x)
  map(seq_along(rwr.tract[[x]]), function(y) max(rwr.tract[[x]][[y]], na.rm = TRUE))
)  
tct.names <- map(seq_along(rwr.tract), function(x) names(rwr.tract[[x]])) %>% unlist()
rwr.mx.df <- data.frame(GEOID = tct.names,
                        mxRWR = unlist(max.rwr))


# Get variance of wildness ------------------------------------------------

wild.tract <- readRDS("D:/Data/IDMTDataArchive/wildness_tract.rds") #wildness data from Aplet et al. 2020
var.wild <- map(seq_along(wild.tract), function(x)
  map(seq_along(wild.tract[[x]]), function(y) sd(wild.tract[[x]][[y]], na.rm = TRUE)^2)
)  
wild.names <-  map(seq_along(wild.tract), function(x) names(wild.tract[[x]])) %>% unlist()
wild.var.df <- data.frame(GEOID = wild.names,
                          varWILD = unlist(var.wild))

# Calculate Theobalds entropy index ---------------------------------------
lu.tract <- readRDS("D:/Data/IDMTDataArchive/LU_extract.rds")

entropy_prop <- function(lu.vect){
  lu.char = as.character(lu.vect[[1]])
  prop.live = mean(stri_detect_regex(lu.char, "^21."), na.rm=TRUE)
  prop.work1 = mean(stri_detect_regex(lu.char, "^22."), na.rm=TRUE)
  prop.work2 = mean(stri_detect_regex(lu.char, "^23."), na.rm=TRUE)
  prop.work3 = mean(stri_detect_regex(lu.char, "^24."), na.rm=TRUE)
  prop.work = sum(prop.work1, prop.work2, prop.work3, na.rm = TRUE)
  prop.play = mean(stri_detect_regex(lu.char, "^4."), na.rm=TRUE)
  prop.shop = mean(stri_detect_regex(lu.char, "^222"), na.rm=TRUE)
  df = data.frame(GEOID = names(lu.vect),
                  plive = if_else(prop.live == 0,1e-19, prop.live),
                  pwork = if_else(prop.work == 0,1e-19, prop.work),
                  pplay =if_else(prop.play == 0,1e-19, prop.play),
                  pshop = if_else(prop.shop == 0,1e-19, prop.shop)
                  
  )
}


end.df.tct <- map(seq_along(lu.tract), function (x) seq_along(lu.tract[x][[1]]) %>% 
                    map(. , function(y) entropy_prop(lu.tract[x][[1]][y])))
entropy.df.tct <-  map(seq_along(end.df.tct), function (x) end.df.tct[[x]] %>% bind_rows()) %>% 
  bind_rows(.) %>% 
  mutate(E = -((plive * log(plive) + pwork * log(pwork) + pplay * log(pplay) + pshop * log(pshop))/log(4)))

entropy.df.tct <- readRDS("D:/Data/IDMTDataArchive/landuseent_tract_df.rds")

# load block group covariate ----------------------------------------------

impervious.bg <- readRDS("D:/Data/IDMTDataArchive/impervious_surf_bg.rds")
mean.imp.surf <- map(seq_along(impervious.bg), function(x)
  map(seq_along(impervious.bg[[x]]), function(y) median(impervious.bg[[x]][[y]], na.rm = TRUE))
)  
bg.names <- map(seq_along(impervious.bg), function(x) names(impervious.bg[[x]])) %>% unlist()

imp.surf.mean.df <- data.frame(GEOID = bg.names,
                               mnImpSurf = unlist(mean.imp.surf))


# make occupancy design matrix --------------------------------------------

df.tct <- med.inc.tct %>% 
  left_join(., ed.lvl.tct) %>% 
  left_join(., entropy.df.tct) %>% 
  left_join(., rwr.mx.df) %>% 
  left_join(., wild.var.df) %>% 
  mutate_all(~if_else(is.na(.x), median(.x, na.rm = TRUE), .x))    
#check to make sure the tract variables are ordered the same as the easement and neighborhood variables
identical(df.tct$GEOID, tct.ease$tctid)
saveRDS(df.tct, "D:/Data/IDMTDataArchive/dftct.rds" )


design.df <- df.tct %>% 
  dplyr::select(., medinctct, percDeg10tct, E, mxRWR, varWILD) %>% 
  mutate_all(~scale(.x))
saveRDS(design.df, "D:/Data/IDMTDataArchive/designdf.rds" )
occ.des.matx <- as.matrix(design.df) #using cty level varying intercept


# Make bg design matrix ---------------------------------------------------
imp.surf.mean.df$GEOID <- as.character(imp.surf.mean.df$GEOID)
imp.surf.mean.df.ordered <- imp.surf.mean.df[order(imp.surf.mean.df$GEOID),] 

bg.df <- st_set_geometry(bg, NULL) %>% 
  dplyr::select(., GEOID, ALAND)
bg.df$ALAND <- log10(as.numeric(bg.df$ALAND)) #log of area for sampling effort

bg.des.df <- imp.surf.mean.df.ordered %>% 
  left_join(., bg.df) %>% 
  mutate_at(c("mnImpSurf", "ALAND"), scale)

identical(as.character(bg.des.df$GEOID), bg.easements$GEOID)  
  
bg.des.matx <- as.matrix(bg.des.df[,c(2:3)]) #using cty level varying intercept


# Set up indices ----------------------------------------------------------

#get tract ids from blockgroup
tct.in.cty <- rethinking::coerce_index(str_sub(tct.ease$tctid, 1, 5))
bg.in.cty <-  rethinking::coerce_index(str_sub(imp.surf.mean.df.ordered$GEOID, 1, 5))
survey_summary <- bg %>% 
  group_by(., STATEFP, COUNTYFP, TRACTCE) %>% 
  summarise(count = n())
n_survey <- as.vector(survey_summary$count)
total_surveys <- nrow(bg)

survey.df <- tibble(site = rep(1:nrow(tct), n_survey),
                    siteID = bg.easements$GEOID) %>% 
  mutate(y = bg.easements$EasePres)

# get start and end indices to extract slices of y for each site
start_idx <- rep(0, nrow(tct))
end_idx <- rep(0, nrow(tct))
for (i in 1:nrow(tct)) {
  if (n_survey[i] > 0) {
    site_indices <- which(survey.df$site == i)
    start_idx[i] <- site_indices[1]
    end_idx[i] <- site_indices[n_survey[i]]
  }
}

any_seen <- rep(0, nrow(tct))
for (i in 1:nrow(tct)) {
  if (n_survey[i] > 0) {
    any_seen[i] <- max(survey.df$y[start_idx[i]:end_idx[i]])
  }
}

identical(any_seen, tct.ease$EasePres)
# Package for stan --------------------------------------------------------

stan_data <- list(
  n_cty = nrow(cty),
  tctincty = tct.in.cty,
  bgincty = bg.in.cty,
  n_site = nrow(tct.ease),
  m_psi = ncol(occ.des.matx),
  X_tct = occ.des.matx,
  total_surveys = nrow(bg.easements),
  m_p = ncol(bg.des.matx),
  X_bg = bg.des.matx, 
  site = survey.df$site,
  y = survey.df$y,
  start_idx = start_idx,
  end_idx = end_idx,
  any_seen = any_seen,
  n_survey = n_survey,
  W_tct = tct.neighbors,
  W_n_tct = sum(tct.neighbors/2),
  W_bg = bg.neighbors,
  W_n_bg = sum(bg.neighbors/2)
)
saveRDS(stan_data, "D:/Data/IDMTDataArchive/IDMTstanData.rds")
