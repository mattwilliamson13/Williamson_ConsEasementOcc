library(sf)
library(tigris)
library(tidyverse)
library(spdep)
library(tidycensus)
library(raster)
library(FedData)
library(velox)
library(tidyverse)
library(magrittr)

options(mc.cores = parallel::detectCores())
options(tigris_use_cache = TRUE)
census_api_key('baaaffb5ed3accd8dfa53c6f827659d43fcdfa21') #get this from the census api webpage see help(census_api_key) for details
orig.data.folder <- "D:/Data/NCED_Data/OriginalData/"
gen.data.folder <- "D:/Data/NCED_Data/GeneratedData/"
rasterOptions(tmpdir = "D:/RastTemp/")

prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
#these objects need to be saved as .rds files for use in IDMT_datagen.R
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

# Get tract level institutional data --------------------------------------
#FROM THEOBALD 2014 PLOS One (had to do this manually because the files weren't downloading or unzipping properly from R)
#These aren't working because the drive is somehow full
download.file("http://csp-inc.org/public/NLUD2010_20140326.zip", "D:/Data/IDMTDataArchive/NLUD2010.zip")
unzip(here::here("DataArchive","NLUD2010.zip"),exdir=here::here("DataArchive"))
unzip(here::here("DataArchive","NLUD2010_W_20140326.zip"),exdir=here::here("DataArchive"))
rasterOptions(tmpdir ="D:/RastTemp/" )

lulc.W <- raster("D:/Theobald_2014_NLUD/NLUD2010_20140326/NLUD2010_W_20140326/nlud2010_w.tif")


lu_extract <- function(stfp, rstr){
  tct.reduced <- tct %>% subset(., STATEFP == stfp)
  lu.crop <- crop(rstr, tct.reduced)
  lu.vx <- velox(lu.crop) 
  tct.lu <- lu.vx$extract(tct.reduced)
  names(tct.lu) <- tct.reduced$GEOID
  return(tct.lu)
} 

tct_lu <- map(seq_along(sfips), function(x) lu_extract(stfp=sfips[x], rstr =lulc.W ))

wildness.rstr <- raster("D:/Data/wild50_2008/wild50_2008")
wild_extract <- function(stfp, rstr){
  tct.reduced <- tct %>% subset(., STATEFP == stfp) %>% st_transform(., st_crs(wildness.rstr))
  wild.crop <- crop(rstr, tct.reduced)
  wild.vx <- velox(wild.crop) 
  wild.lu <- wild.vx$extract(tct.reduced, small = TRUE)
  names(wild.lu) <- tct.reduced$GEOID
  return(wild.lu)
} 

wild_tct <- map(seq_along(sfips), function(x) wild_extract(stfp=sfips[x], rstr =wildness.rstr ))

rwr.raster <- raster("D:/MonumentData/L48_G1G2_1_Sept2013.tif")
rwr_extract <- function(stfp, rstr){
  tct.reduced <- tct %>% subset(., STATEFP == stfp) %>% st_transform(., st_crs(rstr))
  rwr.crop <- crop(rstr, tct.reduced)
  rwr.vx <- velox(rwr.crop) 
  rwr.lu <- rwr.vx$extract(tct.reduced, small=TRUE)
  names(rwr.lu) <- tct.reduced$GEOID
  return(rwr.lu)
} 
rwr_tct <- map(seq_along(sfips), function(x) rwr_extract(stfp=sfips[x], rstr =rwr.raster ))

nlcd_extract <- function(stfp){
  bg.reduced <- bg %>% subset(., STATEFP == stfp) %>% st_transform(. , "+proj=utm +datum=NAD83 +zone=12")%>% as(. , "Spatial")
  nlcd.rd <- get_nlcd(bg.reduced, label = paste0("bg",stfp, "Imp"), dataset = "impervious", raw.dir = "D:/Data/RastTemp/", extraction.dir = "D:/Data/IDMTDataArchive/")
  nlcd.vx <- velox(nlcd.rd)
  bg.reproj <- bg.reduced %>% as(., "sf") %>% st_transform(., st_crs(nlcd.rd))
  rm(bg.reduced)
  bg.nlcd <- nlcd.vx$extract(bg.reproj)
  names(bg.nlcd) <- bg.reproj$GEOID
  return(bg.nlcd)
  gc()
}

bg_nlcd <- map(seq_along(sfips), function(x) nlcd_extract(stfp=sfips[x]))
