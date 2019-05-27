## take an sf object of polygons, estimate the minimum distance necessary to ensure at least 1 neighbor and return the precision matrix.
library(tidyverse)
library(sf)
library(spdep)

# Generate the spatial precision matrix --------------------------------------------
gen_sp_precMatx <- function(geog, rho, tau){
  geog.order <- geog[order(geog$GEOID),] %>% as(., "Spatial");
  coords.geog <- coordinates(geog.order);
  geog.nearnb <- knn2nb(knearneigh(coords.geog, k = 1), row.names = geog.order$GEOID); #estimate distance to first neareset neighbor
  dis.nb.list <- dnearneigh(coords.geog, 0,  max( unlist(nbdists(geog.nearnb, coords.geog))),row.names = geog.order$GEOID); #neighbors within the minimum distance to ensure at least 1 neighbor
  geog.W.matx <- nb2mat(dis.nb.list, style="B", zero.policy = TRUE); #converts list of neighbors to matrix with 0s on diagonal (i.e., no self-neighbros)
  D.geog <- diag(rowSums(geog.W.mat)); #creates a matrix with the number of neighbors as the diagonal
  prec.matx <- tau*(D.geog - rho*geog.W.matx)
}