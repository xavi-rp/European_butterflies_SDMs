###############################################
####      Modelling distributions of       ####
####        European butterflies           ####
###############################################

wd <- "/Users/xavi_rp/Documents/D5_FFGRCC/butterflies_SDMs"
setwd(wd)

library(tidyr)
library(ENMeval)
library(raster)
library(dplyr)
library(dismo)


## Predictors ####

worldclim_path0.5 <- "/Users/xavi_rp/Documents/MinBA_models/wc0.5" #  1km

worldclim_files <- list.files(worldclim_path0.5, full.names = TRUE)
worldclim_files <- worldclim_files[grepl("bil$", worldclim_files) | grepl("tif$", worldclim_files)]
worldclim_files

worldclim_all <- stack(worldclim_files)

eur_coords <- c(-13, 48, 35, 72)
worldclim_all <- crop(worldclim_all, eur_coords)
worldclim_all
plot(worldclim_all[[3]])


## Background points ####

bckgr <- randomPoints(worldclim_all[[1]], n = 10000)

# check coverage
plot(worldclim_all[[1]])
points(bckgr, pch = 20, cex = 0.2)






