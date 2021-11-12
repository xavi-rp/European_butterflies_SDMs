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
library(data.table)


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

worldclim_all_data <- as.data.frame(worldclim_all)
worldclim_all_data <- worldclim_all_data[complete.cases(worldclim_all_data), ]
head(worldclim_all_data)
nrow(worldclim_all_data)


## Background points ####

bckgr <- randomPoints(worldclim_all[[1]], 
                      #n = 10000)
                      n = 15000)
bckgr <- as.data.frame(bckgr)
head(bckgr)
nrow(bckgr)

write.csv(bckgr, "background_points.csv", row.names = FALSE)
bckgr <- read.csv("background_points.csv", header = TRUE)

# check coverage
pdf("background_15k.pdf")
plot(worldclim_all[[1]])
points(bckgr, pch = 20, cex = 0.2)
dev.off()


## Butterflies Occurrences ####

occs_all <- read.csv("sp_records_20211109.csv", header = TRUE)
head(occs_all)
table(occs_all$species)

taxons <- unique(occs_all$sp2)
spcies <- data.frame(taxons, sps = unique(occs_all$species))
spcies

for (t in taxons){
  #print(t)
  t0 <- Sys.time()
  sps <- spcies[spcies$taxons == t, "sps"]
  print(paste0("running... ", sps))
  if(!dir.exists(paste0("models_", t))) {
    dir2save <- paste0("models_", t, "/")
    dir.create(dir2save)
  }
  occs_i <- occs_all[occs_all$sp2 %in% t, c("decimalLongitude", "decimalLatitude")]
  
  occs_i_shp <- SpatialPointsDataFrame(coords = occs_i[, c("decimalLongitude", "decimalLatitude")],
                                       data = data.frame(sp = rep(1, nrow(occs_i))),
                                       proj4string = CRS("+init=EPSG:4326"))
  #names(occs_i_shp) <- t
  occs_i_rstr <- rasterize(occs_i_shp, worldclim_all[[1]], field = "sp", background = 0)
  #names(occs_i_rstr) <- t
  #occs_i_rstr <- occs_i_rstr[[2]]
  occs_i_rstr <- mask(occs_i_rstr, worldclim_all[[1]])
  
  #assign(paste0(t, "_rstr"), occs_i_rstr)
  #print(sum(getValues(occs_i_rstr) == 1, na.rm = T))
  
  
  ## occurrences for training/testing
  sps_data <- stack(occs_i_rstr, worldclim_all) 
  sps_data <- as.data.frame(sps_data)
  sps_data$raster_position <- 1:nrow(sps_data)
  
  # data set for presences
  sps_data_presences <- sps_data[sps_data$layer == 1, ]
  sps_data_presences <- sps_data_presences[complete.cases(sps_data_presences), ]
  rm(sps_data); gc()
  
  # data set for pseudo-absences
  sps_data_absences <- as.data.frame(raster::extract(worldclim_all, bckgr, cellnumbers = TRUE))
  sps_data_absences <- sps_data_absences[!sps_data_absences$cells %in% sps_data_presences$raster_position, ]
  names(sps_data_absences)
  
  ## Running ENMeval (https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0.0-vignette.html)
  modl <- ENMevaluate(occs = sps_data_presences[, names(sps_data_presences) %in% names(worldclim_all)], 
                      envs = NULL, 
                      bg = sps_data_absences[, names(sps_data_absences) %in% names(worldclim_all)], 
                      algorithm = 'maxnet', 
                      partitions = 'block', 
                      tune.args = list(
                        fc = c("L","LQ","LQH","H"),
                        #fc = c("L","LQ"),
                        rm = c(1, 2, 5)
                        #rm = 1:2
                        ),
                      parallel = TRUE,
                      numCores = 7
                      )
  modl
  modl@results
  View(modl@results)
  write.csv(modl@results, file = paste0(dir2save, "ENMeval_results_", t, ".csv"))
  evalplot.stats(e = modl, stats = "or.mtp", color = "fc", x.var = "rm")
  
  occurrences_train <- modl@occs
  occurrences_test <- nrow(modl@occs.testing)  # none because crossvalidation
  background_pts <- nrow(modl@bg)
  
  # selecting optimal model
  results <- eval.results(modl)
  results
  View(results)
  optimal <- results %>% filter(delta.AICc == 0)
  optimal
  
  modl_args <- eval.models(modl)[[optimal$tune.args]]
  modl_args$betas
  str(modl_args)
  
  dev.off()
  pdf(paste0(dir2save, "opt_model_RespCurves", t, ".pdf"))
  plot(modl_args, type = "cloglog")
  # And these are the marginal response curves for the predictor variables wit non-zero 
  # coefficients in our model. We define the y-axis to be the cloglog transformation, which
  # is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
  # (Phillips et al. 2017).
  dev.off()
  
  modl <- modl@models[[optimal$tune.args]]
  gc()
  
  save(modl, file = paste0(dir2save, "opt_model_", t, ".RData"))
  
  # making predictions
  worldclim_all_data_kk <- worldclim_all_data
  sps_predictions_maxent <- predict(object = modl, 
                                    newdata = worldclim_all_data, 
                                    clamp = TRUE,
                                    type = c("cloglog")
                                    )
  sps_predictions_maxent <- as.data.table(sps_predictions_maxent)
  head(sps_predictions_maxent)
  range(sps_predictions_maxent)
  nrow(sps_predictions_maxent)
  
  worldclim_all_data <- as.data.table(as.data.frame(worldclim_all[[1]]))
  worldclim_all_data$raster_position <- 1:nrow(worldclim_all_data)

  worldclim_all_data1 <- worldclim_all_data
  worldclim_all_data1 <- worldclim_all_data1[complete.cases(worldclim_all_data1), ]
  
  worldclim_all_data1[, predictions := sps_predictions_maxent$V1]
  
  
  worldclim_all_data <- merge(worldclim_all_data[, "raster_position", with = FALSE], 
                              worldclim_all_data1[, .SD, .SDcols = c("raster_position", "predictions")], 
                              by = "raster_position", all.x = TRUE)

  sps_preds_rstr <- worldclim_all[[1]]
  sps_preds_rstr <- setValues(sps_preds_rstr, worldclim_all_data$predictions)
  
  #pdf("sps_predictions_maxent_kk.pdf", width = 20, height = 15)
  #par(mfrow = c(1, 2))
  #plot(sps_preds_rstr, zlim = c(0, 1))
  #plot(occs_i_shp, add = TRUE, col = "black")
  #plot(sps_preds_rstr, zlim = c(0, 1))
  #dev.off()
  
  
  #BI_mxnt <- ecospat::ecospat.boyce(fit = sps_preds_rstr,
  #                                  obs = linaria_pres_test_coords, 
  #                                  nclass = 0, 
  #                                  window.w = "default", 
  #                                  res = 100, 
  #                                  PEplot = TRUE)
  
  ## Creating presence/absence map
  # Threshold: minimum presence
  
  threshold1 <- min(extract(sps_preds_rstr, occs_i_shp))
  threshold1 <- quantile(extract(sps_preds_rstr, occs_i_shp), 0.1) # sensitivity = 0.9
  threshold1
  
  threshold2 <- dismo::threshold(evaluate(extract(sps_preds_rstr, occs_i_shp), extract(sps_preds_rstr, bckgr))) # sensitibity default 0.9
  threshold2
  
  a <- c(0, threshold1, 0)
  b <- c(threshold1, 1, 1)
  thr <- rbind(a, b)
  
  sps_preds_rstr_pres_abs <- reclassify(sps_preds_rstr, rcl = thr, filename = '', include.lowest = FALSE, right = TRUE)
  plot(sps_preds_rstr_pres_abs)
  

  pdf(paste0(dir2save, "sps_predictions_maxent_", t, ".pdf"), width = 18, height = 15)
  par(mar = c(6, 8, 6, 8), oma = c(4,0,8,0))
  par(mfrow = c(2, 2))
  plot(sps_preds_rstr, zlim = c(0, 1), main = "Occurrences (GBIF - 1km)", cex.main = 2, cex.sub = 1.5, legend = FALSE)
  plot(occs_i_shp, add = TRUE, col = "black")
  plot(sps_preds_rstr, zlim = c(0, 1), main = "MaxEnt predictions (cloglog)", cex.main = 2, cex.sub = 1.5)
  plot(sps_preds_rstr_pres_abs, main = "Presence-Absence", sub = "Threshold = 10th Percentile training presences", 
       cex.main = 2, cex.sub = 1.5, legend = FALSE)
  title(list(paste0(sps),
             cex = 4), 
        line = 1, outer = TRUE)
  
  dev.off()
  
  print(paste0(t, " run in: ", (Sys.time() - t0)))
}
















