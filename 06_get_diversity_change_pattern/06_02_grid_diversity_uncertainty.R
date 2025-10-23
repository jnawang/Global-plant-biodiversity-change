library(terra)
library(pbapply)
library(parallel)
#
rm(list=ls())
################################################calculate richness change first
#
richness_maxe_path <- '/Users/junnawang/UCDLab/Biodiversity/Results/2024/MAXE_2024_richness/'
#
richness_rfds_path <- '/Users/junnawang/UCDLab/Biodiversity/Results/2024/RFDS_2024_richness/'
#
richness_maxe_rfds_path <- '/Users/junnawang/UCDLab/Biodiversity/Results/2024/MAXE_RFDS_2024_richness/'

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#
# remove inland water cells
template_land <- terra::rast('/Users/junnawang/UCDLab/Biodiversity/data/climate/current/bio1.tif')
#
cal_richness_change <- function(richness_path) {
  richness_files <- list.files(path=richness_path)
  # remove files with _LogOut_ and _PresFut_
  richness_files <- richness_files[!grepl('_LogOut_', richness_files)]
  richness_files <- richness_files[!grepl('_PresFut_', richness_files)]
  richness_files <- richness_files[!grepl('richnessP', richness_files)]
  #
  richness <- terra::rast(file.path(richness_path, richness_files))
  #
  richness_pres <- terra::rast(file.path(richness_path, 'richnessP.tif'))
  #
  #
  richness <- mask(richness, template_land)
  richness_pres <- mask(richness_pres, template_land)
  richness_change <- richness / richness_pres
  
  # remove infinity
  richness_change[is.infinite(richness_change)] <- 0
  
  # set up names to raster layers
  names(richness_change) <- substr(richness_files, 1, nchar(richness_files) - 14)

  #
  return(richness_change)
}
#
# the order of richness files: future suitable, no dispersal, realized, realized_slow, realized_fast
richness_change_maxe <- cal_richness_change(richness_maxe_path)

richness_change_rfds <- cal_richness_change(richness_rfds_path)

richness_change_maxe_rfds <- cal_richness_change(richness_maxe_rfds_path)

###########################################calculate cellwise uncertainty and contributors next
# create a data frame of factors
df_richness_change <- data.frame(SDM=c(rep('MAXE', 96), rep('RFDS', 96), rep('MAXE_RFDS', 96)), 
                         ESM=rep(rep(c('ACCESS', 'ESM', 'MIROC6', 'MPIESM'), each=4), 18), 
                         CO2=rep(c('ssp126', 'ssp245', 'ssp370', 'ssp585'), 72) , 
                         migration=rep(rep(c('full', 'no', 'med', 'slow', 'fast', 'uni'), each=16), 3), change=0)

landcells <- which(!is.na(terra::values(template_land)))
#
uncertainty <- function(change_ratio) {
  df_richness_change$change <- change_ratio
#  print(change_ratio)
  if (sum(is.na(change_ratio)) == 0) {
    fit <- aov(change ~ SDM + ESM + CO2 + migration, df_richness_change, na.action = na.omit)
    return(c(sd(change_ratio, na.rm=T), summary(fit)[[1]][["Sum Sq"]]))     
  } else {
    return(rep(NA, 6))
  }
}
#
mtx_maxe <- terra::values(richness_change_maxe)[landcells, ]
mtx_rfds <- terra::values(richness_change_rfds)[landcells, ]
mtx_maxe_rfds <- terra::values(richness_change_maxe_rfds)[landcells, ]

mtx <- cbind(mtx_maxe, mtx_rfds, mtx_maxe_rfds)
rm(mtx_maxe)
rm(mtx_rfds)
rm(mtx_maxe_rfds)
#
cl <- parallel::makeCluster(detectCores() - 1)  # Use all cores minus 1
parallel::clusterExport(cl, varlist = "df_richness_change")
#
r1 <- pbapply(mtx[1:800000, ], 1, uncertainty, cl = cl)
r2 <- pbapply(mtx[800001:nrow(mtx), ], 1, uncertainty, cl = cl)
r1 <- cbind(r1, r2)
#### this takes 3-5 min. 
parallel::stopCluster(cl)
#
# spatial variation in richness change prediction
richness_change_sd <- template_land
terra::values(richness_change_sd)[landcells] <- r1[1,]
plot(richness_change_sd)
writeRaster(richness_change_sd, file='/Users/junnawang/UCDLab/Biodiversity/Results/2024/Uncertainty_richness/richness_change_std.tif', overwrite=T)
#
richness_change_SDM <- template_land
terra::values(richness_change_SDM)[landcells] <- r1[2,] / colSums(r1[2:5,])
plot(richness_change_SDM, main='Contribution of SDMs')
writeRaster(richness_change_SDM, file='/Users/junnawang/UCDLab/Biodiversity/Results/2024/Uncertainty_richness/richness_change_SDM.tif', overwrite=T)
# 
# this contribution is relatively small
richness_change_ESM <- template_land
terra::values(richness_change_ESM)[landcells] <- r1[3,] / colSums(r1[2:5,])
plot(richness_change_ESM, main='Contribution of earth system models')
writeRaster(richness_change_ESM, file='/Users/junnawang/UCDLab/Biodiversity/Results/2024/Uncertainty_richness/richness_change_ESM.tif', overwrite=T)
#
richness_change_CO2 <- template_land
terra::values(richness_change_CO2)[landcells] <- r1[4,] / colSums(r1[2:5,])
plot(richness_change_CO2, main='Contribution of CO2 scenarios')
writeRaster(richness_change_CO2, file='/Users/junnawang/UCDLab/Biodiversity/Results/2024/Uncertainty_richness/richness_change_CO2.tif', overwrite=T)
#
richness_change_migration <- template_land
terra::values(richness_change_migration)[landcells] <- r1[5,] / colSums(r1[2:5,])
plot(richness_change_migration, main='Contribution of migration')
writeRaster(richness_change_migration, file='/Users/junnawang/UCDLab/Biodiversity/Results/2024/Uncertainty_richness/richness_change_migration.tif', overwrite=T)
#
richness_change_residual <- template_land
terra::values(richness_change_residual)[landcells] <- r1[6,] / colSums(r1[2:6,])
plot(richness_change_residual, main='Contribution of residual')
writeRaster(richness_change_residual, file='/Users/junnawang/UCDLab/Biodiversity/Results/2024/Uncertainty_richness/richness_change_residual.tif', overwrite=T)

