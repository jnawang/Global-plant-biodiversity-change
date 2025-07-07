### code to calcualte warming rate and other environmental factors for each species. 
library(terra)
library(pbapply)
library(gstat)
rm(list=ls())
gc()
# 
SDM.path  <- '../../SDM2024MAXE/Results/'
#species_map <- list.files(path = SDM.path, pattern = '.*tif$')
table    <- read.csv('SDM_perform_merged_good.csv')
spsname  <- table$sps
nsps     <- length(spsname)
#
# read elevation and latitudinal warming velocity files
warm.rate.lat <- terra::rast('maps/warm_vel_lat2024.tif')
warm.rate.ele <- terra::rast('maps/warm_vel_ele2024.tif')

#
mountain    <- terra::rast('maps/k3binary5m.tif')
mountain    <- terra::project(mountain, warm.rate.lat)
HFI         <- terra::rast('maps/HFP2009_5m.tif')
HFI         <- terra::project(HFI, warm.rate.lat)
temp.mean   <- terra::rast('maps/bio1.tif')
temp.mean   <- terra::project(temp.mean, warm.rate.lat)
temp.var    <- terra::rast('maps/bio4.tif')
temp.var    <- terra::project(temp.var, warm.rate.lat)
elev.het    <- terra::rast('maps/elev_heter.tif')
elev.het    <- terra::project(elev.het, warm.rate.lat)

# calculate warming rate within species distribution area
warm_rate_species <- function(i) {
  # cropped mountain area / species area; it will be easier to use species. 
  sps       <- terra::rast(paste0(SDM.path, 'RangePres/', spsname[i], '.tif'))
  sps.fut   <- terra::rast(paste0(SDM.path, 'ESM_245_2081-2100/RangeFut/', spsname[i], '.tif'))     # I use a medium scenario for all. 
  sps[sps.fut==1] <- 1  
  sps[sps == 0]   <- NA
  # change the map to lant-long
  sps       <- trim(sps)
  sps       <- terra::project(sps, mountain)    # why projection changes cell values? because the value is continuous. It is important to have align here
  #
  mount     <- terra::crop(mountain, sps)
  sps.mount <- terra::mask(sps, mount)                     # when do mask, the extent has to match exactly!
  sps.flat  <- terra::mask(sps, mount, inverse=T)
  
  # go through 16 scenarios together: in mount area first, then in flat area. No! Do not use loop, just use multi-layer spatRaster
  cropped   <- terra::crop(warm.rate.ele, sps.mount)
  masked    <- terra::mask(cropped, sps.mount)
  EleVeloT  <- terra::global(masked, fun = 'mean', na.rm = TRUE)[,1]
  cropped   <- terra::crop(warm.rate.lat, sps.flat)
  masked    <- terra::mask(cropped, sps.flat)
  LatVeloT  <- terra::global(masked, fun = 'mean', na.rm = TRUE)[,1]
 
  #
  area.in.mount <- sum(terra::expanse(sps.mount, unit = 'km', transform=TRUE))
  area.on.flat  <- sum(terra::expanse(sps.flat, unit = 'km', transform=TRUE))
  frac.mount    <- area.in.mount/(area.on.flat+area.in.mount)
  
  # for the rest parameters in mount area first, then in flat area. 
  cropped   <- terra::crop(HFI, sps)
  masked    <- terra::mask(cropped, sps.mount)  # this is not efficient, need to resample for every image! 
  HFI.ele <- terra::global(masked, fun = 'mean', na.rm = TRUE)[1,1]
  masked    <- terra::mask(cropped, sps.flat)   # this is not efficient, need to resample for every image! 
  HFI.lat  <- terra::global(masked, fun = 'mean', na.rm = TRUE)[1,1]
  #
  cropped   <- terra::crop(temp.mean, sps)
  masked    <- terra::mask(cropped, sps.mount)
  temp.mean.ele <- terra::global(masked, fun = 'mean', na.rm = TRUE)[1,1]
  masked    <- terra::mask(cropped, sps.flat)
  temp.mean.lat <- terra::global(masked, fun = 'mean', na.rm = TRUE)[1,1]
  # 
  cropped   <- terra::crop(temp.var, sps)
  masked    <- terra::mask(cropped, sps.mount)
  temp.var.ele <- terra::global(masked, fun = 'mean', na.rm = TRUE)[1,1]
  masked    <- terra::mask(cropped, sps.flat)   
  temp.var.lat  <- terra::global(masked, fun = 'mean', na.rm = TRUE)[1,1]
  #
  cropped   <- terra::crop(elev.het, sps)
  masked    <- terra::mask(cropped, sps.mount)
  elev.het.ele <- terra::global(masked, fun = 'mean', na.rm = TRUE)[1,1]
  masked    <- terra::mask(cropped, sps.flat)   
  elev.het.lat <- terra::global(masked, fun = 'mean', na.rm = TRUE)[1,1]     
  
  return(c(spsname[i], area.in.mount, area.on.flat, frac.mount, HFI.ele, HFI.lat, temp.mean.ele, temp.mean.lat, temp.var.ele, temp.var.lat, elev.het.ele, elev.het.lat, EleVeloT, LatVeloT))
}

#
command_args <- commandArgs(trailingOnly = TRUE)
spset        <- as.numeric(command_args[1])

sps.range <- ((spset-1)*10+1):min(spset*10, nsps)

VeloT <- pbsapply(sps.range, warm_rate_species)
write.csv(VeloT, paste0('result/', spset, '.csv'))

