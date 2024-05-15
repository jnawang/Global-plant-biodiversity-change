# calculate species richness distribution based on SDM results
# change this code to do parallel for different scenarios
library(terra, sp)

command_args <- commandArgs(trailingOnly = TRUE)
id.s <- as.numeric(command_args[1])      # vary from 13 to 16
# id.s <- 15

gcms  <- c("ESM")
ssps  <- c(126, 245, 370, 585)
times <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
nssp  <- length(ssps)
ntime <- length(times)
scenarios <- c()
for (itime in 1:ntime) {
  for (issp in 1:nssp) {
    scenarios[(itime-1)*nssp+issp] <- paste(gcms[1], ssps[issp], times[itime], sep = '_')
  }
}

SDM.path <- '../../SDM072023/'
# SDM.path <- 'C:/Users/jnawang/Box/Biodiversity/SDM/'   # TEST LINE

# read the raster for all species
path.present    <- paste0(SDM.path, 'Results/RangePres/')
species.present <- list.files(path.present, pattern = ".tif$")  

path.future     <- paste0(SDM.path, 'Results/', scenarios[id.s], '/RangeFut/')
species.future  <- list.files(path.future, pattern = ".tif$")  

path.futureR     <- paste0(SDM.path, 'Realized/', scenarios[id.s], '/')
species.futureR  <- list.files(path.futureR, pattern = ".tif$")        


# find the species that were simulated in the present and future states
species <- intersect(intersect(species.present, species.future), species.futureR)

# find the well-modeled species
species.SDMgood <- read.csv(paste0(SDM.path, 'merged_stats/SDM_perform_merged.csv'))
species <- intersect(species, paste0(species.SDMgood$sps, '.tif'))
#
nsps    <- length(species)

richness <- terra::rast(paste0(SDM.path, 'preds/bios_pres/bio1.tif'))
richness[!is.na(richness)] <- 0
richness.present  <- richness
richness.future   <- richness
richness.futureR  <- richness
richness.presfut  <- richness
Jaccad.local   <- richness
Jaccad.extinct <- richness
Jaccad.migrate <- richness

#
# I hesitate to use apply or parallel because this is sum. I can use parallel. 
# length(SD_pred)
for (i in 1:nsps) {
  if (i %% 100 == 0) {print(i)}
  present.map <- terra::rast(paste0(path.present, species[i]))
  future.map  <- terra::rast(paste0(path.future,  species[i]))
  futureR.map <- terra::rast(paste0(path.futureR, species[i]))
  # present map
  xy <- terra::xyFromCell(present.map, which(terra::values(present.map)==1))
  id.present <- terra::cellFromXY(richness, xy)
  
  # future map
  xy <- terra::xyFromCell(future.map, which(terra::values(future.map)==1))
  id.future <- terra::cellFromXY(richness, xy)
  
  # future realized map
  xy <- terra::xyFromCell(futureR.map, which(terra::values(futureR.map)==1))
  id.futureR <- terra::cellFromXY(richness, xy)  

  # count number of species for each cell  
  richness.present[id.present] <- richness.present[id.present] + 1
  richness.future[id.future]   <- richness.future[id.future] + 1
  richness.futureR[id.futureR] <- richness.futureR[id.futureR] + 1
  id.presfut <- union(id.present, id.futureR)
  richness.presfut[id.presfut] <- richness.presfut[id.presfut] + 1
  #
  id <- intersect(id.present, id.futureR)
  Jaccad.local[id] <- Jaccad.local[id] + 1
  id <- setdiff(id.present, id.futureR)
  Jaccad.extinct[id] <- Jaccad.extinct[id] + 1
  id <- setdiff(id.futureR, id.present)
  Jaccad.migrate[id] <- Jaccad.migrate[id] + 1
}

# for testing
# par(mfrow = c(2,4))
# plot(richness.present)
# plot(richness.future)
# plot(richness.futureR)
# plot(richness.presfut)
# plot(Jaccad.local)
# plot(Jaccad.extinct)
# plot(Jaccad.migrate)

#
writeRaster(richness.present, paste0(SDM.path, 'richness_more/richnessP_', scenarios[id.s], '.tif'), overwrite=T)
writeRaster(richness.future,  paste0(SDM.path, 'richness_more/richnessF_', scenarios[id.s], '.tif'), overwrite=T)
writeRaster(richness.futureR, paste0(SDM.path, 'richness_more/richnessR_', scenarios[id.s], '.tif'), overwrite=T)
writeRaster(richness.presfut, paste0(SDM.path, 'richness_more/richness_presfut_', scenarios[id.s], '.tif'), overwrite=T)
writeRaster(Jaccad.local,     paste0(SDM.path, 'richness_more/Jaccad_local_',   scenarios[id.s], '.tif'), overwrite=T)
writeRaster(Jaccad.extinct,   paste0(SDM.path, 'richness_more/Jaccad_extinct_', scenarios[id.s], '.tif'), overwrite=T)
writeRaster(Jaccad.migrate,   paste0(SDM.path, 'richness_more/Jaccad_migrate_', scenarios[id.s], '.tif'), overwrite=T)