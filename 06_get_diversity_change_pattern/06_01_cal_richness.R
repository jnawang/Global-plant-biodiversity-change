# calculate species richness distribution based on SDM results
# change this code to do parallel for different scenarios
 
library(terra, sp)
#
command_args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(command_args[1])      # can vary from 0 to 32; 0: present; 1:16: future_suit; 17:32: future_realized; 33:48: future realized975; 49:60: future realized025;
#
# corresponding scenarios
if (id > 0) {
  id.s <- (id - 1) %% 16 + 1
}
#
gcms  <- c("ACCESS", "ESM", "MIROC6", "MPIESM")
ssps  <- c(126, 245, 370, 585)
# times <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
times <- c("2081-2100")
ngcm  <- length(gcms)
nssp  <- length(ssps)
ntime <- length(times)
scenarios <- c()
for (igcm in 1:ngcm) {
  for (itime in 1:ntime) {
    for (issp in 1:nssp) {
      scenarios[(igcm-1)*ntime*nssp+(itime-1)*nssp+issp] <- paste(gcms[igcm], ssps[issp], times[itime], sep = '_')
    }
  }
}

#
SDM.path <- '../../SDM2024MAXE_RFDS/'
# SDM.path <- 'C:/Users/jnawang/Box/Biodiversity/SDM/'   # TEST LINE

# find the file path
if (id==0) {
  path <- paste0(SDM.path, 'Results/RangePres95s/')
} else if (id >= 1 & id <= 16) {
  path <- paste0(SDM.path, 'Results/', scenarios[id.s], '/RangeFut95s/')
  print(path)
} else if (id >= 17 & id <= 32) {
  path <- paste0(SDM.path, 'Realized/', scenarios[id.s], '/')
  print(path)
} else if (id >= 33 & id <= 48) {
  path <- paste0(SDM.path, 'Realized975/', scenarios[id.s], '/')
  print(path)
} else if (id >= 49 & id <= 64) {
  path <- paste0(SDM.path, 'Realized025/', scenarios[id.s], '/')
  print(path)
} else if (id >= 65 & id <= 80) {
  path <- paste0(SDM.path, 'RealizedUniversal/', scenarios[id.s], '/')
  print(path)
} else if (id >= 81 & id <= 112) {
  path.pres <- paste0(SDM.path, 'Results/RangePres95s/')
  path.fut  <- paste0(SDM.path, 'Results/', scenarios[id.s], '/RangeFut95s/')
}

# find the well-modeled species
species.SDMgood <- read.csv('../SDM_perform_merged_good2025.csv')
species <- species.SDMgood$sps
#
nsps    <- length(species)
#
richness <- terra::rast('../../SDM2024MAXE/preds/bios_pres/bio1.tif')
richness[!is.na(richness)] <- 0
#
if (id >= 81 & id <= 112) {
  for (i in 1:nsps) {
    if (i %% 1000 == 0) {print(i)}
    map.pres <- terra::rast(paste0(path.pres, species[i], '.tif'))
    map.fut  <- terra::rast(paste0(path.fut, species[i], '.tif'))
    #
    if (id >= 81 & id <= 96) {
      icell.sps <- union(which(terra::values(map.pres)==1), which(terra::values(map.fut)==1))
    } else if (id >= 97 & id <= 112) {
      icell.sps <- intersect(which(terra::values(map.pres)==1), which(terra::values(map.fut)==1))
    }
    xy  <- terra::xyFromCell(map.pres, icell.sps)
    icell <- terra::cellFromXY(richness, xy)
    #
    richness[icell] <- richness[icell] + 1
  }
} else {
  for (i in 1:nsps) {
    if (i %% 1000 == 0) {print(i)}
    map <- terra::rast(paste0(path, species[i], '.tif'))
    xy <- terra::xyFromCell(map, which(terra::values(map)==1))
    icell <- terra::cellFromXY(richness, xy)
    #
    richness[icell] <- richness[icell] + 1
  }
}

#
if (id==0) {
  writeRaster(richness, paste0(SDM.path, 'richness/richnessP.tif'), overwrite=T)
} else if (id >= 1 & id <= 16) {
  writeRaster(richness, paste0(SDM.path, 'richness/richnessF_', scenarios[id.s], '.tif'), overwrite=T)
} else if (id >= 17 & id <= 32) {
  writeRaster(richness, paste0(SDM.path, 'richness/richnessR_', scenarios[id.s], '.tif'), overwrite=T)
} else if (id >= 33 & id <= 48) {
  writeRaster(richness, paste0(SDM.path, 'richness/richnessR975_', scenarios[id.s], '.tif'), overwrite=T)
} else if (id >= 49 & id <= 64) {
  writeRaster(richness, paste0(SDM.path, 'richness/richnessR025_', scenarios[id.s], '.tif'), overwrite=T)
} else if (id >= 65 & id <= 80) {
  writeRaster(richness, paste0(SDM.path, 'richness/RealizedUniversal_', scenarios[id.s], '.tif'), overwrite=T)
} else if (id >= 81 & id <= 96) {
  writeRaster(richness, paste0(SDM.path, 'richness/richness_PresFut_', scenarios[id.s], '.tif'), overwrite=T)
} else if (id >= 97 & id <= 112) {
  writeRaster(richness, paste0(SDM.path, 'richness/richnessO_', scenarios[id.s], '.tif'), overwrite=T)
}
