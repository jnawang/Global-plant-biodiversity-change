library('Matrix')
library('tidyverse')
library(terra)
rm(list = ls())
gc()
#
# SDM.path <- 'C:/Users/jnawang/Box/Biodiversity/SDM/'
SDM.path <- '../../SDM072023/'
template <- terra::rast(paste0(SDM.path, 'preds/bios_pres/bio1.tif'))

id.cells    <- which(!is.na(terra::values(template)))
n.cells.eff <- length(id.cells)                                      # 1759642  total number of effective cells

command_args <- commandArgs(trailingOnly = TRUE)
id.chunk <- as.numeric(command_args[1])    # this should be received from argument variable

max.cell.chunk <- 6400                     # I need about 275 cores
icell.start    <- (id.chunk-1)*max.cell.chunk + 1
icell.end      <- min(id.chunk*max.cell.chunk, n.cells.eff)
ncell          <- icell.end - icell.start + 1
cells          <- id.cells[icell.start:icell.end]
#  
xy             <- terra::xyFromCell(template, cells)

# read the raster for all species
# present
path.present    <- paste0(SDM.path, 'Results/RangePres/')

# I only consider one scenario so far
# future
path.future     <- paste0(SDM.path, 'Results/ESM_370_2081-2100/RangeFut/')

# future realized
path.futureR    <- paste0(SDM.path, 'Realized/ESM_370_2081-2100/')

# find the well-modeled species
species.SDMgood  <- read.csv(paste0(SDM.path, 'merged_stats/SDM_perform_merged.csv'))
species.realized <- read.csv(paste0(SDM.path, 'merged_stats/area_shrink_merged.csv'))

species <- intersect(species.SDMgood$sps, species.realized$Species)
nsps    <- length(species)
print(nsps)

# create a Large logical matrix to store all the distribution data    # I need to test if this works 
getClass('lMatrix')
dist.present <- Matrix(FALSE, nrow = ncell, ncol = nsps, sparse = T)
dist.future  <- Matrix(FALSE, nrow = ncell, ncol = nsps, sparse = T)
dist.futureR <- Matrix(FALSE, nrow = ncell, ncol = nsps, sparse = T)


t1 <- Sys.time()
print('Before reading any maps')
print(t1)
for (i.sps in 1:nsps) {
  # get distribution value of certain chunk of cells [10 thousand]
  # i.sps <- 47      # testing lines
  if (i.sps %% 1000 == 0) {print(i.sps)}
  present.map <- terra::rast(paste0(path.present, species[i.sps], '.tif'))
  future.map  <- terra::rast(paste0(path.future,  species[i.sps], '.tif'))
  futureR.map <- terra::rast(paste0(path.futureR, species[i.sps], '.tif'))
  
  # I need to coerce 0 and 1 into FALSE and TRUE, and put it into a sparse matrix 
  # I want to know how large the matrix is.
  tmp <- terra::extract(present.map, xy)[,1]==1
  if (sum(tmp, na.rm = T) > 0) {
    dist.present[which(tmp), i.sps] <- TRUE      # something wrong here; should check if value == 1
  }
  tmp <- terra::extract(future.map, xy)[,1]==1
  if (sum(tmp, na.rm = T) > 0) {
    dist.future[which(tmp), i.sps]  <- TRUE
  }
  tmp <- terra::extract(futureR.map, xy)[,1]==1
  if (sum(tmp, na.rm = T) > 0) {
    dist.futureR[which(tmp), i.sps]  <- TRUE
  }
}
t2 <- Sys.time()
print('After reading all maps')

print(t2-t1)
print(object.size(dist.present))
print(object.size(dist.future))
print(object.size(dist.futureR))


# open two binary files
output.present = file(paste0('map2bin/', id.chunk, "present.dat"), "wb")
output.future  = file(paste0('map2bin/', id.chunk, "ESM_370_2081-2100.dat"), "wb")
output.futureR = file(paste0('map2bin/', id.chunk, "Realized_ESM_370_2081-2100.dat"), "wb")

# for three files
writeBin(as.integer(ncell), output.present)
writeBin(as.integer(nsps), output.present)
writeBin(as.integer(ncell), output.future)
writeBin(as.integer(nsps), output.future)
writeBin(as.integer(ncell), output.futureR)
writeBin(as.integer(nsps), output.futureR)

# try to save dist.present and dist.future to binary files; every row we record number of non-zero numbers in each other and the column number
for (icell in 1:ncell) {
  # print("==============")
  # print(cells[icell])
  # print(icell)

  # here is what we need to read
  # write out three files
  # the first binary file
  sps.id <- which(dist.present[icell, ])
  writeBin(cells[icell], output.present)
  #print(length(sps.id))
  writeBin(length(sps.id), output.present)
  if(length(sps.id)>0){
    #print(sps.id)
    writeBin(sps.id, output.present)
  }

  # the second binary file
  sps.id <- which(dist.future[icell, ])
  writeBin(cells[icell], output.future)
  # print(length(sps.id))
  writeBin(length(sps.id), output.future)
  if(length(sps.id)>0){
    # print(sps.id)
    writeBin(sps.id, output.future)
  }

  # the third binary file
  sps.id <- which(dist.futureR[icell, ])
  writeBin(cells[icell], output.futureR)
  # print(length(sps.id))
  writeBin(length(sps.id), output.futureR)
  if(length(sps.id)>0){
    # print(sps.id)
    writeBin(sps.id, output.futureR)
  }
}

close(output.present)
close(output.future)
close(output.futureR)