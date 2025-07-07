################################Introduction of the script################################
##this script uses an efficient methods to determine if a grid is reachable
##Author: Junna Wang, 3/21/2025
#-----workflow of this script-----
# 1. read raster and plant trait files
# 2. create output folders
# 3. loop through each climate scenarios
# 4. for each scenario, get its predicted shift rates
# 5. read grid distance files (created by other scripts)
# 6. determine if each future grid is reachable
#########################
library(terra)
library(tidyverse)     
library(dismo)
library(pbapply)
library(data.table)
rm(list = ls())
gc()

##########################################################Functions##################################################
#----------------------------Function 1: divide global grids into chunks-------
grid_info <- function(r, block_size_row, block_size_col) {
  nrow_r <- nrow(r)
  ncol_r <- ncol(r)
  total_cells <- nrow_r * ncol_r
  #
  n_chunks_row <- ceiling(nrow_r / block_size_row)
  n_chunks_col <- ceiling(ncol_r / block_size_col)
  #
  row_indices <- rep(1:nrow_r, each=ncol_r)
  col_indices <- rep(1:ncol_r, nrow_r)
  #
  chunk_row_id <- (row_indices - 1) %/% block_size_row + 1
  chunk_col_id <- (col_indices - 1) %/% block_size_col + 1
  chunk_id <- (chunk_row_id - 1) * n_chunks_col + chunk_col_id
  #
  grid_info <- data.table(
    cell_id = 1:total_cells,
    row = row_indices,
    col = col_indices,
    chunk_id = chunk_id
  )
  
  # remove grids in waters
  landcells <- which(!is.na(terra::values(r)))
  
  grid_info <- grid_info[cell_id %in% landcells]
  #
  return(grid_info)
}
#----------------------------end of Function 1----------------------

#----------------------------Function 2: Find reachable grids-------
identify_reachable_grids <- function(idc.tbd, dt.idc.tbd.global, dt.idc.present.global, dt.distance.ele, dt.distance.lat, max.dispersal.dist, d.lat, d.ele) {
  # r.present: the present species range map
  # r.future: the future species range map
  # dt.distance: data table saving distance between origin grid and target grid
  # return grids in dt.distance that are reachable: idc.reachable;
  
  # browser()
  if (!is.na(d.ele)) {
    dt.distance.ele.subset <- dt.distance.ele[distance <= max.dispersal.dist & dif_ele <= d.ele]
  } else {
    dt.distance.ele.subset <- dt.distance.ele[distance <= max.dispersal.dist]
  }

  if (!is.na(d.lat)) {
    dt.distance.lat.subset <- dt.distance.lat[distance <= max.dispersal.dist & dif_lat <= d.lat]
  } else {
    dt.distance.lat.subset <- dt.distance.lat[distance <= max.dispersal.dist]
  }

  idc.reachable.global.ele <- c()
  if (nrow(dt.distance.ele.subset) > 0) {
    dt.distance.ele.subset <- dt.distance.ele.subset[dt.idc.tbd.global, nomatch = 0, on = "origin"]

    if (nrow(dt.distance.ele.subset) > 0) {
      dt.distance.ele.subset <- dt.distance.ele.subset[dt.idc.present.global, nomatch = 0, on = "target"]
    }

    if (nrow(dt.distance.ele.subset) > 0) {
      idc.reachable.global.ele <- dt.distance.ele.subset[, .SD[1], by=origin][, origin]
    }

    # a simpler way, but much slower!
    # idc.reachable.global.ele <- dt.distance.ele[origin %in% idc.tbd.global, .(flag = any(
    #     target %in% idc.present.global
    #   )), by = origin][flag == TRUE, origin]

  }

  idc.reachable.global.lat <- c()
  if (nrow(dt.distance.lat.subset) > 0) {
    dt.distance.lat.subset <- dt.distance.lat.subset[dt.idc.tbd.global, nomatch = 0, on = "origin"]

    if (nrow(dt.distance.lat.subset) > 0) {
      dt.distance.lat.subset <- dt.distance.lat.subset[dt.idc.present.global, nomatch = 0, on = "target"]
    }

    if (nrow(dt.distance.lat.subset) > 0) {
      idc.reachable.global.lat <- dt.distance.lat.subset[, .SD[1], by=origin][, origin]
    }
  }

  # get reachable global ID.
  idc.reachable.global <- c(idc.reachable.global.ele, idc.reachable.global.lat)

  # get reachable local ID.
  if (length(idc.reachable.global) > 0) {
    idc.reachable <- idc.tbd[match(idc.reachable.global, dt.idc.tbd.global$origin)]
  } else {
    idc.reachable <- c()
  }
  return(idc.reachable)
}
###########################old code, work but not efficient for large dataset#############################  
#   # select based on max dispersal distance
#   dt.distance.subset <- dt.distance[distance <= max.dispersal.dist]
#   dt.distance.subset <- dt.distance.subset[target %in% idc.present.global]
#   dt.distance.subset <- dt.distance.subset[origin %in% idc.tbd.global]
# 
#   # the tbd cells satisfy max dispersal distance
#   idc.tbd.global.subset <- unique(dt.distance.subset$origin)
# 
#   if (length(idc.tbd.global.subset) > 0) {
#     dt.mountain <- data.table(origin=idc.tbd.global.subset, mountain=terra::values(r.mountain)[idc.tbd.global.subset])
#     dt.distance.subset <- merge(dt.distance.subset, dt.mountain, by='origin', all.x = TRUE)
#     #
#     dt.distance.subset.ele <- subset(dt.distance.subset, mountain==1)
#     if (!is.na(d.ele) & nrow(dt.distance.subset.ele) > 0) {
#       ele.origin <- terra::extract(r.elev, dt.distance.subset.ele$origin)
#       ele.target <- terra::extract(r.elev, dt.distance.subset.ele$target)
#       dt.distance.subset.ele$dif_ele <- ele.origin - ele.target
#       # select elevational reachable
#       dt.distance.subset.ele <- dt.distance.subset.ele[dif_ele <= d.ele]
#     }
# 
#     dt.distance.subset.lat <- subset(dt.distance.subset, is.na(mountain))
#     if (!is.na(d.lat) & nrow(dt.distance.subset.lat) > 0) {
#       lat.origin <- terra::yFromCell(r.template, dt.distance.subset.lat$origin)
#       lat.target <- terra::yFromCell(r.template, dt.distance.subset.lat$target)
#       dt.distance.subset.lat$dif_lat <- lat.origin - lat.target
#       # select elevational reachable
#       dt.distance.subset.lat <- dt.distance.subset.lat[dif_lat <= d.lat]
#     }
# 
#     # get reachable global ID.
#     idc.reachable.global <- unique(c(dt.distance.subset.ele$origin, dt.distance.subset.lat$origin))
# 
#     # get reachable local ID.
#     idc.reachable <- idc.tbd[match(idc.reachable.global, idc.tbd.global)]
# 
#     return(idc.reachable)
#   } else {
#     # no reachable grids
#     return(c())
#   }
# }
###########################end of old code, work but not efficient for large dataset#############################

#----------------------------end of Function 2---------------------

##################################################################################################################################
# First of all, define a universal shift rates
shift.rate.lat.universal <- 0.7          # 0.7 km/year; unit is km/year
shift.rate.ele.universal <- 3.7        # 3.7 m/year 

##################################################################################################################################
###### Read files shared by all species
# plant name and traits
# median speed
shift.rates.plant <- read.csv('predict_shift_rate2025.csv')
plants            <- read.csv('plant_taxonomy_nloc_gf_mage_disp.csv')
#
# global map files: mountain, elevation, and global barrier
r.mountain  <- terra::rast('k3binary5m.tif')
r.elev      <- terra::rast('elev.tif')

# read a template to index global ID of new cells, newly added. 
r.template <- terra::rast('bio1.tif')

###### get block information ######
block_size_row = 20
block_size_col = 50
chunk_grid_mapping <- grid_info(r.template, block_size_row, block_size_col)
chunk_index <- sort(unique(chunk_grid_mapping$chunk_id))

cell_chunk_rows <- read.csv("cell_chunk_rows.csv")
cell_chunk_rows$chunk_id <- chunk_index[cell_chunk_rows$chunk_id]

chunk_grid_mapping <- merge(chunk_grid_mapping, cell_chunk_rows, by=c("cell_id", "chunk_id"), all.x = TRUE)
setkey(chunk_grid_mapping, cell_id)

###### get species distribution model information ######
# on my laptopt, 1 SDM
SDMs <- c("SDM")
SDM.path.output  <- file.path("/Users/junnawang/UCDLab/Biodiversity", SDMs)

# 3 SDMs, which are stored in the folders
# SDMs <- c("SDM2024MAXE_RFDS", "SDM2024MAXE", "SDM2024RFDS")
# SDM.path.output <- file.path('..', SDMs)

nSDMs <- length(SDMs)

###### set species' present range folder and species' future range folder ######
# remember to give one folder name for each SDM, because SDMs may use different thresholding method
# on my laptopt, 1 SDM
RangePres <- c("RangePres")
RangeFut <- c("RangeFut")

# # on cluster, 3 SDMs, should be consistent with the order of SDMs
# RangePres <- c("RangePres95s", "RangePres", "RangePres95s")
# RangeFut <- c("RangeFut95s", "RangeFut", "RangeFut95s")

###### get earth system models and CO2 emission information ######
# 16 future scenarios;
gcms  <- c("ACCESS", "ESM", "MIROC6", "MPIESM")   # DO NOT CHANGE THE ORDERS OF GCMS!
ssps  <- c(126, 245, 370, 585)
# times <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
# years <- c(rep(45,length(ssps)), rep(65,length(ssps)), rep(85,length(ssps)), rep(105,length(ssps)))
########only consider one time frame#######################
times <- c("2081-2100")
years <- c(rep(95, length(ssps)*length(gcms)))
futvars   <- expand.grid(ssps, gcms, times)
scenarios <- paste(futvars[,2],futvars[,1],futvars[,3],sep="_")    
###DOUBLE CHECK AND MAKE SURE SCENARIOS: ACCESS_126, ACCESS_245, ACCESS_370, ACCESS_585; ESM_126, ESM_245, ESM_370, ESM_585; 
###MIROC6_126, MIROC6_245, MIROC6_370, MIROC6_585; MPIESM_126, MPIESM_245, MPIESM_370, MPIESM_585;
nscenarios <- length(scenarios)

###### get 4 migration velocity information ######
migrations <- c('Realized975', 'Realized025', 'Realized', 'RealizedUniversal') # DO NOT CHANGE THE ORDERS OF MIGRATIONS!
nshift.rates <- length(migrations)

##################################################################################################################################
##################################################################################################################################
###create a function for a species
get.realized.map <- function(name) {
  # print species name first
  print(name)
  name <- 'Rhagodia preissii'          # 'Cornus florida'    # 'Eriophorum angustifolium'  'Oreopanax flaccidus'

  ############################################Read species dispersal characteristics##############################################
  plant <- plants[plants$sps==name, ]
  if (plant$disp == 'Autochor' | plant$disp == 'Ombrochor') {
    max.dispersal <- 16    # unit = m; max dispersal as 
  } else {
    if (plant$disp == 'Zoochor' | plant$disp == 'Hemerochor') {
      max.dispersal <- 5740
    } else {
      max.dispersal <- 4200
    }
  }
  
  # assign dispersal capacities for each scenario
  d.lat <- matrix(data=NA, nrow=nscenarios, ncol=nshift.rates)
  d.ele <- matrix(data=NA, nrow=nscenarios, ncol=nshift.rates)
  max.dispersal.dist <- matrix(data=NA, nrow=nscenarios, ncol=nshift.rates)
  #
  shift.rates <- shift.rates.plant[shift.rates.plant$Species==name, ]
  for (iscenario in 1:nscenarios) {
    for (ishift.rate in 1:nshift.rates) {
      if (ishift.rate <= 3) {
        # use predict shift rates
        if (nrow(shift.rates) == 0) {
          d.lat[iscenario, ishift.rate] <- NA                       # no prediction
          d.ele[iscenario, ishift.rate] <- NA
        } else {
          if (ishift.rate == 1) {           # fast velocity
            d.lat[iscenario, ishift.rate] <- shift.rates$Lat975[shift.rates$scenario==iscenario] * years[iscenario] * 1000  # unit is meter
            d.ele[iscenario, ishift.rate] <- shift.rates$Ele975[shift.rates$scenario==iscenario] * years[iscenario]         # unit m           
          } else if (ishift.rate == 2) {    # slow velocity
            d.lat[iscenario, ishift.rate] <- shift.rates$Lat025[shift.rates$scenario==iscenario] * years[iscenario] * 1000  # unit is meter
            d.ele[iscenario, ishift.rate] <- shift.rates$Ele025[shift.rates$scenario==iscenario] * years[iscenario]         # unit m                     
          } else if (ishift.rate == 3) {    # medium velocity
            d.lat[iscenario, ishift.rate] <- shift.rates$Lat50[shift.rates$scenario==iscenario] * years[iscenario] * 1000  # unit is meter
            d.ele[iscenario, ishift.rate] <- shift.rates$Ele50[shift.rates$scenario==iscenario] * years[iscenario]         # unit m                     
          }
        }
        
        # maximum dispersal range
        max.dispersal.dist[iscenario, ishift.rate] <- max.dispersal*years[iscenario]/plant$AgeOfMaturity_years
        if (!is.na(d.lat[iscenario, ishift.rate])) {
          max.dispersal.dist[iscenario, ishift.rate] <- max(max.dispersal.dist[iscenario, ishift.rate], d.lat[iscenario, ishift.rate]) 
        }
        
      } else if (ishift.rate == 4) {
        # use universal shift rates
        d.lat[iscenario, ishift.rate] <- shift.rate.lat.universal * years[iscenario] * 1000
        d.ele[iscenario, ishift.rate] <- shift.rate.ele.universal * years[iscenario]
        
        # make maximum dispersal range = universal shift rates
        max.dispersal.dist[iscenario, ishift.rate] <- max(d.lat[iscenario, ishift.rate], 0.0)
      }
    }
  }

  # browser()
  ############################################################Create folders############################################
  folders <- c(migrations, 'area_shrink', "realized_map")
  for (iSDM in 1:nSDMs) {
    for (ifolder in 1:length(folders)) {
      if(!dir.exists(file.path(SDM.path.output[iSDM], folders[ifolder]))) {
        dir.create(file.path(SDM.path.output[iSDM], folders[ifolder]), recursive = TRUE)
      }
      #
      if (ifolder > nshift.rates) {next}
      subfolders <- scenarios
      for (subfolder in subfolders) {
        if(!dir.exists(file.path(SDM.path.output[iSDM], folders[ifolder], subfolder))) {
          dir.create(file.path(SDM.path.output[iSDM], folders[ifolder], subfolder), recursive = TRUE)
        }
      }
    } 
  }
  
  ##############################################read in current and future range rasters############################################
  # nSDMs r.present
  # nSDMs * nscenarios r.future
  r.present.all <- vector("list", nSDMs)
  r.future.all  <- vector("list", nSDMs * nscenarios)
  
  # initialize potential future cell id list
  # nSDMs * nscenarios elements in a vector list
  idc.tbd.all <- vector("list", nSDMs * nscenarios)
  idc.tbd.all <- lapply(idc.tbd.all, function(x) numeric(0))
  
  # nSDMs * nscenarios elements in a data table list
  idc.tbd.global.all <- vector("list", nSDMs * nscenarios)
#  idc.tbd.global.all <- lapply(idc.tbd.global.all, function(x) numeric(0))

  # nSDMs elements in a data table list
  idc.present.global.all <- vector("list", nSDMs)
#  idc.present.global.all <- lapply(idc.present.global.all, function(x) numeric(0))
  
  for (iSDM in 1:nSDMs) {
    print(paste0('SDM: ', SDMs[iSDM]))
    
    # read present range
    r.present <- terra::rast(file.path(SDM.path.output[iSDM],'Results', RangePres[iSDM], paste0(name, '.tif')))
    r.present[r.present == 0] <- NA
    
    r.present.all[[iSDM]] <- r.present
    idc.present <- which(!is.na(terra::values(r.present)))
    xy <- xyFromCell(r.present, idc.present)
    
    # create present range grid global ID data table and put it into a list
    dt.idc.present.global <- data.table(target = cellFromXY(r.template, xy))
    setkey(dt.idc.present.global, target)
    idc.present.global.all[[iSDM]] <- dt.idc.present.global
    
    for (iscenario in 1:nscenarios) {
      print(paste0('Scenario:', iscenario))
      
      # read future range
      r.future  <- terra::rast(file.path(SDM.path.output[iSDM], 'Results', scenarios[iscenario], RangeFut[iSDM], paste0(name, '.tif')))   #Asplenium bulbiferum    # Cornus florida  
      r.future[r.future == 0] <- NA    
      
      ifuture <- (iSDM - 1) * nscenarios + iscenario
      r.future.all[[ifuture]] <- r.future
      
      # get potential future grid and put it into a list
      idc.future <- which(!is.na(terra::values(r.future)))
      idc.tbd <- setdiff(idc.future, idc.present)
      idc.tbd.all[[ifuture]] <- idc.tbd     
      
      # create potential future grid global id table and put it into a list
      xy <- xyFromCell(r.present, idc.tbd)      
      dt.idc.tbd.global <- data.table(origin = cellFromXY(r.template, xy))  
      setkey(dt.idc.tbd.global, origin)
      idc.tbd.global.all[[ifuture]] <- dt.idc.tbd.global
    }
  }
  
  # browser()
  
  # get the max tbd grid index across all scenarios, and put it into a data table
  dt.idc.tbd.global.max <- unique(data.table::rbindlist(idc.tbd.global.all))
  setkey(dt.idc.tbd.global.max, origin)
  
  # get the max present grid index across all scenarios, and put it into a data table
  dt.idc.present.global.max <- unique(data.table::rbindlist(idc.present.global.all))
  setkey(dt.idc.present.global.max, target)
  
  # the max maximum dispersal distance across all scenarios
  max.dispersal.dist.max <- max(max.dispersal.dist)
  
  # initialize reachable cell id list
  # nSDMs * nscenarios * nshift.rates r.future.real
  idc.reachable.all <- vector("list", nSDMs * nscenarios * nshift.rates)
  idc.reachable.all <- lapply(idc.reachable.all, function(x) numeric(0))
  
  if (nrow(dt.idc.tbd.global.max) > 0) {

    # find all the distance files that will be used
    chunk_grid_mapping_subset <- chunk_grid_mapping[dt.idc.tbd.global.max, nomatch = 0]
    chunk_grid_files <- unique(chunk_grid_mapping_subset$chunk_id)
    chunk_grid_files <- match(chunk_grid_files, chunk_index)
    print(length(chunk_grid_files))
    
    # browser()
    
    dt.distance <- data.table(origin=integer(), target=integer(), distance=double())
    time0 <- Sys.time()
    
    for (ichunk in 1:length(chunk_grid_files)) {
      
      # read in distance files dt.distance and determine if a cell is reachable      
      system.time(dt.distance0 <- readRDS(file.path("/Users/junnawang/Downloads/chunks_select/", paste0(chunk_grid_files[ichunk], '.rds'))))
      
      
      
      dt.distance03 <- read_fst('/Users/junnawang/Downloads/2183.fts', as.data.table = TRUE)
      
      dt.distance0 <- read_fst('/Users/junnawang/Downloads/880.fts', as.data.table = TRUE)
      chunk_grid_mapping_subset_1chunk <- chunk_grid_mapping[chunk_id==chunk_index[880]]
      chunk_grid_mapping_subset_1chunk <- chunk_grid_mapping_subset_1chunk[!is.na(first_appearance)]
      
      dt.distance0 <- read_fst(file.path('..', 'grid_distance_earth', 'chunks2', paste0(chunk_grid_files[ichunk], '.fts')), as.data.table = TRUE)

      # select useful rows in this chunk file
      
      irows <- apply(chunk_grid_mapping_subset_1chunk[, c("first_appearance", "last_appearance")], 1, function(x) { return(seq(from=x[1], to=x[2])) })
      irows <- unlist(irows)
      dt.distance0 <- dt.distance0[irows]
      dt.distance0$origin <- irows
      print(head(irows))

      
      
      
      
      
      
            
      # time00 <- Sys.time()
#      subset useful dt.distance only; the order of subsetting strongly affect computational velocity!!!
      # dt.distance0 <- dt.distance0[origin %in% idc.tbd.global.max & 
      #                                target %in% idc.present.global.max &
      #                                distance <= max.dispersal.dist.max]
      # browser()
      # use three chains for the three filters, which is slower, because dt.distance0 is not keyed 
      dt.distance0 <- dt.distance0[distance <= max.dispersal.dist.max]
      
      setkey(dt.distance0, target)
      dt.distance0 <- dt.distance0[dt.idc.present.global.max, nomatch = 0]
      
      setkey(dt.distance0, origin)
      dt.distance0 <- dt.distance0[dt.idc.tbd.global.max, nomatch=0] 

      # setindex(dt.distance0, origin)
      # dt.distance0 <- dt.distance0[.(origin.unique),
      #                              nomatch=0,
      #                              on=.(origin)]
      # 
      # setindex(dt.distance0, target)
      # dt.distance0 <- dt.distance0[.(target.unique),
      #                              nomatch=0,
      #                              on=.(target)]
      
      # dt.distance0 <- dt.distance0[target %in% target.unique]      
      # dt.distance0 <- dt.distance0[origin %in% origin.unique]

      dt.distance <- rbind(dt.distance, dt.distance0)
      
      # time01 <- Sys.time()
      # print(paste0('subset', time01-time00))
      
      # wait until there are 3 million elements in dt.distance
      if (nrow(dt.distance) < 3000000 & ichunk < length(chunk_grid_files)) { next }
      
      # differentiate mountain and non-mountain grids
      dt.distance[, mountain := terra::values(r.mountain)[origin]]
      
      # subset mountain and non-mountain grids
      dt.distance.ele <- dt.distance[mountain==1]
      
      # add one column dif_ele to dt.distance.ele
      ele.origin <- terra::extract(r.elev, dt.distance.ele$origin)
      ele.target <- terra::extract(r.elev, dt.distance.ele$target)
      dt.distance.ele[, dif_ele := ele.origin - ele.target]
      
      # subset non-mountain grids       
      dt.distance.lat <- dt.distance[is.na(mountain)]

      lat.origin <- terra::yFromCell(r.template, dt.distance.lat$origin)
      lat.target <- terra::yFromCell(r.template, dt.distance.lat$target)
      
      # differentiate northern and southern hemisphere when calculate latitudinal differences
      
     # browser()
      
      # add one column dif_lat, here need to consider different ways of calculating lat difference between south and north hemispheres.
      dt.distance.lat[, dif_lat := sign(lat.origin) * (lat.origin - lat.target)]
      
      time1 <- Sys.time()
      print(paste0('Time needed to accumulate 3 million useful grid distance: ', time1 - time0))
      
      # browser()
      # start the computation
      for (iSDM in 1:nSDMs) {
        
        for (iscenario in 1:nscenarios) {
          ifuture <- (iSDM - 1) * nscenarios + iscenario
          
          # if no future grids outside the current range (i.e., tbd grids), skip this loop
          if (nrow(idc.tbd.global.all[[ifuture]]) == 0) { next }
          
          # should be removed because it takes time to do the cross-section
          # # if no future tbd grids in this dt.distance file, skip this loop
          # if (sum( unique(dt.distance$origin) %in% idc.tbd.global.all[[ifuture]]) == 0) { next }
          
          for (ishift.rate in 1:4) {
            
            # calculate index ireal: the index of future realized range
            ireal <- (iSDM - 1) * nscenarios * nshift.rates + (iscenario - 1) * nshift.rates + ishift.rate
            
            # call reachable function to get local id of grid that is reachable
            idc.reachable <- identify_reachable_grids(idc.tbd.all[[ifuture]], 
                                                      idc.tbd.global.all[[ifuture]], 
                                                      idc.present.global.all[[iSDM]], 
                                                      dt.distance.ele, dt.distance.lat, 
                                                      max.dispersal.dist[iscenario, ishift.rate], 
                                                      d.lat[iscenario, ishift.rate], 
                                                      d.ele[iscenario, ishift.rate])
            
            # merge with reachable id from prior files
            if (length(idc.reachable) > 0) {
              idc.reachable.all[[ireal]] <- c(idc.reachable.all[[ireal]], idc.reachable)
            }
          }
        }
      } 
      time2 <- Sys.time()
      print(time2 - time1)
      print(ichunk)
      
      # reset the data table for next-round reading of grid distance 
      dt.distance <- data.table(origin=integer(), target=integer(), distance=double())
      time0 <- Sys.time()
    }
  }

  # browser()      
  # get r.future.real files
  for (iSDM in 1:nSDMs) {
    r.present <- r.present.all[[iSDM]]
    
    # initialize realized maps
    r.future.real.all <- terra::rast()
    
    # initialize the area shrink data frame
    area.shrink <- data.frame(matrix(0, nrow=1, ncol=1+nscenarios*(2+nshift.rates)))
    colnames(area.shrink) <- c('presA', paste0('suitAF', 1:16), paste0('sameAF', 1:16), paste0('realAF_fast', 1:16), paste0('realAF_slow', 1:16), paste0('realAF_med', 1:16), paste0('realAF_univ', 1:16))  # Present area, future suitable area, future realized area. 
    area.shrink[1,1]      <- sum(terra::values(terra::cellSize(r.present))[!is.na(terra::values(r.present))])
    
    for (iscenario in 1:nscenarios) {
      
      ifuture <- (iSDM - 1) * nscenarios + iscenario
      r.future <- r.future.all[[ifuture]]
      
      for (ishift.rate in 1:4) {
        
        # calculate index ireal: the index of future realized range
        ireal <- (iSDM - 1) * nscenarios * nshift.rates + (iscenario - 1) * nshift.rates + ishift.rate
        
        # use the reachable grid results
        r.future.real <- r.future
        idc.unreachable <- setdiff(idc.tbd.all[[ifuture]], idc.reachable.all[[ireal]])
        
        # set unreachable grid to NA
        if (length(idc.unreachable) > 0) {
          r.future.real[idc.unreachable] <- NA
        }
        
        # calculate realized future habitat area, and put it to area.shrink file
        r.final <- r.future.real
        terra::values(r.final)[intersect(which(!is.na(terra::values(r.present))), which(!is.na(terra::values(r.future.real))))] <- 0
        terra::values(r.final)[setdiff(which(!is.na(terra::values(r.present))), which(!is.na(terra::values(r.future.real))))] <- -1
        #    plot(r.final)
        
        # browser()
        
        # calculate shrink areas
        if (ishift.rate == 1) {
          area.shrink[1,1+iscenario]            <- sum(terra::values(terra::cellSize(r.future))[!is.na(terra::values(r.future))]) / area.shrink[1,1]
          area.shrink[1,1+iscenario+nscenarios] <- sum(terra::values(terra::cellSize(r.final))[which(terra::values(r.final)==0)]) / area.shrink[1,1]
        }
        area.shrink[1,1+iscenario+nscenarios*(ishift.rate+1)] <- sum(terra::values(terra::cellSize(r.future.real))[!is.na(terra::values(r.future.real))]) / area.shrink[1,1]
        
        r.future.real.all <- c(r.future.real.all, r.future.real)
      }
    }
    
#    write out realized future range
    writeRaster(r.future.real.all, file.path(SDM.path.output[iSDM], 'realized_map', paste0(name, '.tif')), overwrite=TRUE)
    
    # write our area.shrink for each SDM
    write.csv(area.shrink, file.path(SDM.path.output[iSDM], 'area_shrink', paste0(name, '.csv')), row.names = F)
  }      

}
#----------------------------end of function 3--------------------------------

#####################################################################################################################################
###########################################################start to call the function################################################
####get modelled species name; 
species <- read.csv('SDM_perform_merged_good2025.csv')
species <- species$sps

# read species ID from command
command_args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(command_args[1])
# id <- id + 40000

# one command line calculate 8 species: ((id-1)*8+1):min((id*8),nsp)   
pblapply(species[id], get.realized.map)

