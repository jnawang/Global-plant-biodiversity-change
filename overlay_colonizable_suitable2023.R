library(raster)
library(tidyverse)     
library(dismo)
library(pbapply)
rm(list = ls())
gc()
###############################################################################################################################################
# a new climate analogue function
`ensemble.analogue1` <- function(targetdata, out, normdata, weightdata, z=2)
{
  .BiodiversityR <- new.env()
  # predict distance to target
  for (i in 1:ncol(out)) {
    out[,i] <- as.numeric(out[,i]) - as.numeric(targetdata[i])
    out[,i] <- abs(out[,i])
    out[,i] <- as.numeric(out[,i]) / as.numeric(normdata[i])
    out[,i] <- (out[,i]) ^ z
    out[,i] <- as.numeric(out[,i]) * as.numeric(weightdata[i])
  }
  p  <- rowSums(out)
  z2 <- 1/z
  p <- p ^ z2
  p <- round(p, digits=8)
  return(p)
}
#######end of the environmental analogue function

# the second function to divide a raster into patches with certain distance
patchify <- function(x, distance) {
  require(rgeos)
  require(sp)
  require(sf)
  require(raster)
  #
  p4s = sp::proj4string(x)           # Hopefully, this is not used at all! 
  p <- raster::rasterToPolygons(x,fun=function(x){x>0},na.rm=TRUE,dissolve=TRUE)         # this step takes times; does n = 8 or 16 matter?
  p <- sp::disaggregate(as(sf::st_as_sf(p), 'Spatial'))
  pproj <- sp::spTransform(p, p4s)
  bproj <- rgeos::gBuffer(pproj, width=distance)
  pproj$patch <- sp::over(pproj, sp::disaggregate(bproj))
  b <- sp::spTransform(pproj, sp::CRS(p4s))                  # why have this step?
  pout <- raster::aggregate(b, by='patch')
  pout$patch <- base::as.factor(pout$patch)                  # polygon
  rout <- raster::rasterize(pout, x)                         # raster
  return(rout)
}
#######end of the patchy function

##################################################################################################################################
################################################Read files shared by all species##################################################
SDM.path <- '../SDM072023/'
#
# plant name and traits
# median speed
shift.rates.plant <- read.csv('Final_prediction_shift_rate_median.csv')
plants            <- read.csv('plant_taxonomy_nloc_gf_mage_disp.csv')
#
# global map files: mountain, elevation, and global barrier
r.elevation <- raster('elev.tif')
r.mountain  <- raster('k3binary5m.tif')
r.barrier   <- raster('global_barrier.tif')           # cell values of barriers are 3

# read a template to index global ID of new cells, newly added. 
r.template <- raster(paste0(SDM.path, "preds/bios_pres/bio1.tif"))

##################################################################################################################################
##################################################################################################################################
###create a function for a species
get.realized.map <- function(name) {
  # print species name first
  print(name)
  #######################################################Read Raster files##########################################################
  # read species name and current distribution
  # name <- 'Oreopanax flaccidus'          # 'Cornus florida'    # 'Eriophorum angustifolium'  'Oreopanax flaccidus'
  ## read species present range
  r.present <- raster(paste0(SDM.path,'Results/RangePres/', name, '.tif'))
  r.present[r.present == 0] <- NA
  
  ########crop mountain, elevation, and global barrier maps to species ranges
  r.ele   <- crop(r.elevation, r.present)           
  r.mount <- crop(r.mountain, r.present)            # make the extent smaller;
  r.bar   <- resample(r.barrier, r.present)
  
  #########################################Read relative importance of environmental variables######################################
  stats <- read.csv(paste0(SDM.path, 'Results/stats/', name, '_stats.csv'))
  stats <- data.frame(variable = stats[3:19,1], ri = as.numeric(stats[3:19,3]))
  stats <- stats[which(stats$ri>1.0), ]     # only consider the variables whose importance > 1 %
  stats$category <- 'bios'
  asoil <- list.files(paste0(SDM.path, "preds/soil/"), pattern = ".tif$")
  asoil <- gsub('.{4}$', '', asoil)
  stats$category[stats$variable %in% asoil] <- 'soil'
  
  #######read current environmental data
  stats$path <- ''
  for (i in 1:nrow(stats)) {
    if (stats$category[i] == 'soil') {
      stats$path[i] <- paste0(SDM.path, "preds/soil/", stats$variable[i], '.tif')
    } else {
      stats$path[i] <- paste0(SDM.path, "preds/bios_pres/", stats$variable[i], '.tif')
    }
  }
  
  current.stack <- lapply(stats$path, raster)
  current.stack <- stack(current.stack)
  current.stack <- stack(mask(crop(current.stack, r.present), r.present))
  
  ####get weight data and norm data
  weightdata <- stats$ri/sum(stats$ri)
  normdata   <- raster::cellStats(current.stack, stat="sd")
  #######################################################Read species' characteristics##############################################
  lifespan <- 1   # this can be improved further
  plant <- plants[plants$sps==name, ]
  if (plant$disp == 'Autochor' | plant$disp == 'Ombrochor') {
    self.dispersal <- TRUE
    max.dispersal <- 16    # unit = m; max dispersal as 
  } else {
    self.dispersal <- FALSE
    if (plant$disp == 'Zoochor' | plant$disp == 'Hemerochor') {
      max.dispersal <- 5740
    } else {
      max.dispersal <- 4200
    }
  }
  
  # self-defining parameters
  # how many points for a patch? 2 how many edges for a patch? 3 how many cells for a patch for patch separation
  max.patch.point <- 500
  max.edge.point  <- 200
  patch.dist      <- 50000         # maximum distance to define a patch (50 km)
  ############################################################Step 0#############################################################
  #################################################calculate max climate distance within one patch###############################
  ##########################################Step 0.1 divide current species range into different patches#########################
  # patch.present <- clump(r.present, directions=8, gaps=F)
  patch.present <- patchify(r.present, patch.dist)
  patch.freq.p  <- freq(patch.present)
  # use representative points to represent each patch
  patch.point.present <- patch.present
  ##########Step 0.2 get maximum distance within a patch using current map
  npatch.present <- max(raster::values(patch.present), na.rm = T)
  ncell.eff      <- sum(patch.freq.p[1:(npatch.present), 2])      # total effective cells.
  # loop through every patch
  dist.thrd <- 0
  for (ipatch in 1:npatch.present) {
    if (patch.freq.p[ipatch,2] >= ncell.eff/npatch.present) {
      print(ipatch)
      patch  <- patch.present
      patch[patch!=ipatch] <- NA              # single patch
      # patch boundary; something wrong with this. 
      boundary <- boundaries(patch, type='inner')
      xy.bd <- xyFromCell(boundary, which(raster::values(boundary)==1))
      nbd   <- nrow(xy.bd)
      if (nbd > max.edge.point) {
        id_bd    <- sample(1:nbd, max.edge.point, replace=FALSE)
      } else {
        id_bd   <- 1:nbd
      }
      # select random points within one patch
      patch.point <- patch                        # single patch points
      if (patch.freq.p[ipatch,2] > max.patch.point) {
        patch.point.id <- dismo::randomPoints(patch, n=max.patch.point, cellnumbers=T, excludep=T)
        patch.point.NAid <- setdiff(which(!is.na(raster::values(patch))), patch.point.id)
        patch.point[patch.point.NAid] <- NA
        patch.point.present[patch.point.NAid] <- NA
      }
      # extract target data
      targetdata  <- data.frame(raster::extract(current.stack, xy.bd[id_bd,]))    
      # extract background data
      bgdata <- raster::values(current.stack)[!is.na(raster::values(patch.point)), ]
      if (class(bgdata)[1]=="numeric") {
        if (length(weightdata) == 1) {      # only one important variable
          bgdata <- as.matrix(bgdata)
          analogue.sd <- pbapply(targetdata, 1, ensemble.analogue1, out=bgdata, normdata=normdata, weightdata=weightdata, z=2)
        } else {                        # only one bg point
          bgdata <- t(as.matrix(bgdata))
          analogue.sd <- pbapply(targetdata, 1, ensemble.analogue1, out=bgdata, normdata=normdata, weightdata=weightdata, z=2)
          analogue.sd <- t(as.matrix(analogue.sd))
        }
      } else {
        analogue.sd <- pbapply(targetdata, 1, ensemble.analogue1, out=bgdata, normdata=normdata, weightdata=weightdata, z=2)
      }
      dist.thrd <- max(max(analogue.sd,na.rm=T), dist.thrd)
      print(c(ipatch,dist.thrd))
    }
  }
  # location, patch ID, and values of the representative present patch
  xy.present  <- xyFromCell(r.present, which(!is.na(raster::values(patch.point.present))))            # current representative cells
  patch.id    <- raster::values(patch.point.present)[!is.na(raster::values(patch.point.present))]
  targetdata  <- data.frame(raster::extract(current.stack, xy.present))
  #########################################################sTART TO WORK ON SCENARIOS####################################
  # future scenarios;
  gcms  <- c("ESM")
  ssps  <- c(126, 245, 370, 585)
  times <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
  years <- c(rep(45,length(ssps)), rep(65,length(ssps)), rep(85,length(ssps)), rep(105,length(ssps)))
  futvars   <- expand.grid(gcms, ssps, times)
  scenarios <- paste(futvars[,1],futvars[,2],futvars[,3],sep="_")
  ############################################################Create folders############################################
  folders <- c(scenarios, 'shift_direction', 'area_shrink', 'shift_direction_cell')
  for (folder in folders) {
    if(!dir.exists(paste0(SDM.path, 'Realized/', folder))){
      dir.create(paste0(SDM.path, 'Realized/', folder), recursive = TRUE)
    }
  }
  ######################################################Prepare for two output data frames###############################
  # the first data frame with species name, future patch center, future_percent, source patch center, source_percent
  shift.direction <- data.frame(scenarios=character(), patch.x.f=double(), patch.y.f=double(), percent.f=double(), patch.x.s=double(), patch.y.s=double(), percent.s=double())
  # the second data frame with 1+16+16+16 columns
  area.shrink     <- data.frame(matrix(0, nrow=1, ncol=49))
  colnames(area.shrink) <- c('presA', paste0('suitAF', 1:16), paste0('realAF', 1:16), paste0('sameAF', 1:16))  # Present area, future suitable area, future realized area. 
  area.shrink[1,1]      <- sum(raster::values(area(r.present))[!is.na(raster::values(r.present))])
  # the third data frame with scenario, globalID, angle, averaged_angle
  shift.direction.cell <- data.frame(iscenario=character(), GID = integer(), distance = double(), angle=double(), anglePavg=double())
  #
  shift.rates <- shift.rates.plant[shift.rates.plant$Species==name, ]
  #
  iresult   <- 0.0      # for writing out results
  #
  # loop through every scenario
  for (iscenario in 13:length(scenarios)) {
    # print the start of each scenario
    print(paste0('Scenario:', iscenario))
    # 
    r.future  <- raster(paste0(SDM.path, 'Results/', scenarios[iscenario], '/RangeFut/', name, '.tif'))   #Asplenium bulbiferum    # Cornus florida  
    r.future[r.future == 0] <- NA    
    #
    # check if a species has 0 suitable area
    if (sum(!is.na(raster::values(r.future)))==0) {
      writeRaster(r.future, paste0(SDM.path, 'Realized/', scenarios[iscenario], '/', name, '.tif'), overwrite=TRUE)
      iresult                       <- iresult + 1
      shift.direction[iresult, 1]   <- scenarios[iscenario]
      shift.direction[iresult, 2:7] <- 0.0
      next
    }
    
    # future environments: update future bios path
    for (i in 1:nrow(stats)) {
      if (stats$category[i] == 'bios') {
        stats$path[i] <- paste0(SDM.path, "preds/bios_fut/", stats$variable[i], '_', scenarios[iscenario], '.tif')
      }
    }
    future.stack <- lapply(stats$path, raster)
    future.stack <- stack(future.stack)
    future.stack <- stack(mask(crop(future.stack, r.future), r.future))
    #
    # Assign predicted species migration to this scenario; double check the meaning of scenario?
    if (nrow(shift.rates) == 0) {
      d.lat <- NA                       # no prediction
      d.ele <- NA
    } else {
      d.lat <- shift.rates$LatShiftR[shift.rates$scenario==iscenario] * years[iscenario] * 1000  # unit is meter
      d.ele <- shift.rates$EleShiftR[shift.rates$scenario==iscenario] * years[iscenario]         # unit m    
    }
    #
    # maximum dispersal range
    max.dispersal.dist <- max.dispersal*years[iscenario]/plant$AgeOfMaturity_years
    if (!is.na(d.lat)) {
      max.dispersal.dist <- max(max.dispersal.dist, abs(d.lat))                                  # not sure if I should have the abs(d.lat) here.
    }
    ############################################################Step 1#############################################################
    #################################################Remove cells that's beyond species migration prediction#######################
    ####calculate distance between current representative cells and every future cell
    idc.future <- which(!is.na(raster::values(r.future)))
    bgdata <- raster::values(future.stack)[idc.future, ]
    if (class(bgdata)[1]=="numeric") {
      if (length(weightdata) == 1) {      # only one important variable
        bgdata <- as.matrix(bgdata)
        analogue.sd <- pbapply(targetdata, 1, ensemble.analogue1, out=bgdata, normdata=normdata, weightdata=weightdata, z=2)
      } else {                        # only one bg point
        bgdata <- t(as.matrix(bgdata))
        analogue.sd <- pbapply(targetdata, 1, ensemble.analogue1, out=bgdata, normdata=normdata, weightdata=weightdata, z=2)
        analogue.sd <- t(as.matrix(analogue.sd))
      }
    } else {
      analogue.sd <- pbapply(targetdata, 1, ensemble.analogue1, out=bgdata, normdata=normdata, weightdata=weightdata, z=2)
    }
    # please note that the index of the analogue.sd is not the index of the raster map. 
    #####this step takes a few minutes
    
    r.future.real <- r.future
    raster::values(r.future.real) <- NA
    # process one patch each time
    for (idp in 1:npatch.present) {
      #
      r.bar.patch <- r.bar
      # apply criteria 1: grid distance from the source cell < maxmum dispersal capacity
      # not sure if I want to use this line. 
      # if (self.dispersal) {
      #   r.self.dispersal <- merge(r.present, r.future)
      #   r.bar.patch[is.na(r.self.dispersal)] <- 3
      # }
      #
      raster::values(r.bar.patch)[which(raster::values(patch.present)==idp)] <- 1
      r.d <- gridDistance(r.bar.patch, origin=1, omit=3)                          # Origin: present species area; omit: barrier
      cells.noblock <- which(raster::values(r.d) <= max.dispersal.dist)
      
      # apply criteria 2: from this patch
      descendant <- lapply(which(patch.id==idp), function(x) idc.future[which(analogue.sd[,x] <= dist.thrd)])
      descendant <- unique(unlist(descendant))
      
      # combine the two criteria
      descendant <- intersect(descendant, cells.noblock)
      
      # differentiate mountain and non-mountain cells
      descendant.lat <- intersect(which(is.na(raster::values(r.mount))), descendant)
      descendant.ele <- setdiff(descendant, descendant.lat)
      
      # apply criteria 3: latitudinal and elevational distance < predicted shift distance [reachable based on y and prediction]
      if (length(descendant.lat) > 0) {
        if (is.na(d.lat)) {
          raster::values(r.future.real)[descendant.lat] <- 1
        } else {
          # lat range of this patch
          lat.range <- range(yFromCell(r.present, which(raster::values(patch.present)==idp)), na.rm = T)
          # expand the range based predicated shift distance
          # both north
          if (lat.range[1] >= 0) {
            if (d.lat >= 0) {
              lat.range[2] = lat.range[2] + d.lat
            } else {
              lat.range[1] = lat.range[1] + d.lat
            }
          } else if(lat.range[2] <= 0) {         # both south
            if (d.lat >= 0) {
              lat.range[1] = lat.range[1] - d.lat
            } else {
              lat.range[2] = lat.range[2] - d.lat
            }  
          } else {                    # cross equator, expand, we donot consider shrink (because shrink is controlled by habitat suit)
            if (d.lat > 0) {
              lat.range[2] = lat.range[2] + d.lat
              lat.range[1] = lat.range[1] - d.lat
            }
          }
          # to show if y or elevation within this range
          temp <- yFromCell(r.future, descendant.lat)
          raster::values(r.future.real)[descendant.lat[which(temp>=lat.range[1] & temp<=lat.range[2])]] <- 1
        }      
      }
      #
      if (length(descendant.ele) > 0) {
        if (is.na(d.ele)) {
          raster::values(r.future.real)[descendant.ele] <- 1
        } else {
          ele.range <- range(raster::values(r.ele)[which(raster::values(patch.present) == idp)], na.rm = T)
          if (d.ele >= 0) {
            ele.range[2] = ele.range[2] + d.ele   # expand upper limit
          } else {
            ele.range[1] = ele.range[1] + d.ele   # expand lower limit
          }
          temp <- raster::values(r.ele)[descendant.ele]
          raster::values(r.future.real)[descendant.ele[which(temp>=ele.range[1] & temp<=ele.range[2])]] <- 1
        }      
      }
    }
    
    # may be deleted in the future
    # par(mfrow = c(2,4))
    # plot(r.present)
    # plot(r.future)
    # plot(r.future.real)
    # plot(mask(r.future, r.mount))
    # plot(mask(r.future.real, r.mount))
    # plot(mask(r.future, r.mount, inverse = TRUE))
    # plot(mask(r.future.real, r.mount, inverse = TRUE))
    
    # output the future species range
    writeRaster(r.future.real, paste0(SDM.path, 'Realized/', scenarios[iscenario], '/', name, '.tif'), overwrite=TRUE)
    #
    # check if a species has 0 realized area
    if (sum(!is.na(raster::values(r.future.real)))==0) {
      iresult                     <- iresult + 1
      shift.direction[iresult, 1] <- scenarios[iscenario]
      shift.direction[iresult, 2:7] <- 0.0
      next
    }
    #
    ############################################################Step 3#############################################################
    #################################################migration direction between 1-3 dominant patches##############################
    # a merged map with present only range = -1, present and future range = 0, and future only range = 1
    r.final <- r.future.real    
    raster::values(r.final)[intersect(which(!is.na(raster::values(r.present))), which(!is.na(raster::values(r.future.real))))] <- 0
    raster::values(r.final)[setdiff(which(!is.na(raster::values(r.present))), which(!is.na(raster::values(r.future.real))))] <- -1
#    plot(r.final)
    #####################get patches of future habitat
    patch.future <- clump(r.future.real, directions=8, gaps=F)
    # patch.future <- patchify(r.future.real, patch.dist)
    patch.freq.f <- freq(patch.future)
    # patch.freq.f <- patch.freq.f %>% data.frame %>% dplyr::filter(!is.na(value)) %>% arrange(desc(count)) %>% summarize(value, count, percent=count/sum(count))
    patch.freq.f <- data.frame(patch.freq.f)
    patch.freq.f <- patch.freq.f[which(!is.na(patch.freq.f$value)), ]
    patch.freq.f <- patch.freq.f[order(patch.freq.f$count, decreasing = TRUE), ]
    patch.freq.f$percent <- patch.freq.f$count / sum(patch.freq.f$count)
    
    ### fill out the data frame of species shift direction
    total.percent <- 0.0
    for (i in 1:3) {
      total.percent = total.percent + patch.freq.f$percent[i]
      # how many source patches we need? 
      idc <- which(raster::values(patch.future)==patch.freq.f$value[i])
      #
      ndescendant <- c()
      for (idp in 1:npatch.present) {
        # how many future cells are from the patch idp
        descendant <- lapply(which(patch.id==idp), function(x) idc.future[which(analogue.sd[which(idc.future %in% idc), x] <= dist.thrd)])
        ndescendant[idp] <- length(unlist(descendant))    # I do not need to use unique here.
      }
      ndescendant <- ndescendant / sum(ndescendant)
      # find maximum two source patches for each future patch
      iresult    <- iresult + 1
      xy.future  <- colMeans(xyFromCell(patch.future, idc))
      xy.present <- colMeans(xyFromCell(patch.present, which(raster::values(patch.present)==which.max(ndescendant))))
      shift.direction[iresult, 1]   <- scenarios[iscenario]
      shift.direction[iresult, 2:7] <- c(xy.future[1], xy.future[2], patch.freq.f$percent[i], xy.present[1], xy.present[2], max(ndescendant))
      # use an arrow to show the migration direction      
      # arrows(shift.direction[iresult,5], shift.direction[iresult,6], shift.direction[iresult,2], shift.direction[iresult,3])
      if (max(ndescendant) < 0.5) {
        # have second source patches
        iresult <- iresult + 1
        temp    <- sort(ndescendant, decreasing = T, index.return = T)
        xy.present <- colMeans(xyFromCell(patch.present, which(raster::values(patch.present)==temp$ix[2])))
        shift.direction[iresult, 1]   <- scenarios[iscenario]
        shift.direction[iresult, 2:7] <- c(xy.future[1], xy.future[2], patch.freq.f$percent[i], xy.present[1], xy.present[2], temp$x[2])
        # use an arrow to show the migration direction      
        # arrows(shift.direction[iresult,5], shift.direction[iresult,6], shift.direction[iresult,2], shift.direction[iresult,3])
      }
      if (total.percent > 0.5) {break}
    }
    #
    # Alternative method to calculate migration direction for each cell
    xy.new  <- xyFromCell(r.final, which(raster::values(r.final)==1))
    if (nrow(xy.new) > 0) {                                                  # add limitation
      bc.pres <- boundaries(r.present, type="inner")                         # edge cells with value = 1
      xy.bc.pres <- xyFromCell(r.present, which(raster::values(bc.pres)==1))
      # I have to set a limitation to migration direction, they only move polarward!
      id.source  <- pbsapply(1:nrow(xy.new), function(x) {
        # north hemisphere: choose bc points with lower y values
        if (xy.new[x,2] >=0) {
          id.selected <- which(xy.bc.pres[,2] < xy.new[x,2])
        } else {
          # south hemisphere: choose bc points with higher y values
          id.selected <- which(xy.bc.pres[,2] > xy.new[x,2])
        }
		# if no suitable points found, use all the boundary points
        if (length(id.selected) == 0) {
          id.selected <- 1:nrow(xy.bc.pres)
        }
        xy.bc.selected <- xy.bc.pres[id.selected,]
        return(id.selected[which.min(pointDistance(xy.new[x,], xy.bc.selected, lonlat=F, allpairs=F))])
      })
      if (nrow(xy.new) == 1) {                                              # deal with the case with only one cell 
        xy.source <- t(as.matrix(xy.bc.pres[id.source,]))
      } else {
        xy.source <- xy.bc.pres[id.source,]
      }
      # for testing purpose, I just plot many figures to look at this. 
      # for (ii in 1:(as.integer(nrow(xy.new)/1000) + 1)) { #as.integer(nrow(xy.new)/500)
      #   plot(r.final)
      #   istart = (ii-1)*1000+1
      #   iend   = min(ii*1000, nrow(xy.new))
      #   selected = seq(istart, iend, 1)
      #   arrows(xy.source[selected,1], xy.source[selected,2], xy.new[selected,1], xy.new[selected,2], length=0.1)
      # }
      # end of testing purpose
      distance  <- pointDistance(xy.new, xy.source, lonlat=F, allpairs=F) / 1000                     # convert this to km
      angle     <- atan2(xy.new[,2] - xy.source[,2], xy.new[,1] - xy.source[,1]) / 3.14159 * 180     # unit: degree. 
      # what if I get an average angle within a patch?
      r.new <- r.final
      raster::values(r.new)[which(raster::values(r.new)!=1)] <- NA
      # plot(r.new)
      patch.new <- clump(r.new, directions=8, gaps=F)
      # plot(patch.new)
      idp <- raster::extract(patch.new, xy.new)
      GID <- cellFromXY(r.template, xy.new)
      df  <- as.data.frame(cbind(idp, GID, distance, angle))
      df.tmp <- df %>% group_by(idp) %>% summarize(anglePavg = mean(angle), angleSD = sd(angle), num=n())
      df     <- left_join(df, df.tmp[,1:2], by='idp')
      df$idp <- iscenario
      colnames(df)[1] <- "iscenario"
      #
      shift.direction.cell <- rbind(shift.direction.cell, df)
    }
    #
    #
    # calculate statistics about species range shift
    # including current area, future area, future area with dispersal limitation; why not using the value == 1?
    area.shrink[1,1+iscenario]  <- sum(raster::values(area(r.future))[!is.na(raster::values(r.future))]) / area.shrink[1,1]
    area.shrink[1,17+iscenario] <- sum(raster::values(area(r.future.real))[!is.na(raster::values(r.future.real))]) / area.shrink[1,1]
    area.shrink[1,33+iscenario] <- sum(raster::values(area(r.final))[which(raster::values(r.final)==0)]) / area.shrink[1,1]    
  }
  
  # write csv files
  write.csv(shift.direction,      paste0(SDM.path, 'Realized/shift_direction/', name, '.csv'), row.names = F)
  write.csv(area.shrink,          paste0(SDM.path, 'Realized/area_shrink/', name, '.csv'), row.names = F)
  write.csv(shift.direction.cell, paste0(SDM.path, 'Realized/shift_direction_cell/', name, '.csv'), row.names = F)
}

#####################################################################################################################################
###########################################################start to call the function################################################
####get modelled species name; 
species <- list.files(paste0(SDM.path, "Results/RangePres/"), pattern = ".tif$")
species <- gsub('.{4}$', '', species)
# get species name from bioshift rate result
# nsp     <- nrow(shift.rates.plant) / 16
# species <- shift.rates.plant$Species[(1:nsp)*16]

# read species ID from command
command_args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(command_args[1])
# id <- id + 40000

# one command line calculate 8 species: ((id-1)*8+1):min((id*8),nsp)   
pblapply(species[id], get.realized.map)




