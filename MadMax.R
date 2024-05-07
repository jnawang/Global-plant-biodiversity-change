
# Check = Condition for running SDMs with files exist
# proj.vars is a list
MadMax <- function(occ, spname, fit.vars, proj.vars, proj.names, output_dir = "Results", check = F){
  
  dirs <- c("Maxnet","th","stats",
            "LogOutPres","RangePres","Occ")
  dirs.fut <- c("LogOutFut","RangeFut")
  
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  lapply(proj.names, function(i){
    dir.test <- paste(output_dir,i,sep="/")
    if(!dir.exists(dir.test)){
      dir.create(dir.test, recursive = TRUE)
      lapply(dirs.fut, function(x){
        dir.create(paste0(dir.test,"/",x), recursive = TRUE)
      }
      )
    }
  })
  lapply(dirs, function(i){
    dir.test <- paste(output_dir,i,sep="/")
    if(!dir.exists(dir.test)){
      dir.create(dir.test, recursive = TRUE)
    }
  })
  
  require(SDMtune)
  require(dismo)
  require(sp)
  require(raster)
  require(terra)
  
  # check?
  if(check){
    mycheck <- sapply(dirs, function(x){
      my.files <- list.files(paste0(output_dir,"/",x))
      grep(spname,my.files)
    })
    mycheck <- any(is.na(mycheck>0))
  }
  if(!check){
    mycheck <- TRUE
  }
  
  if(mycheck){
    
    # make a spatial points object for the occurrences
    occs.sp <- sp::SpatialPoints(occ,proj4string = CRS(terra::crs(fit.vars, proj=T)))
    # create two buffers around occurrences; one for background, the other for future range
    # buffer uses unit in meter. 5 degrees (~550km) = 550000m; future range 20 degree (~2200 km)
    bb.buf <- rgeos::gBuffer(occs.sp, byid = TRUE, 
                             width = 10*550*1000,)
    fut.buf <- rgeos::gBuffer(occs.sp, byid = TRUE, 
                             width = 10*550*1000,)
    # aggregate buffers
    bb.buf  <- terra::aggregate(bb.buf, dissolve = TRUE)
    fut.buf <- terra::aggregate(fut.buf, dissolve = TRUE)
    # crop and mask environmental rasters of the present by the extent to make the background extent
    bg.ext.pres <- terra::crop(fit.vars, bb.buf)
    bg.ext.pres <- terra::mask(bg.ext.pres, terra::vect(bb.buf))
    # 
    # # generate up to 20,000 random points within the background extent of the present
    # # https://github.com/jamiemkass/ENMeval/issues/26       How to create the bias probability file?
    # # dismo: has the function randomPoints
    # bg.points.pres <- raster::sampleRandom(x=bg.ext.pres,
    #                                        size=20000,
    #                                        na.rm=T, #removes the 'Not Applicable' points
    #                                        sp=T) # return spatial points
    sample.bias <- raster::raster(paste0(output_dir, "/sample_density.tif"))
    sample.bias <- raster::crop(sample.bias, bb.buf)
    sample.bias <- raster::mask(sample.bias, bb.buf)
    
    bg.points.pres <- try(dismo::randomPoints(sample.bias, 20000, occs.sp, excludep=TRUE, prob=TRUE,
                 cellnumbers=FALSE, warn=2), silent = T)
    if (class(bg.points.pres)[1] == "try-error") {
      # reduce number of sampling points
      bg.points.pres <- try(dismo::randomPoints(sample.bias, sum(!is.na(raster::values(sample.bias)))-length(occs.sp)-200, occs.sp, excludep=TRUE, prob=TRUE,
                                                cellnumbers=FALSE, warn=2), silent = T)
      if (class(bg.points.pres)[1] == "try-error") {
        # if it still has an error, we use random sampling without considering bias information. 
        # bg.points.pres <- raster::sampleRandom(x=bg.ext.pres, size=20000, na.rm=T, sp=T)
		bg.points.pres <- dismo::randomPoints(sample.bias, 20000, occs.sp, excludep=TRUE, prob=FALSE, cellnumbers=FALSE, warn=2)
      }
    }

    ##################
    # Run MaxEnt using SDMtune
    # predata
    # Create SWD object
    sdmdata <- prepareSWD(species = spname, 
                          p = occ, 
                          a = coordinates(bg.points.pres),
                          env = bg.ext.pres)

    # Split presence locations in training (80%) and testing (20%) datasets
    sdmdata_ %<-% trainValTest(sdmdata, test = 0.2, only_presence = TRUE)
    train <- sdmdata_[[1]]
    test <- sdmdata_[[2]]
    
    # Train Maxent model with cross-validation
    folds <- randomFolds(train, k = 4, only_presence = TRUE, seed = 125)
    cv_model <- SDMtune::train(method = "Maxent", data = train, folds = folds)
    
    # Search best tuning parameters: fc and reg
    h <- list(reg = seq(0.2, 5, 0.2), fc = c("l", "lq", "lh", "lqp", "lqph"))
    cv_model1 <- randomSearch(cv_model, hypers = h, metric = "auc", pop = 15, interactive = F,
                           seed = 65466)
    
    # final model
    sdmtune_model <- train(method = "Maxent", data = cv_model1@models[[1]]@data,
                    fc = cv_model1@results[1, 1], reg = cv_model1@results[1, 2])
    print('AUC of train data')
    print(auc(sdmtune_model))
    print('AUC of test data')
    print(auc(sdmtune_model, test = test))

    # importance of env variables
    vi <- maxentVarImp(sdmtune_model)
    print(vi)
#    plotResponse(sdmtune_model, var = vi[1,1], type = "cloglog",
#                 only_presence = TRUE, marginal = TRUE, fun = mean, color = "blue")

    # Get threshold
    thr <- SDMtune::thresholds(sdmtune_model, type = "cloglog", test = test)
    # Using Maximum training sensitivity plus specificity
    th <- thr$`Cloglog value`[3]

    # Predict the current distribution
    # Create SWD object
    data_full <- prepareSWD(species = spname, 
                            p = occ, 
                            a = coordinates(bg.points.pres),
                            env = bg.ext.pres)
    
    sdm_pred_pres = SDMtune::predict(sdmtune_model, 
                                     data = bg.ext.pres, 
                                     type = "cloglog")
    
    # Junna simplified the code, not using parallel
	sdm_pred_fut <- list()
	for (i in 1:length(proj.names)) {
        bg.ext.fut <- terra::crop(proj.vars[[i]], fut.buf)
        bg.ext.fut <- terra::mask(bg.ext.fut, terra::vect(fut.buf))
        sdm_pred_fut[[i]] <- SDMtune::predict(sdmtune_model,
#                         data = proj.vars[[i]],
                         data = bg.ext.fut,
                         type = "cloglog")
	}

	
    ## Calculate species distribution area
#    sdarea_pres <- sum(raster::values(area(sdm_pred_pres))[which(raster::values(sdm_pred_pres>th))])
#    sdarea_fut <- rep(0.0, length(sdm_pred_fut))
#    sdshift_fut <- rep(0.0, length(sdm_pred_fut))
#    for (i in 1:length(sdm_pred_fut)) {
#      sdarea_fut[i] <- sum(raster::values(area(sdm_pred_pres))[which(raster::values(sdm_pred_fut[[i]]>=th))])
#      sdshift_fut[i] <- 1.0 - sum(raster::values(area(sdm_pred_pres))[which(raster::values(sdm_pred_fut[[i]]>=th) & raster::values(sdm_pred_pres>=th))]) / sdarea_pres
#    }

    ## get the distribution matrix: this matrix includes all the occurence and distribution information
#    sdmatrix <- matrix(NA, nrow = 14, ncol = ncell(bios_pres[[1]]))
#    sdmatrix[1,] <- getValues(sdm_pred_pres)
#    sdmatrix[2:7,] <- getValues(sdm_pred_fut[[1:6]])
#    sdmatrix[8,] <- getValues(sdm_pred_pres) >= th
#    sdmatrix[9:14,] <- getValues(sdm_pred_pres[[1:6]]) >= th


    ############################
    # save outputs
    ## model
    saveRDS(
      sdmtune_model,
      paste0(output_dir,"/Maxnet/",spname,"_mod.RDS"))
    ## thresholds
    write.csv(
      data.frame(sps = spname,
                 th = thr),
      paste0(output_dir,"/th/",spname,"_thr.csv"),
      row.names = F)
    ## model stats
    ## current area, area of 6 scenarios, change percent of 6 scenarios
    write.table(
      data.frame(sps = spname,
                 nOC = nrow(occ),
                 AUC      = SDMtune::auc(sdmtune_model),
                 AUC_test = SDMtune::auc(sdmtune_model,test),
                 TSS      = SDMtune::tss(sdmtune_model),
                 AICc     = SDMtune::aicc(sdmtune_model, bg.ext.pres)),
      file = paste0(output_dir,"/stats/",spname,"_stats.csv"),
      sep = ",", row.names = F)
    # write area of species range, ratio of species range shift
#    range_fut <- rbind(sdarea_fut, sdshift_fut)
#    colnames(range_fut) <- proj.names
#    write.table(range_fut, file=paste0(output_dir,"/stats/",spname,"_stats.csv"),
#                sep = ",", row.names = F, append=TRUE)
    # write variable importance
    write.table(vi, file=paste0(output_dir,"/stats/",spname,"_stats.csv"),
                sep = ",", row.names = F, append=TRUE)


    # logistic output
#    terra::writeRaster(
#      sdm_pred_pres,
#      paste0(output_dir, "/LogOutPres/",spname,".tif"), overwrite=TRUE)
    # range
    terra::writeRaster(
      sdm_pred_pres >= th,
      paste0(output_dir,"/RangePres/",spname,".tif"), overwrite=TRUE)

    # occurrence
    png(file=paste0(output_dir,"/Occ/",spname,".png"), width=600, height=350)
    plot(sdm_pred_pres)
    plot(occs.sp, add=TRUE, pch = '.')
    dev.off()

    # Maps future
    lapply(1:length(proj.names), function(x){

      # logistic output future
#      terra::writeRaster(
#        sdm_pred_fut[[x]],
#        paste0(output_dir,"/",proj.names[x],"/LogOutFut/",spname,".tif"), overwrite=TRUE)    
      # range future
      terra::writeRaster(
        sdm_pred_fut[[x]] >= th,
        paste0(output_dir,"/",proj.names[x],"/RangeFut/",spname,".tif"), overwrite=TRUE)
    }
    )
  }
}

