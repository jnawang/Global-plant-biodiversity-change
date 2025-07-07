
# Check = Condition for running SDMs with files exist
# proj.vars is a list
MadMax <- function(occ, spname, fit.vars, proj.vars, proj.names, output_dir = "Results", check = F){
  
  dirs <- c("maxe_rfds","th","stats",
            "LogOutPres","RangePres", "RangePresNoomit", "RangePresEss", "RangePres95s", "Occ")
  dirs.fut <- c("LogOutFut","RangeFut", "RangeFutNoomit", "RangeFutEss", "RangeFut95s")
  
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
    } else {
      lapply(dirs.fut, function(x){
        dir.test.th <- paste0(dir.test,"/",x)
        if (!dir.exists(dir.test.th)) {
          dir.create(dir.test.th, recursive = TRUE)
        }
      })
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
    occs.sp <- sp::SpatialPoints(occ, proj4string = CRS(terra::crs(fit.vars, proj=T)))
    # create two buffers around occurrences; one for background, the other for future range
    # buffer uses unit in meter. 5 degrees (~550km) = 550000m; future range 20 degree (~2200 km)
    bb.buf  <- terra::buffer(occs.sp, width = 4*550*1000)
    fut.buf <- bb.buf
    # aggregate buffers
    bb.buf  <- terra::aggregate(bb.buf, dissolve = TRUE)
    fut.buf <- terra::aggregate(fut.buf, dissolve = TRUE)
    # crop and mask environmental rasters of the present by the extent to make the background extent
    bg.ext.pres <- terra::crop(fit.vars, bb.buf)
    bg.ext.pres <- terra::mask(bg.ext.pres, terra::vect(bb.buf))
    # 
    # generate up to 20,000 random points within the background extent of the present
    # https://github.com/jamiemkass/ENMeval/issues/26       How to create the bias probability file?
    # dismo: has the function randomPoints, so sampleRandom function is no longer used. 
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
    # browser()
    #
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
#    train <- sdmdata                if we assume th = 0.5; we use this line, but using th=0.5 affects richness calculation.         
#    browser()
    # Train Maxent model with cross-validation
    #    folds <- randomFolds(train, k = 2, only_presence = TRUE, seed = 125)
    # spatially blocked cross validation, added on 9/7/2024
    sf_df <- sf::st_as_sf(cbind(train@coords, pa = train@pa),
                          coords = c("X", "Y"),
                          crs = terra::crs(fit.vars,
                                           proj = TRUE))
    spatial_folds <- cv_spatial(x = sf_df,
                                column = "pa",
                                rows_cols = c(10, 10),      # I can set size based on correlation distance.
                                k = 4,
                                hexagon = FALSE,             # this is the shape of the block
                                selection = "random",
                                plot = TRUE)               # "systematic"; 'random'
    cv <- 'Spatial_CV'
    ## make sure both train and test sets have observed occurrence data! if not, use random method.
    if (sum(spatial_folds$records==0)>0) {
      spatial_folds <- randomFolds(train, k = 4, only_presence = TRUE, seed = 125)
      cv <- 'Random_CV'
      # random forests do not allow 0 record; I have to deal with this.
    }
    print(cv)
#    browser()
    ####################################################################################################################
    #######################################################RUN Maxent first#############################################
    cv_model <- SDMtune::train(method = "Maxent", data = train, folds = spatial_folds)
    # Search best tuning parameters: fc and reg
    h <- list(reg = seq(0.2, 5, 0.2), fc = c("l", "lq", "lh", "lqp", "lqph"))       # , "lqpth"
    cv_model1 <- randomSearch(cv_model, hypers = h, metric = "auc", pop = 20, interactive = F,
                              seed = 65466)
    #
    # final maxent model
    sdmtune_model <- train(method = "Maxent", data = cv_model1@models[[1]]@data,
                           fc = cv_model1@results[1, 1], reg = cv_model1@results[1, 2])
    #
    # Predict the current distribution
    sdm_pred_pres_maxe = SDMtune::predict(sdmtune_model, 
                                     data = bg.ext.pres, 
                                     type = "cloglog")
    #
    prob_occ <- na.omit(terra::extract(sdm_pred_pres_maxe, occ)[,2])
    prob_all <- na.omit(terra::extract(sdm_pred_pres_maxe, coordinates(bg.points.pres))[,1])
    AUC_maxe <- enmSdmX::evalAUC(prob_occ, prob_all)
    TSS_maxe <- enmSdmX::evalTSS(prob_occ, prob_all)
    CBI_maxe <- enmSdmX::evalContBoyce(prob_occ, prob_all)
    print(CBI_maxe)
    ###########
    #######################################################RUN Random Forest second#############################################
    if (file.exists(paste0("../SDM2024RFDS/",output_dir,"/stats/",spname,"_stats.csv"))) {
      stats     <- read.csv(paste0("../SDM2024RFDS/",output_dir,"/stats/",spname,"_stats.csv"))
      ivarselec <- which(names(fit.vars) %in% stats$sps[3:nrow(stats)])
      print(ivarselec)
    } else {
      ivarselec <- 1:15      # use all variables if this species is not modelled before. 
    }
    ##
    training   <- cbind(train@data[,ivarselec], train@pa)
    nenvar     <- ncol(training) - 1
    ##
    colnames(training)[nenvar+1] <- 'occ'
    training$occ <- as.factor(training$occ)
    prNum <- as.numeric(table(training$occ)["1"])    # number of presences
    prNum <- min(prNum, as.integer(nrow(training)/3))
    # the sample size in each class; the same as presence number
    smpsize <- c("0" = prNum, "1" = prNum)
    rf_downsample <- randomForest(formula = occ ~.,
                                  data = training,
                                  ntree = 500,
                                  sampsize = smpsize,
                                  replace = TRUE)
    ##
    #######################################################################################################
    # Predict the current distribution by rf_downsample
    sdm_pred_pres_rfds <- bg.ext.pres[[1]]
    present_env   <- as.data.frame(terra::values(bg.ext.pres))
    # predict
    terra::values(sdm_pred_pres_rfds) <- predict(rf_downsample, present_env[,ivarselec], type = "prob")[,2]
    #
    # calculate CBI
    prob_occ <- na.omit(terra::extract(sdm_pred_pres_rfds, occ)[,2])
    prob_all <- na.omit(terra::extract(sdm_pred_pres_rfds, coordinates(bg.points.pres))[,1])
    AUC_rfds <- enmSdmX::evalAUC(prob_occ, prob_all)
    TSS_rfds <- enmSdmX::evalTSS(prob_occ, prob_all)
    CBI_rfds <- enmSdmX::evalContBoyce(prob_occ, prob_all)
    print(CBI_rfds)    
    ####
    ##########################################ensemble average of current prediction########################
    if (CBI_rfds > 0 & CBI_maxe > 0) {
      weight_rfds <- CBI_rfds
      weight_maxe <- CBI_maxe
    } else if (CBI_rfds > 0) {
      weight_rfds <- CBI_rfds
      weight_maxe <- 0      
    } else if (CBI_maxe > 0) {
      weight_rfds <- 0
      weight_maxe <- CBI_maxe
    } else {
      weight_rfds <- CBI_rfds
      weight_maxe <- 0
    }
    
    # weighted average
    sdm_pred_pres <- (sdm_pred_pres_rfds * weight_rfds + sdm_pred_pres_maxe * weight_maxe) / (weight_rfds + weight_maxe)
    
    # performance of ensemble average
    prob_occ <- na.omit(terra::extract(sdm_pred_pres, occ)[,2])
    prob_all <- na.omit(terra::extract(sdm_pred_pres, coordinates(bg.points.pres))[,1])
    AUC_maxe_rfds <- enmSdmX::evalAUC(prob_occ, prob_all)
    TSS_maxe_rfds <- enmSdmX::evalTSS(prob_occ, prob_all)
    CBI_maxe_rfds <- enmSdmX::evalContBoyce(prob_occ, prob_all)
    
    # if we get bad ensemble average results, use the best model of the two
    if (CBI_maxe_rfds < 0.7 & CBI_maxe_rfds < max(CBI_rfds, CBI_maxe)) {
      # use the best models
      if (CBI_rfds > CBI_maxe) {
        weight_rfds <- CBI_rfds
        weight_maxe <- 0
        AUC_maxe_rfds <- AUC_rfds
        TSS_maxe_rfds <- TSS_rfds
        CBI_maxe_rfds <- CBI_rfds
      } else {
        weight_rfds <- 0
        weight_maxe <- CBI_maxe
        AUC_maxe_rfds <- AUC_maxe
        TSS_maxe_rfds <- TSS_maxe
        CBI_maxe_rfds <- CBI_maxe
      }
      
      # recalculate current distribution
      sdm_pred_pres <- (sdm_pred_pres_rfds * weight_rfds + sdm_pred_pres_maxe * weight_maxe) / (weight_rfds + weight_maxe)
    }
    
    ####
    #################################################get threshold values###################################
    # browser()
    data <- sdmdata
    n_p <- sum(data@pa == 1)
    n_a <- sum(data@pa == 0)
    #
    pred_rfds <- predict(rf_downsample, data@data[,ivarselec], type = "prob")[,2]
    pred_maxe <- SDMtune::predict(sdmtune_model, data = data, type = "cloglog") 
    
    # get ensemble average    
    pred <- (pred_rfds * weight_rfds + pred_maxe * weight_maxe) / (weight_rfds + weight_maxe)
    
    # calculate threshold
    p_pred <- pred[1:n_p]
    a_pred <- pred[(n_p+1):(n_p+n_a)]
    #
    e <- evaluate(p=p_pred, a=a_pred)
    thr <- threshold(e, sensitivity=0.95)
    thr <- pmax(thr, 0.05)        # the lower limit of occurrence probability is 0.05. 
    th_mss <- as.numeric(thr[2])
    th_noomit <- as.numeric(thr[3])
    th_ess <- as.numeric(thr[5])
    th_95s <- as.numeric(thr[6])
    # 
    ##############
    ###############################################predict for the future#####################################
    ####only work on ESM cases to save space, 3/9/2025: 
    # ESM.id <- which(grepl('^ESM_', proj.names))     # based on the order the variables are stored. 
    # proj.names <- proj.names[ESM.id]
    # print(proj.names)
    ####
    sdm_pred_fut <- list()
    for (i in 1:length(proj.names)) {
      bg.ext.fut <- terra::crop(proj.vars[[i]], fut.buf)
      bg.ext.fut <- terra::mask(bg.ext.fut, terra::vect(fut.buf))
      #### predict by maxent
      sdm_pred_fut[[i]] <- SDMtune::predict(sdmtune_model, data = bg.ext.fut, type = "cloglog")
      #
      #### predict by rfds                    
      future_env       <- as.data.frame(terra::values(bg.ext.fut))
      future_prob_rfds <- predict(rf_downsample, future_env[,ivarselec], type = "prob")[,2]
      #
      # combine the two results
      terra::values(sdm_pred_fut[[i]]) <- (terra::values(sdm_pred_fut[[i]]) * weight_maxe + future_prob_rfds * weight_rfds) / (weight_rfds + weight_maxe)
    }
    ####################################################################################
    #
    ############################
    # save outputs
    ## model, do not save these models because they are too large.
    # saveRDS(
    #   rf_downsample,
    #   paste0(output_dir,"/maxe_rfds/",spname,"_mod_rfds.RDS"))
    # saveRDS(
    #   sdmtune_model,
    #   paste0(output_dir,"/maxe_rfds/",spname,"_mod_maxe.RDS"))    
    
    ## thresholds
    write.csv(
      data.frame(sps = spname,
                 th = thr),
      paste0(output_dir,"/th/",spname,"_thr.csv"),
      row.names = F)
    
    ## model stats
    write.table(
      data.frame(sps = spname,
                 nOC = nrow(occ),
                 AUC_maxe = AUC_maxe,
                 TSS_maxe = TSS_maxe,
                 CBI_maxe = CBI_maxe,                 
                 AUC_rfds = AUC_rfds,
                 TSS_rfds = TSS_rfds,
                 CBI_rfds = CBI_rfds,                 
                 AUC = AUC_maxe_rfds,
                 TSS = TSS_maxe_rfds,
                 CBI = CBI_maxe_rfds),
      file = paste0(output_dir,"/stats/",spname,"_stats.csv"),
      sep = ",", row.names = F)
    
    # logistic output
    terra::writeRaster(
      sdm_pred_pres,
      paste0(output_dir, "/LogOutPres/",spname,".tif"), overwrite=TRUE)
    # range
    terra::writeRaster(
      sdm_pred_pres >= th_mss,
      paste0(output_dir,"/RangePres/",spname,".tif"), overwrite=TRUE)
    # more present ranges
    # terra::writeRaster(
    #   sdm_pred_pres >= th_noomit,
    #   paste0(output_dir,"/RangePresNoomit/",spname,".tif"), overwrite=TRUE)
    #
    terra::writeRaster(
      sdm_pred_pres >= th_ess,
      paste0(output_dir,"/RangePresEss/",spname,".tif"), overwrite=TRUE)
    #
    terra::writeRaster(
      sdm_pred_pres >= th_95s,
      paste0(output_dir,"/RangePres95s/",spname,".tif"), overwrite=TRUE)
    #
    # Maps future
    lapply(1:length(proj.names), function(x){

      # logistic output future
      if (grepl('^ESM_', proj.names[x])) {
        terra::writeRaster(
          sdm_pred_fut[[x]],
          paste0(output_dir,"/",proj.names[x],"/LogOutFut/",spname,".tif"), overwrite=TRUE)
      }
      # range future
      terra::writeRaster(
        sdm_pred_fut[[x]] >= th_mss,
        paste0(output_dir,"/",proj.names[x],"/RangeFut/",spname,".tif"), overwrite=TRUE)
      # more range future
      # terra::writeRaster(
      #   sdm_pred_fut[[x]] >= th_noomit,
      #   paste0(output_dir,"/",proj.names[x],"/RangeFutNoomit/",spname,".tif"), overwrite=TRUE)
      # 
      terra::writeRaster(
        sdm_pred_fut[[x]] >= th_ess,
        paste0(output_dir,"/",proj.names[x],"/RangeFutEss/",spname,".tif"), overwrite=TRUE)
      # 
      terra::writeRaster(
        sdm_pred_fut[[x]] >= th_95s,
        paste0(output_dir,"/",proj.names[x],"/RangeFut95s/",spname,".tif"), overwrite=TRUE)
      # 
    }
    )
  }
}

