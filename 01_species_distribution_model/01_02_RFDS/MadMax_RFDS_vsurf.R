
# Check = Condition for running SDMs with files exist
# proj.vars is a list
MadMax <- function(occ, spname, fit.vars, proj.vars, proj.names, output_dir = "Results", check = F){
  
  dirs <- c("RFDS","th","stats",
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
    
    # RF method does not need this.
    # # Split presence locations in training (80%) and testing (20%) datasets
    # sdmdata_ %<-% trainValTest(sdmdata, test = 0.2, only_presence = TRUE)
    # train <- sdmdata_[[1]]
    # test <- sdmdata_[[2]]
    train <- sdmdata
    
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
    #    browser()
    if (sum(spatial_folds$records==0)>0) {
      spatial_folds <- randomFolds(train, k = 4, only_presence = TRUE, seed = 125)
      cv <- 'Random_CV'
      # random forests do not allow 0 record; I have to deal with this.
    }
    print(cv)
    ####################################################
    ## Revised to include RSURF, 11/16/2024
    varselec <- VSURF(x=train@data, y=train@pa, ncores = detectCores(), parallel = TRUE)
    if (varselec$nums.varselect[3]==0) {
      if (varselec$nums.varselect[2]==0) {
        if (varselec$nums.varselect[1]==0) {
          ivarselec <- 1:ncol(train@data)
        } else {
          ivarselec <- varselec$varselect.thres
        }
      } else {
        ivarselec <- varselec$varselect.interp
      }
    } else {
      ivarselec <- varselec$varselect.pred
    }
    print(varselec$nums.varselect)
    #      browser()
    ##
    training   <- cbind(train@data[,ivarselec], train@pa)
    nenvar     <- ncol(training) - 1
    ##################################end of revision + L231###############    
    colnames(training)[nenvar+1] <- 'occ'
    training$occ <- as.factor(training$occ)
    prNum <- as.numeric(table(training$occ)["1"])# number of presences
    prNum <- min(prNum, as.integer(nrow(training)/3))
    # the sample size in each class; the same as presence number
    smpsize <- c("0" = prNum, "1" = prNum)
    rf_downsample <- randomForest(formula = occ ~., 
                                  data = training, 
                                  ntree = 500, 
                                  sampsize = smpsize,
                                  replace = TRUE)
    # check model performance, AUC of cross validation
    prob <- predict(rf_downsample, train@data, type = "prob")[,2]
    prob_occ <- prob[train@pa==1]
    prob_all <- prob[train@pa==0]
    AUC <- enmSdmX::evalAUC(prob_occ, prob_all)
    TSS <- enmSdmX::evalTSS(prob_occ, prob_all)
    # # for the test data
    # prob <- predict(rf_downsample, test@data, type = "prob")[,2]
    # prob_occ <- prob[test@pa==1]
    # prob_all <- prob[test@pa==0]
    # AUC_test <- enmSdmX::evalAUC(prob_occ, prob_all)
    # for CV validation
    # information in spatialfold
    prob_occ_train <- c()
    prob_all_train <- c()
    prob_occ_test  <- c()
    prob_all_test  <- c()      
    for (k in 1:4) {
      # deal with different CV methods
      if (cv=='Spatial_CV') {
        id.train <- spatial_folds$folds_list[[k]][[1]]
        id.test  <- spatial_folds$folds_list[[k]][[2]]          
      } else if (cv=='Random_CV') {
        id.train <- spatial_folds$train[,k]
        id.test  <- spatial_folds$test[,k]
      }
      # train a model
      # the sample size in each class; the same as presence number
      prNum <- as.numeric(table(training[id.train,]$occ)["1"]) # number of presences
      prNum <- min(prNum, as.integer(nrow(training)/3))
      smpsize <- c("0" = prNum, "1" = prNum)        
      cv_model <- randomForest(formula = occ ~., 
                               data = training[id.train,], 
                               ntree = 500, 
                               sampsize = smpsize,
                               replace = TRUE)
      pred.train <- predict(cv_model, training[id.train, 1:nenvar], type = "prob")[,2]
      prob_occ_train <- c(prob_occ_train, pred.train[training$occ[id.train]=='1'])
      prob_all_train <- c(prob_all_train, pred.train[training$occ[id.train]=='0'])
      # test a model
      pred.test  <- predict(cv_model, training[id.test, 1:nenvar], type = "prob")[,2]
      prob_occ_test <- c(prob_occ_test, pred.test[training$occ[id.test]=='1'])
      prob_all_test <- c(prob_all_test, pred.test[training$occ[id.test]=='0'])
    }
    AUC_CV_train <- enmSdmX::evalAUC(prob_occ_train, prob_all_train)
    AUC_CV_test  <- enmSdmX::evalAUC(prob_occ_test, prob_all_test)
    #
    # get relative importance of env variables
    vi <- importance(rf_downsample) 
    vi <- vi / sum(vi) * 100    # relative importance
    print(vi)
    #  
    ################calculate the threshold, using the same method from SDMtune, Junna 2/22/2025. 
    # browser()
    data <- sdmdata
    n_p <- sum(data@pa == 1)
    n_a <- sum(data@pa == 0)
    #
    pred <- predict(rf_downsample, data@data[,ivarselec], type = "prob")[,2]
    p_pred <- pred[1:n_p]
    a_pred <- pred[(n_p + 1):(n_p + n_a)]    
    #
    e <- evaluate(p=p_pred, a=a_pred)
    thr <- threshold(e, sensitivity=0.95)
    thr <- pmax(thr, 0.05)        # the lower limit of occurrence probability is 0.05. 
    th_mss <- as.numeric(thr[2])
    th_noomit <- as.numeric(thr[3])
    th_ess <- as.numeric(thr[5])
    th_95s <- as.numeric(thr[6])
    ######################################################the end of calculating thresholds#############################################
    #    
    # Predict the current distribution
    sdm_pred_pres <- bg.ext.pres[[1]]
    present_env   <- as.data.frame(terra::values(bg.ext.pres))      
    # predict
    terra::values(sdm_pred_pres) <- predict(rf_downsample, present_env[,ivarselec], type = "prob")[,2]
    #
    # calculate CBI
    prob_occ <- na.omit(terra::extract(sdm_pred_pres, occ)[,2])
    prob_all <- na.omit(terra::extract(sdm_pred_pres, coordinates(bg.points.pres))[,1])
    CBI <- enmSdmX::evalContBoyce(prob_occ, prob_all)
    print(CBI) 
    #
    sdm_pred_fut <- list()
    for (i in 1:length(proj.names)) {
      bg.ext.fut <- terra::crop(proj.vars[[i]], fut.buf)
      bg.ext.fut <- terra::mask(bg.ext.fut, terra::vect(fut.buf))
      # print the environmental results      
      #      writeRaster(bg.ext.fut, 'C:/Users/jnawang/Box/Biodiversity/SDM/Results/test/bg.ext.fut.tif', overwrite=T)
      sdm_pred_fut[[i]] <- bg.ext.fut[[1]]        # initialize it
      future_env        <- as.data.frame(terra::values(bg.ext.fut))
      terra::values(sdm_pred_fut[[i]]) <- predict(rf_downsample, future_env[,ivarselec], type = "prob")[,2]
    }
    ####################################################################################
    #
    ############################
    # save outputs
    ## model
    saveRDS(
      rf_downsample,
      paste0(output_dir,"/RFDS/",spname,"_mod.RDS"))
    # thresholds
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
                 AUC = AUC,
                 AUC_test = 0.0,
                 AUC_CV_train = AUC_CV_train,
                 AUC_CV_test  = AUC_CV_test,
                 TSS = TSS,
                 CBI = CBI),
      file = paste0(output_dir,"/stats/",spname,"_stats.csv"),
      sep = ",", row.names = F)
    # write area of species range, ratio of species range shift
#    range_fut <- rbind(sdarea_fut, sdshift_fut)
#    colnames(range_fut) <- proj.names
#    write.table(range_fut, file=paste0(output_dir,"/stats/",spname,"_stats.csv"),
#                sep = ",", row.names = F, append=TRUE)
    # write variable importance
    write.table(vi, file=paste0(output_dir,"/stats/",spname,"_stats.csv"),
                sep = ",", row.names = T, append=TRUE)


    # logistic output
   terra::writeRaster(
     sdm_pred_pres,
     paste0(output_dir, "/LogOutPres/",spname,".tif"), overwrite=TRUE)
    # range
    terra::writeRaster(
      sdm_pred_pres >= th_95s,
      paste0(output_dir,"/RangePres/",spname,".tif"), overwrite=TRUE)

    # occurrence
    png(file=paste0(output_dir,"/Occ/",spname,".png"), width=600, height=350)
    plot(sdm_pred_pres)
    plot(occs.sp, add=TRUE, pch = '.')
    dev.off()

    # Maps future
    lapply(1:length(proj.names), function(x){

      # logistic output future
     terra::writeRaster(
       sdm_pred_fut[[x]],
       paste0(output_dir,"/",proj.names[x],"/LogOutFut/",spname,".tif"), overwrite=TRUE)
      # range future
      terra::writeRaster(
        sdm_pred_fut[[x]] >= th_95s,
        paste0(output_dir,"/",proj.names[x],"/RangeFut/",spname,".tif"), overwrite=TRUE)
    }
    )
  }
}

