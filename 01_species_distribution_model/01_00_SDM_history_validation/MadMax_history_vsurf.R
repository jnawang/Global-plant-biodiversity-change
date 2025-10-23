
# Check = Condition for running SDMs with files exist
# proj.vars is a list
MadMax_history <- function(occ, occ.his, method, spname, fit.vars, proj.vars, proj.names, output_dir = "Results", check = F){
  
  dirs <- c("Maxnet","th","stats_hist",
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
    bb.buf <- rgeos::gBuffer(occs.sp, byid = TRUE, 
                             width = 2200*1000,)
    fut.buf <- rgeos::gBuffer(occs.sp, byid = TRUE, 
                             width = 2200*1000,)
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
    
    # browser()
    
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
    ##################
    # Run MaxEnt using SDMtune
    # predata
    # Create SWD object
    sdmdata <- prepareSWD(species = spname, 
                          p = occ, 
                          a = coordinates(bg.points.pres),
                          env = bg.ext.pres)

    # Split presence locations in training (80%) and testing (20%) datasets
    sdmdata_ <- trainValTest(sdmdata, test = 0.2, only_presence = TRUE)
    train <- sdmdata_[[1]]
    test <- sdmdata_[[2]]
    # browser()
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
    ## make sure both train and test sets have observed occurrence data! if not, use random method.
    cv <- 'Spatial_CV'
    if (sum(spatial_folds$records==0)>0) {
      spatial_folds <- randomFolds(train, k = 4, only_presence = TRUE, seed = 125)
      cv <- 'Random_CV'
      # random forests do not allow 0 record; I have to deal with this.
    }
    print(cv)
    ## done, 9/7/2024
    ### I will try four different methods
#    browser()
    if (method=='Maxent') {
	  cv_model <- SDMtune::train(method = method, data = train, folds = spatial_folds)
	  print('1')
#      cv_model <- SDMtune::train(method = method, data = train, folds = spatial_folds)
#	  sdmtune_model <- SDMtune::train(method = method, data = train)
      # 
#     # Search best tuning parameters: fc and reg
      h <- list(reg = seq(0.2, 5, 0.2), fc = c("l", "lq", "lh", "lqp", "lqph"))
#	  print('to random search')
      cv_model1 <- SDMtune::randomSearch(cv_model, hypers = h, metric = "auc", pop = 15, interactive = FALSE, 
                             seed = 65466)
      print('done random search')
      # final model
      sdmtune_model <- train(method = method, data = cv_model1@models[[1]]@data,
                      fc = cv_model1@results[1, 1], reg = cv_model1@results[1, 2])      
    } else if (method=='ANN') {
      cv_model <- train(method=method, data = train, size = 10, rang = 0.1, decay = 5e-4, maxit = 200, folds = spatial_folds)
      sdmtune_model <- train(method=method, data = train, size = 10, rang = 0.1, decay = 5e-4, maxit = 200)    # They do not need cross-validation
    } else if (method=="RF") {
      ## Train a Random Forest model
      cv_model <- train(method=method, data = train, ntree = 300, folds = spatial_folds)
      sdmtune_model <- train(method=method, data = train, ntree = 300)
    } else if (method=="RF_dsamp") {
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
      ##
      training   <- cbind(train@data[,ivarselec], train@pa)
      nenvar     <- ncol(training) - 1    
      ## revision ends here, but at another location
      colnames(training)[nenvar+1] <- 'occ'
      training$occ <- as.factor(training$occ)
      prNum <- as.numeric(table(training$occ)["1"]) # number of presences
      # the sample size in each class; the same as presence number
	  # make sure random size is not too big, < one third of training data;
      prNum   <- min(prNum, as.integer(nrow(training)/3))
      smpsize <- c("0" = prNum, "1" = prNum)
      rf_downsample <- randomForest(formula = occ ~., 
                                    data = training, 
                                    ntree = 500, 
                                    sampsize = smpsize,
                                    replace = TRUE)
    } else if (method=="BRT") {
      ## Train a Boosted Regression Tree model
      cv_model <- train(method=method, data = train, n.trees = 300, shrinkage = 0.001, folds = spatial_folds)
      sdmtune_model <- train(method=method, data = train, n.trees = 300, shrinkage = 0.001)
    }
	#
    # Calculate metrics to assess performance!
    if (method != "RF_dsamp"){
      AUC <- auc(sdmtune_model)
      AUC_test <- SDMtune::auc(sdmtune_model, test = test)
      AUC_CV_train <- SDMtune::auc(cv_model)
      AUC_CV_test  <- SDMtune::auc(cv_model, test=T)
      TSS <- SDMtune::tss(sdmtune_model)
    } else {
      # I need to calculate these metrics on my own
      # for the training data
      prob <- predict(rf_downsample, train@data, type = "prob")[,2]
      prob_occ <- prob[train@pa==1]
      prob_all <- prob[train@pa==0]
      AUC <- enmSdmX::evalAUC(prob_occ, prob_all)
      TSS <- enmSdmX::evalTSS(prob_occ, prob_all)
      # for the test data
      prob <- predict(rf_downsample, test@data, type = "prob")[,2]
      prob_occ <- prob[test@pa==1]
      prob_all <- prob[test@pa==0]
      AUC_test <- enmSdmX::evalAUC(prob_occ, prob_all)
      # for CV validation
      # information in spatialfold
      prob_occ_train <- c()
      prob_all_train <- c()
      prob_occ_test <- c()
      prob_all_test <- c()      
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
        prNum    <- as.numeric(table(training$occ[id.train])["1"]) # number of presences
        prNum    <- min(prNum, as.integer(nrow(training)/3))
		    smpsize  <- c("0" = prNum, "1" = prNum)
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
    }
    # browser()

    # Get threshold
    if (method=="RF" | method=="RF_dsamp") {
      th <- 0.5
    } else {
      thr <- SDMtune::thresholds(sdmtune_model, type = "cloglog", test = test)
      # Using Maximum training sensitivity plus specificity
      th <- thr$`Cloglog value`[3]      
    }
    print(th)
    
    # Predict the current distribution
    # Create SWD object
    if (method!="RF_dsamp") {
      sdm_pred_pres <- SDMtune::predict(sdmtune_model, 
                                       data = bg.ext.pres, 
                                       type = "cloglog")
    } else {
      # initialize
      sdm_pred_pres <- bg.ext.pres[[1]]
      present_env   <- as.data.frame(terra::values(bg.ext.pres))      
      # predict
      terra::values(sdm_pred_pres) <- predict(rf_downsample, present_env, type = "prob")[,2]
    }
#    plot(sdm_pred_pres)
#    plot(occs.sp, add=T)
    #
    # calculate Continuous Boyce Index (CBI), Junna 9/7/2024
    prob_occ <- na.omit(terra::extract(sdm_pred_pres, occ)[,2])
    prob_all <- na.omit(terra::extract(sdm_pred_pres, coordinates(bg.points.pres))[,1])
    CBI <- enmSdmX::evalContBoyce(prob_occ, prob_all)
    print(CBI) 
    # done! Junna 9/7/2024
    
#    browser()
    # do not use parallel
    sdm_pred_fut <- list()
    for (i in 1:length(proj.names)) {
      bg.ext.fut <- terra::crop(proj.vars[[i]], fut.buf)
      bg.ext.fut <- terra::mask(bg.ext.fut, terra::vect(fut.buf))
      # print the environmental results      
      if (method!="RF_dsamp") {
        sdm_pred_fut[[i]] <- SDMtune::predict(sdmtune_model,
                                              #                         data = proj.vars[[i]],
                                              data = bg.ext.fut,
                                              type = "cloglog")
      } else {
        sdm_pred_fut[[i]] <- bg.ext.fut[[1]]        # initialize it
        future_env   <- as.data.frame(terra::values(bg.ext.fut))
        terra::values(sdm_pred_fut[[i]]) <- predict(rf_downsample, future_env[,ivarselec], type = "prob")[,2]
      }
      # evaluate historical model here, using historical occ. 
      prob_occ <- na.omit(terra::extract(sdm_pred_fut[[i]], occ.his)[,2])
      prob_all <- na.omit(terra::extract(sdm_pred_fut[[i]], coordinates(bg.points.pres))[,1])
      # different evaluation matrics
      AUC.his <- enmSdmX::evalAUC(prob_occ, prob_all)
      TSS.his <- enmSdmX::evalTSS(prob_occ, prob_all)
      CBI.his <- enmSdmX::evalContBoyce(prob_occ, prob_all)
      print(paste(AUC.his, TSS.his, CBI.his, sep='_'))
    }
#    plot(sdm_pred_fut[[1]])


    # browser()
    ############################
    # save outputs
    ## model stats
    ## current area, area of 6 scenarios, change percent of 6 scenarios
    write.table(
      data.frame(sps = spname,
                 nOC = nrow(occ), 
                 nc.cur = sum(terra::values(sdm_pred_pres)>=th, na.rm=T), 
                 nc.his = sum(terra::values(sdm_pred_fut[[1]])>=th, na.rm=T), 
                 AUC      = AUC,
                 AUC_test = AUC_test,
                 AUC_CV_train = AUC_CV_train,
                 AUC_CV_test  = AUC_CV_test,
				         TSS      = TSS,
                 CBI      = CBI,
                 AUC.his  = AUC.his,
                 TSS.his  = TSS.his,  
                 CBI_his  = CBI.his),
      file = paste0(output_dir,"/stats_hist/",method,"_",spname,"_stats.csv"),
      sep = ",", row.names = F)
    
    # Maps future
    lapply(1:length(proj.names), function(x){
      outdir <- paste0(output_dir,"/",proj.names[x])
      #
      # range future
      raster::writeRaster(
        sdm_pred_fut[[x]] >= th,
        paste0(output_dir,"/",proj.names[x],"/RangeFut/",method,"_",spname,".tif"), overwrite=TRUE)
    }
    )
  }
}

