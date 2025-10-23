rm(list=ls())
gc()

########################
# Setup
list.of.packages <- c("raster", "terra", "blockCV", "enmSdmX", "randomForest", "VSURF",
                      "SDMtune","dismo",
                      "tidyverse",
                      "foreach","doParallel","pbapply","doFuture","future","doRNG")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

sapply(list.of.packages, require, character.only = TRUE)

########################
# Source code
source("MadMax_RFDS_vsurf_predict_ESM.R")

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
spset <- as.numeric(command_args[1])
#spset <- 28475


########################
# Load occurrences

## species to model
# The csv is species name and its occurrence
sps.current <- read.csv('../SDM2024MAXE/Current_species.csv')
sptogo <- sps.current[spset, 1]
#

# load current occ
occ_cur <- readRDS('../SDM2024MAXE/sp_occr_1960-2020.Rds')            # 1960-2020
occ_cur <- occ_cur %>% filter(sps==sptogo)
pts.moll <- occ_cur[, 2:3]
#
########################
# load predictors
# soil and aridity
asoil <- list.files("../SDM2024MAXE/preds/soil", pattern = ".tif$", full.names = T)
asoil <- terra::rast(asoil)

# present
bios_pres <- list.files("../SDM2024MAXE/preds/bios_pres", pattern = ".tif$", full.names = T)
bios_pres <- terra::rast(bios_pres)

bios_pres <- c(bios_pres,asoil)

# projections variables to the future
# we choose ensemble average of 16 scenarios
gcms  <- c("ESM", "ACCESS", "MIROC6", "MPIESM")
ssps  <- c(126, 245, 370, 585)
# times <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
times <- c("2081-2100")
futvars <- expand.grid(gcms, ssps, times)
futvars$go <- paste(futvars[,1],futvars[,2],futvars[,3],sep="_")

bios_fut <- list()
for(i in 1:length(futvars$go)){
  bios_fut. <- list.files("../SDM2024MAXE/preds/bios_fut", 
                          pattern = paste0('_', futvars$go[i],".tif$"), 
                          full.names = T)
  bios_fut. <- terra::rast(bios_fut.)
  # fix names
  mypattern <- paste0("_",futvars$go[i])
  mypattern <- gsub("-",".",mypattern)
  names(bios_fut.) <- gsub(mypattern,"",names(bios_fut.))
  bios_fut. <- c(bios_fut.,asoil)
  
  bios_fut[[i]] <- bios_fut.
}
names(bios_fut) <- futvars$go
rm(asoil)

########################
# Run SDM

cat("\nRunning model for", sptogo,"\n")

#### Using doFuture

m <-try({
  
  MadMax(spname = sptogo,
         occ = pts.moll,
         fit.vars = bios_pres,
         proj.vars = bios_fut, # list with predictors for future scenarios
         proj.names = futvars$go, # vector names for proj.vars in the same order as proj.vars
         output_dir = "Results",
         check = F)
  
},silent = T)

if(class(m) == "try-error"){
  error_dir <- paste0("Error")
  if(!dir.exists(error_dir)){
    dir.create(error_dir, recursive = TRUE)
  }
  write.csv(data.frame(sps = sptogo,
                       error = m[1]),
            paste0(error_dir,"/",spset,".csv"))
  
  print(m)
}
