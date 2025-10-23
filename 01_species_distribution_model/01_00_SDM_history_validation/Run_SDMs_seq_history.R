rm(list=ls())
gc()

########################
list.of.packages <- c("raster", "terra", "blockCV", "enmSdmX", "randomForest", "VSURF",
                      "SDMtune", "dismo", "tidyverse", "parallel")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

sapply(list.of.packages, require, character.only = TRUE)

########################
# Source code
source("MadMax_history_vsurf.R")

########################
# Args
command_args <- commandArgs(trailingOnly = TRUE)
spset <- as.numeric(command_args[1])

########################
# Load occurrences

## species to model
# The RDS is just the name of these species
sps.history <- read.csv('History_species_1880-1920.csv')
sptogo <- sps.history[spset, 1]
# sptogo <- readRDS("sptorun.RDS")[spset]


print(Sys.time())

############################################################################################
# load current occ
occ_cur <- readRDS('sp_occr_1960-2020.RDS')            # 1960-2020
occ_cur <- occ_cur %>% filter(sps==sptogo)
pts.moll <- occ_cur[, 2:3]


# load historic occ
occ_his <- readRDS('sp_occr_1880-1920.RDS')
occ_his <- occ_his %>% filter(sps==sptogo)
occ_his <- occ_his[, 2:3]
############################################################################################
# load predictors

# soil and aridity
asoil <- list.files("preds/soil", pattern = ".tif$", full.names = T)
asoil <- terra::rast(asoil)

# present
bios_pres <- list.files("preds/bios_pres", pattern = ".tif$", full.names = T)
bios_pres <- terra::rast(bios_pres)

bios_pres <- c(bios_pres,asoil)

# projections variables to the future
# we choose ensemble average of 16 scenarios
gcms  <- c("1900")
futvars <- data.frame(go=gcms)

bios_fut <- list()
for(i in 1:length(futvars$go)){
  bios_fut. <- list.files("preds/bios_fut", 
                          pattern = paste0(futvars$go[i],".tif$"), 
                          full.names = T)
  bios_fut. <- terra::rast(bios_fut.)
  # fix names
  mypattern <- paste0("_",futvars$go[i])
  mypattern <- gsub("-",".",mypattern)
  names(bios_fut.) <- gsub(mypattern,"",names(bios_fut.))
  names(bios_fut.) <- c("bio1", "bio12", "bio15", "bio17", "bio19", "bio2", "bio4", "build", "crop", "moisture")
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
  
  MadMax_history(spname = sptogo,
         occ = pts.moll,
         occ.his = occ_his, 
         method = 'BRT',              # 'Maxent', "ANN", "BRT", "RF_dsamp"
         fit.vars = bios_pres,      
         proj.vars = bios_fut,        # list with predictors for future scenarios
         proj.names = futvars$go,     # vector names for proj.vars in the same order as proj.vars
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

print(Sys.time())
