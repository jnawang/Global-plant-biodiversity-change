library(terra)
library(picante)

# setwd("D:/rSpace/beta_diversity")
# setwd('C:/Users/bnulr/Desktop/JunnaWang/beta_diversity/')
setwd('/home/jnawang/SDM2024MAXE_RFDS_process/novel_community/map2bin/')

read.filerange <- function(from, to) {
  data <- list()
  for (filenum in from:to) {
    read.file <- file(paste(filenum, "present.dat", sep = ""), "rb")      # present.dat  Realized_ESM_245_2081-2100.dat  Junna
    two_int <- readBin(read.file, integer(), n = 2)
    batch <- two_int[1] # 6400
    for (i in 1:batch) {
      two_int <- readBin(read.file, integer(), n = 2)
      # print(c(filenum, two_int))
      cel_idx <- two_int[1]
      num <- two_int[2]
      if (num > 0) {
        idices <- readBin(read.file, integer(), n = num)
        data[[cel_idx]] <- idices
      } else {
        data[[cel_idx]] <- {0}
      }
    }
    close(read.file)
  }
  return(data)
}


caculate.file <- function(ifile = 177) {
#  ifile <- 177
  # read phylogenetic tree
  # tree   <- readRDS("C:/Users/jnawang/Box/Biodiversity/pre_post_processing/plants_phy.RDS")
  tree   <- readRDS("/home/jnawang/SDM2024MAXE_RFDS_process/richness/traits/plants_phy.RDS")
  # change tree label
  tree$tip.label <- gsub("_"," ", tree$tip.label)
  # read species list  
  # sps_id <- read.csv('C:/Users/jnawang/Box/Biodiversity/pre_post_processing/traits_phylo/model_sps_id.csv')
  sps_id <- read.csv('/home/jnawang/SDM2024MAXE_RFDS_process/richness/traits/model_sps_id.csv')
  #
  bio1 <- terra::values(rast("/home/jnawang/SDM2024MAXE/preds/bios_pres/bio1.tif"))
  # bio1  <- terra::values(terra::rast("C:/Users/jnawang/Box/Biodiversity/SDM/preds/bios_pres/bio1.tif"))
  landcells <- which(!is.na(bio1))
  ncell.file <- 6400
  id.start   <- (ifile - 1) * ncell.file + 1
  id.end     <- min(ifile*ncell.file, length(landcells))
  cells.file <- landcells[id.start:id.end]
  # read species in these grids
  data <- read.filerange(ifile, ifile)        # this is a list. 
  # output
  output <- data.frame(GID=cells.file, pd=0.0)
  # I need to figure out the index
  for (i in 1:length(cells.file)) {
    i.cel <- cells.file[i]
    if (i%%100==0) {print(i)}
      sps_cell <- data[[i.cel]]
      # calculate PD for each grid
      if (length(sps_cell) > 1 | length(sps_cell) == 1 & sps_cell[1] > 0) {
        sample <- matrix(1, nrow=1, ncol=length(sps_cell))
        colnames(sample) <- sps_id$sps[sps_cell]
        subtree <- prune.sample(sample,tree)
        output$pd[i] <- pd(sample, subtree)$PD[1]  # this is too slow. 
      }
  }
#  write.csv(output, paste("C:/Users/jnawang/Downloads/pd-", ifile, "present.csv", sep = ""), row.names = FALSE)  # present.csv  Realized_ESM_245_2081-2100.csv Junna
  write.csv(output, paste("/home/jnawang/SDM2024MAXE_RFDS_process/richness/phylogenetic/pd-", ifile, "present.csv", sep = ""), row.names = FALSE)
}

# caculate.file(1)
#

# in total 244 files
command_args <- commandArgs(trailingOnly = TRUE)
id.f <- as.numeric(command_args[1])
caculate.file(id.f)

# what's the range of beta diversity?
