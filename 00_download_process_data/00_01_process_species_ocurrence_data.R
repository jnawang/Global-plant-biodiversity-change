library(terra)
library('tidyverse')
library('ggplot2')
library('lubridate')
rm(list=ls())
#
setwd('/Users/junnawang/UCDLab/Biodiversity/data/occurrence')
#
occ_all <- readRDS('sp_occr_records_all.Rds')
occ_all <- occ_all %>% filter(!is.na(date_collected)) %>% 
  distinct(.keep_all = TRUE)  # remove duplicate points
#
# add moll coordinate
vect_coords <- terra::vect(cbind(occ_all[,3], occ_all[,2]), crs = "EPSG:4326")

# Transform to Mollweide projection (EPSG:54009)
vect_mollweide <- terra::project(vect_coords, "+proj=moll +datum=WGS84")
occ_all[,6:7] <- crds(vect_mollweide)
colnames(occ_all)[c(1, 6:7)] <- c('sps', 'x', 'y')
occ_all_write <- occ_all %>% distinct(sps, cell, .keep_all = TRUE)
saveRDS(occ_all_write[,c(1, 6:7)], file = "sp_occr_all.Rds")
#
#
occ_old <- occ_all %>% filter(date_collected < as.Date('1920-12-31') & date_collected > as.Date('1880-1-1'))  # I try to use this data
#
# thin the data that are in one grid; easy, we get cell id for each location and then remove duplicates
occ_old <- occ_old %>% distinct(sps, cell, .keep_all = TRUE)
#
# find species with > 20 data
occ_num_sps <- occ_old %>% group_by(sps) %>% summarise(nocc=n_distinct(cell)) %>% arrange(desc(nocc)) %>% filter(nocc>20)
#
write.csv(occ_num_sps, 'History_species.csv', row.names = F)
#
# find species in our data base
# species_modeled <- read.csv('/Users/junnawang/UCDLab/Biodiversity/pre_post_processing/traits_phylo/model_sps_id.csv')
species_modeled <- readRDS('plants_list.RDS')
occ_num_sps <- occ_num_sps %>% filter(sps %in% species_modeled) # 5557 species
occ_old <- occ_old %>% filter(sps %in% occ_num_sps$sps)
#
# write the occ in 1880-1920
saveRDS(occ_old[,c(1, 6:7)], file = "sp_occr_1880-1920.Rds")
#
##################
# get current data
occ_cur <- occ_all %>% filter(date_collected >= as.Date('1960-1-1') & date_collected <= as.Date('2024-1-1'))  # I try to use this data
#
# thin the data that are in one grid; easy, we get cell id for each location and then remove duplicates
occ_cur <- occ_cur %>% distinct(sps, cell, .keep_all = TRUE)
occ_num_sps <- occ_cur %>% group_by(sps) %>% summarise(nocc=n_distinct(cell)) %>% filter(nocc>=10)  #  %>% arrange(desc(nocc))
occ_num_sps <- occ_num_sps %>% filter(sps %in% species_modeled)
#
write.csv(occ_num_sps, 'Current_species.csv', row.names = F)
#
occ_cur <- occ_cur %>% filter(sps %in% occ_num_sps$sps)
#
# histogram of species occurrence data
# hist(occ_num_sps$nocc, bin=100)
#
# write the occ in 1970-2020
saveRDS(occ_cur[,c(1, 6:7)], file = "sp_occr_1960-2020.Rds")




