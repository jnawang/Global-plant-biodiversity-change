########a function to calculate distance between any two land grids on earth#######
##author: Junna Wang, 3/20/2025
##comment: saving the data will take lots of space, and reading it will take time, 
##comment: so making the data as small as possible
##use data table, rather than data frame to be efficient.
###################################################################################

library(raster)
library(terra)
library(data.table)
library(fst)

rm(list = ls())
gc()

######get command argument
command_args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(command_args[1])

# missed.chunks <- read.csv('missed_chunk_id.csv')
# id <- missed.chunks[id, 1]

print(id)
# id = 2005      # 15612.05, 2005, 6505

#################################################################################
#----------------------------a function to divide global grids into chunks-------
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
#----------------------------end of function

# use the function
r.template <- terra::rast("../SDM2024MAXE/preds/bios_pres/bio1.tif")
# r.template <- terra::rast("/Users/junnawang/UCDLab/Biodiversity/SDM/preds/bios_pres/bio1.tif")

block_size_row = 20
block_size_col = 50
chunk_grid_mapping <- grid_info(r.template, block_size_row, block_size_col)
setkey(chunk_grid_mapping, cell_id)

chunk_index <- sort(unique(chunk_grid_mapping$chunk_id))  # 2306 chunks in total

#################################################################################
# current working chunk_id
chunk_id_work <- chunk_index[id]

# get cell id of this working chunk
cell_id_work <- unique(chunk_grid_mapping$cell_id[chunk_grid_mapping$chunk_id==chunk_id_work])

######read raster files
r.barrier <- terra::rast('global_barrier.tif')           # cell values of barriers are 3 terra::rast
r.elev <- terra::rast('elev.tif')

r.empty <- r.template
terra::values(r.empty) <- NA

#################################################################################
######calculate grid distance to the origin grid
df.distance <- data.table(origin=integer(), target=integer(), distance=double())

time1 <- Sys.time()
for (icell in cell_id_work) {
  # set origin cell
  r.empty[icell] <- 1
#  plot(r.empty)
  
  b <- terra::buffer(r.empty, width=2200000, background=NA)  # 2500 km buffer
  # plot(b)
  b.trim <- terra::trim(b)
  # plot(b.trim)
  
  r.barrier.grid <- r.barrier
  r.barrier.grid[icell] <- 1
  r.barrier.grid <- terra::crop(r.barrier.grid, b.trim)
  
  # plot(r.barrier.grid)
  xy.origin <- terra::xyFromCell(r.template, icell)
  # points(xy.origin, col="red", pch=16, cex=1)
  
  grid.distance <- raster::gridDistance(raster::raster(r.barrier.grid), origin=1, omit=3)
  grid.distance <- terra::rast(grid.distance)
  # plot(grid.distance)

# get xy values of non-NA grids, distance < 2200 km, not themselves
  target.cell <- which(!is.na(terra::values(grid.distance)) & 
                         terra::values(grid.distance) < 1200000 & 
                         terra::values(grid.distance) > 0)
  
  # only if there are surrounding land cells    
  if (length(target.cell) > 0) {
    xy <- xyFromCell(grid.distance, target.cell)
    
    # ensure that these grids are all on lands
    values <- terra::extract(r.template, xy)
    
    # ensure that these grids are not at too-low elevation
    ele.target <- terra::extract(r.elev, xy)
    ele.origin <- terra::extract(r.elev, xy.origin)
    
    # get low land cell index
    id.low.land <- which(!is.na(values) & ele.origin$elev - ele.target$elev < 1200) # below 1200 m (maximum elevation migration)
    
    target.cell <- target.cell[id.low.land]      
    xy <- xy[id.low.land, , drop = FALSE]
    target.gid <- terra::cellFromXY(r.template, xy)
    
    # put it into the data table df.distance
    df <- data.table(origin=icell, 
                     target=target.gid, 
                     distance=terra::values(grid.distance)[target.cell])
    
    df.distance <- rbind(df.distance, na.omit(df))
  }

  # print(icell)
}
time2 <- Sys.time()
print(time2 - time1)

# check this df size
# object.size(df.distance)

result <- df.distance[, .(first_appearance = min(.I), last_appearance = max(.I)), by = origin]
write.csv(result, file=paste0('index/', id, '.csv'), row.names = F)


# try rds mode
# saveRDS(df.distance, paste0('chunks/', id, '.rds'))
# 
# system.time(readin <- readRDS(paste0('chunks/', id, '.rds')))

write_fst(df.distance[,2:3], paste0('chunks2/', id, '.fts'), compress = 0)

# system.time(readin <- read_fst(paste0('chunks/', id, '.fts')))

############################################below is the section, testing if my calculation above is correct##########################
# # reading this file is slow, linearly related to file size. 
# system.time(df_test <- readRDS(paste0('chunks/', id, '.rds')))
# 
# id.test <- unique(df_test$origin)[1]
# xy.origin <- xyFromCell(r.template, id.test)
# df_test_sub <- subset(df_test, origin==id.test)
# 
# r.empty.test <- r.empty
# terra::values(r.empty.test)[df_test_sub$target] <- df_test_sub$distance
# plot(r.empty.test)
# points(xy.origin, col="red", pch=16, cex=1)

