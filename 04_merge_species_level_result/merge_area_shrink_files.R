# combine area shrink and shift direction
SDM.path <- '../../SDM2024MAXE/Realized/'
# SDM.path <- 'C:/Users/jnawang/Box/Biodiversity/SDM/Realized/'

files.area.shrink     <- list.files(path=paste0(SDM.path, 'area_shrink/'), pattern = '.csv$')
files.area.shrink     <- gsub('.csv$', '', files.area.shrink)


data.area.from.file <- function(file) {
  data <- read.csv(paste0(SDM.path, 'area_shrink/', file, '.csv'))
  data <- cbind(Species=file, data)
  return(data)
}

lst.data <- lapply(files.area.shrink, data.area.from.file)
data.area.shrink <- do.call(rbind, lst.data)
write.csv(data.area.shrink, 'area_shrink_merged_MAXE.csv', row.names = F)
