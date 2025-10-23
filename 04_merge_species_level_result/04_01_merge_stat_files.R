rm(list = ls())
# combine area shrink and shift direction
SDM.path <- '../../SDM2024MAXE/Results/stats/'
# SDM.path <- 'C:/Users/jnawang/Downloads/stats/'

# merge performance data into a large file
files  <- list.files(path=SDM.path, pattern = '.csv$', full.names = TRUE)


data.perform.from.file <- function(file) {
  data <- read.csv(file)
  return(data[1,])
}

lst.data     <- lapply(files,  data.perform.from.file)
data.perform <- do.call(rbind, lst.data)

# only choose AUC_test > 0.75, use AUC_test or AUC?
# data.perform <- data.perform[data.perform$AUC_test>=0.75,]

write.csv(data.perform, 'SDM_perform_merged.csv', row.names = F)
################


# merge relative importance into another file
#file <- files[1]
data.ri.from.file <- function(file) {
  data <- read.csv(file)
  n  <- nrow(data)
  ri <- data[3:n, 3][order(data[3:n, 1])]
  out <- data.frame(t(c(data[1,1], data[1,3], ri)))
  colnames(out) <- c("sps", "AUC", sort(data[3:n, 1]))
  return(out)     # 
}

lst.data <- lapply(files,  data.ri.from.file)
data.ri  <- do.call(rbind, lst.data)

# only choose AUC_test > 0.75, use AUC_test or AUC?
# data.ri <- data.ri[data.ri$AUC_test>=0.75,]
write.csv(data.ri, 'SDM_ri_merged.csv', row.names = F)
