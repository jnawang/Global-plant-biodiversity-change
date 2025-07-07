library(terra)
library(betapart)

# setwd("D:/rSpace/beta_diversity")
# setwd('C:/Users/bnulr/Desktop/JunnaWang/beta_diversity/')
setwd('/home/jnawang/SDM2024MAXE_RFDS_process/novel_community/map2bin/')        # use absolute dir

read.filerange <- function(from, to) {
    data <- list()
    for (filenum in from:to) {
        read.file <- file(paste(filenum, "Realized_ESM_245_2081-2100.dat", sep = ""), "rb")      # Realized_ESM_245_2081-2100.dat   present.dat  Junna
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
            }
        }
        close(read.file)
    }
    return(data)
}

read.cel_indices <- function(filenum) {
    data <- list()
    read.file <- file(paste(filenum, "Realized_ESM_245_2081-2100.dat", sep = ""), "rb")   # Realized_ESM_245_2081-2100.dat     present.dat   Junna
    two_int <- readBin(read.file, integer(), n = 2)
    batch <- two_int[1] # 6400
    for (i in 1:batch) {
        two_int <- readBin(read.file, integer(), n = 2)
        cel_idx <- two_int[1]
        num <- two_int[2]
        data <- c(data, cel_idx)
        if (num > 0) {
            readBin(read.file, integer(), n = num)
        }
    }
    close(read.file)
    return(unlist(data))
}

# abstract.fstidx <- function() {
#     bin_files <- list.files(pattern = "\\.dat$")
#     # Create an empty list to store results
#     result_list <- list()
# 
#     # Loop through each file
#     for (bin_file in bin_files) {
#         num <- as.integer(sub("present\\.dat", "", bin_file))
#         third_integer <- readBin(file(bin_file, "rb"), integer(), n = 3)[3]
#         result_list[[as.character(num)]] <- third_integer
#     }
#     result_df <- data.frame(ifile = as.integer(names(result_list)), fstidx = unlist(result_list))
# 
#     # Write the result to a CSV filethird_value
#     write.csv(result_df, "fstidx.csv", row.names = FALSE)
# 
#     # Print the result
#     print(result_df)
# }
# abstract.fstidx()


find_boundary_file <- function(data_frame, N) {
    # Find the last row where col2 is less than or equal to N
    last_row_index <- max(which(data_frame$fstidx <= N))

    # Return the corresponding col1 value
    if (last_row_index > 0) {
        return(data_frame$ifile[last_row_index])
    } else {
        return(1)
    }
}

conert2matrix <- function(lst) {
    distinct_values <- sort(unique(unlist(lst)))
    list_size <- length(lst)
    m <- matrix(0, nrow = list_size, ncol = length(distinct_values))
    for (i in 1:list_size) {
        m[i, match(lst[[i]], distinct_values)] <- 1
    }
    colnames(m) <- distinct_values
    return(m)
}

# for example, working on calculating beta diversity for the cells in file 10
# decide how many files we need to read in
# find the first cell and the last cell
caculate.file <- function(ifile = 33) {
    bio1 <- terra::values(rast("/home/jnawang/SDM2024MAXE/preds/bios_pres/bio1.tif"))
    # dimensions  : 1761, 4333, 1 (row, col)
    landcells <- which(!is.na(bio1))
    nlandcell <- length(landcells)
    nrol <- 1761
    ncol <- 4333
    ncell.file <- 6400
    nneighbor <- 3
    threshold <- 25
    #
    nfile  <- ceiling(length(landcells)/ncell.file)
    fstidx <- data.frame(ifile=1:nfile, fstidx = landcells[1+(0:(nfile-1))*ncell.file])
    print(fstidx)
    #
    cell.min <- landcells[(ifile - 1) * ncell.file + 1] - nneighbor - nneighbor * ncol
    cell.max <- landcells[min(ifile * ncell.file + 1, nlandcell)] + nneighbor + nneighbor * ncol    
    
    # read these files into a large integer list with index as cell.id
    ############################# your code here.

    data <- read.filerange(find_boundary_file(fstidx, cell.min), find_boundary_file(fstidx, cell.max))
    object.size(data)
    # loop through every cell in the file 10
    # for each cell, find its 11*11 neighbor cells
    # example code
    result <- data.frame()
    cel_indices <- read.cel_indices(ifile)
    for (i.cel in 1:length(cel_indices)) {
        focal.cell <- cel_indices[i.cel]
        cells.region <- vector()
        for (i in -nneighbor:nneighbor) {
            cells.region <- c(cells.region, (focal.cell - nneighbor):(focal.cell + nneighbor) + ncol * i)
        }
        # remove cells below the last land cell
        cells.region <- cells.region[cells.region <= landcells[nlandcell]]
        
        # print(cells.region)

        valid_cells <- cells.region[!sapply(cells.region, function(idx) {
            return(is.null(data[[idx]]))
        })]
        if (length(valid_cells) < threshold) {
            result <- rbind(result, data.frame(
                cell = focal.cell, validnum = length(valid_cells),
                # JTU = NA, JNE = NA, JAC = NA    # for Jaccard method
                SIM = NA, SNE = NA, SOR = NA
            ))
        } else {
            valid_cells.data <- lapply(valid_cells, function(x) {
                data[[x]]
            })
            if (length(valid_cells) == threshold) {
                valid_cells.m <- conert2matrix(valid_cells.data)
                # valid_cells.beta <- beta.multi(valid_cells.m, index.family = "jaccard")  # for Jaccard method
                valid_cells.beta <- beta.multi(valid_cells.m, index.family = "sorensen")
            } else {
                mat <- sapply(1:100, function(x) {
                    sample.data <- valid_cells.data[sample(seq_len(length(valid_cells)), threshold)]
                    valid_cells.m <- conert2matrix(sample.data)
                    # unlist(beta.multi(valid_cells.m, index.family = "jaccard"))   # for Jaccard method
                    unlist(beta.multi(valid_cells.m, index.family = "sorensen"))
                })
                mat <- t(mat)
                df <- as.data.frame(mat)
                valid_cells.beta <- as.list(colMeans(df))
            }
            result <- rbind(result, data.frame(
                cell = focal.cell, validnum = length(valid_cells),
                # JTU = valid_cells.beta$beta.JTU, JNE = valid_cells.beta$beta.JNE, JAC = valid_cells.beta$beta.JAC  # for Jaccard method
                SIM = valid_cells.beta$beta.SIM, SNE = valid_cells.beta$beta.SNE, SOR = valid_cells.beta$beta.SOR
            ))
            # print(paste(i.cel, focal.cell, length(valid_cells), valid_cells.beta$beta.JAC))  # for Jaccard method
            print(paste(i.cel, focal.cell, length(valid_cells), valid_cells.beta$beta.SOR))
        }
        # if (i.cel > 100) {
        #     break
        # }
    }
    write.csv(result, paste("/home/jnawang/SDM2024MAXE_RFDS_process/richness/beta/beta-", ifile, "Realized_ESM_245_2081-2100.csv", sep = ""), row.names = FALSE)     # Realized_ESM_245_2081-2100.dat     present.dat   Junna  
}

# in total 244 files
command_args <- commandArgs(trailingOnly = TRUE)
id.f <- as.numeric(command_args[1])
caculate.file(id.f)

