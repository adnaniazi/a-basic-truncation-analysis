rm(list = ls())
mac_linux = Sys.info()[['sysname']]
#import packages
library(foreach)
library(doParallel)

if (mac_linux == 'Linux'){
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("rhdf5")
    setwd("/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/a-basic-truncation-analysis/mean-median-across-all-conditions-analysis")
    source('extract_read_data_rhdf5.R')
    bc08_fast5_list = list.files(path = "/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode08_all_reads/",
                                 pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
    bc09_fast5_list = list.files(path = "/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode09_all_reads/",
                                 pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
    bc10_fast5_list = list.files(path = "/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode10_all_reads/",
                                 pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
}

if (mac_linux == 'Darwin'){
    setwd("/Users/adnaniazi/Documents/phd/code/a-basic-truncation-analysis/mean-median-across-all-conditions-analysis")
    source('extract_read_data_hdf5r.R')
    bc08_fast5_list = list.files(path = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode08_all_reads/",
                                 pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
    bc09_fast5_list = list.files(path = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode09_all_reads/",
                                 pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
    bc10_fast5_list = list.files(path = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode10_all_reads/",
                                 pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
}


source('find_mean_median.R')
source('find_mean_median_parallel.R')

ls3 = find_mean_median_parallel(bc08_fast5_list, mac_linux)


# make an dataframe
N <- max(length(bc08_fast5_list), length(bc09_fast5_list), length(bc10_fast5_list))
df <- data.frame(bc08_mean = numeric(N), bc08_median = numeric(N),
                 bc09_mean = numeric(N), bc09_median = numeric(N),
                 bc10_mean = numeric(N), bc10_median = numeric(N))

# fill the data frame
i = 0
for (filepath in bc09_fast5_list){
    i = i + 1
    print(i)
    tryCatch({
        if (mac_linux == 'Linux'){
            read_data <- extract_read_data_rhdf5(filepath)
        }

        if (mac_linux == 'Darwin'){
            read_data <- extract_read_data_hdf5r(filepath)
        }

        ans <- find_mean_median(read_data)
        df$bc09_mean[[i]] = ans$mean_norm_read_data
        df$bc09_median[[i]] = ans$median_norm_read_data
    },
    error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
    })
}
##################
### BARCODE 08 ###
##################

#import packages
library(foreach)
library(doParallel)

# Calculate the number of cores
no_cores <- detectCores() - 20
# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#start time
strt<-Sys.time()

#loop
ls<-foreach(filepath = bc09_fast5_list) %dopar% {
    tryCatch({
        if (mac_linux == 'Linux'){
            read_data <- extract_read_data_rhdf5(filepath)
        }
        if (mac_linux == 'Darwin'){
            read_data <- extract_read_data_hdf5r(filepath)
        }

        ans <- find_mean_median(read_data)

    },
    error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")

        ans <- list(mean_norm_read_data = 0, median_norm_read_data = 0)
    })

}
stopCluster(cl)
ls

##################
### BARCODE 09 ###
##################

#import packages
library(foreach)
library(doParallel)

# Calculate the number of cores
no_cores <- detectCores() - 20
# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#start time
strt<-Sys.time()

#loop
ls2<-foreach(filepath = bc09_fast5_list) %dopar% {
    tryCatch({
        if (mac_linux == 'Linux'){
            read_data <- extract_read_data_rhdf5(filepath)
        }
        if (mac_linux == 'Darwin'){
            read_data <- extract_read_data_hdf5r(filepath)
        }

        ans <- find_mean_median(read_data)

    },
    error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")

        ans <- list(mean_norm_read_data = 0, median_norm_read_data = 0)
    })

}
stopCluster(cl)
ls2

##################
### BARCODE 10 ###
##################

# Calculate the number of cores
no_cores <- detectCores() - 20
# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)
#loop
ls3<-foreach(filepath = bc10_fast5_list) %dopar% {
    tryCatch({
        if (mac_linux == 'Linux'){
            read_data <- extract_read_data_rhdf5(filepath)
        }
        if (mac_linux == 'Darwin'){
            read_data <- extract_read_data_hdf5r(filepath)
        }

        ans <- find_mean_median(read_data)

    },
    error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")

        ans <- list(mean_norm_read_data = 0, median_norm_read_data = 0)
    })

}
stopCluster(cl)
ls3

