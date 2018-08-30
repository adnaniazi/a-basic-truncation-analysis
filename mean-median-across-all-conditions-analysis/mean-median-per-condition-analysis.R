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
ls<-foreach(filepath = bc08_fast5_list) %dopar% {
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

####################################
### FURTHER PROCESSING OF LISTS ###
####################################
save(ls, ls2, ls3, file = 'mean_median_data.RData')

bc08_mean <- unlist(sapply(ls, function(x) x[1]))
bc08_median <- unlist(sapply(ls, function(x) x[2]))
bc09_mean <- unlist(sapply(ls2, function(x) x[1]))
bc09_median <- unlist(sapply(ls2, function(x) x[2]))
bc10_mean <- unlist(sapply(ls3, function(x) x[1]))
bc10_median <- unlist(sapply(ls3, function(x) x[2]))

names(bc08_mean) <- NULL
names(bc09_mean) <- NULL
names(bc10_mean) <- NULL

names(bc08_median) <- NULL
names(bc09_median) <- NULL
names(bc10_median) <- NULL

# make a dataframe of all this data
maxlength <- max(length(bc08_mean), length(bc09_mean), length(bc10_mean))
length(bc08_mean) <- maxlength
length(bc09_mean) <- maxlength
length(bc10_mean) <- maxlength

length(bc08_median) <- maxlength
length(bc09_median) <- maxlength
length(bc10_median) <- maxlength

# replace all zeros by NAs
bc08_mean[bc08_mean==0] <- NA
bc09_mean[bc09_mean==0] <- NA
bc10_mean[bc10_mean==0] <- NA

bc08_median[bc08_median==0] <- NA
bc09_median[bc09_median==0] <- NA
bc10_median[bc10_median==0] <- NA

df = data.frame(bc08_mean, bc09_mean, bc10_mean,
                bc08_median, bc09_median, bc10_median)

####################################
############## PLOTTING OF DF ######
####################################
library(ggplot2)
df2 = stack(df)

ggplot(data = df2, aes(x = ind, y=values))+
    geom_violin()





