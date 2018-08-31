rm(list = ls())
mac_linux = Sys.info()[['sysname']]
#import packages
library(foreach)
library(doParallel)
library(rhdf5)
library(dplyr)
library(tidyverse)

if (mac_linux == 'Linux'){
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("rhdf5")
    source('extract_read_data_rhdf5.R')
    bc08_fast5_list = list.files(path = "/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode08_all_reads/",
                                 pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
    bc09_fast5_list = list.files(path = "/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode09_all_reads/",
                                 pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
    bc10_fast5_list = list.files(path = "/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode10_all_reads/",
                                 pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
    setwd('/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/a-basic-truncation-analysis/model-making-analysis')
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
    setwd('/Users/adnaniazi/Documents/phd/code/a-basic-truncation-analysis/model-making-analysis')
}

source('extract_read_data_rhdf5.R')

min_reads = min(length(bc08_fast5_list), length(bc09_fast5_list), length(bc10_fast5_list))
bc08_fast5_list_truncated <- bc08_fast5_list[1:min_reads]
bc09_fast5_list_truncated <- bc09_fast5_list[1:min_reads]
bc10_fast5_list_truncated <- bc10_fast5_list[1:min_reads]


## barcode-09 processing
########################

# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#start time
strt<-Sys.time()

#loop
ls<-foreach(filepath = bc10_fast5_list) %dopar% {
    tryCatch({
        if (mac_linux == 'Linux'){
            read_df <- extract_read_data_rhdf5(filepath)
        }
        if (mac_linux == 'Darwin'){
            read_df <- extract_read_data_hdf5r(filepath)
        }

        read_df

    },
    error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")

        read_df <- NA
    })

}
stopCluster(cl)
ls










extract_read_data_rhdf5('/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode10_all_reads/0/sars_HP_Z240_Tower_Workstation_20180710_FAH87711_MN21607_sequencing_run_197_lig_79579_read_11_ch_46_strand.fast5')


final_df$length = as.character(final_df$length)
final_df_copy$length = as.character(final_df_copy$length)


c <- merge(final_df, final_df_copy, by = "model_state", all = TRUE)
c <- c %>% unite(pdata_combined, c(pdata.x, pdata.y))
c <- c %>% unite(tdata_combined, c(tdata.x, tdata.y))
c <- c %>% unite(rdata_combined, c(rdata.x, rdata.y))
c <- c %>% unite(length_combined, c(length.x, length.y))







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
