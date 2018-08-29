#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")

# for linux
#setwd("/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/partial_code/a-basic-truncation-analysis/mean-median-across-all-conditions-analysis")

# for MAC
setwd("/Users/adnaniazi/Documents/phd/code/a-basic-truncation-analysis/mean-median-across-all-conditions-analysis")

# problems with execution on mac
#source('extract_read_data_rhdf5.R')
source('extract_read_data_hdf5r.R')
source('find_mean_median.R')


# Read all the filenames
bc08_fast5_list = list.files(path = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode08_all_reads/",
                             pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
bc09_fast5_list = list.files(path = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode09_all_reads/",
                             pattern = '*.fast5', full.names = TRUE, recursive = TRUE)
bc10_fast5_list = list.files(path = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode10_all_reads/",
                             pattern = '*.fast5', full.names = TRUE, recursive = TRUE)

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
        read_data <- extract_read_data_hdf5r(filepath)
        ans <- find_mean_median(read_data)
        df$bc09_mean[[i]] = ans$mean_norm_read_data
        df$bc09_median[[i]] = ans$median_norm_read_data  
    }, 
    error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
    })
}






