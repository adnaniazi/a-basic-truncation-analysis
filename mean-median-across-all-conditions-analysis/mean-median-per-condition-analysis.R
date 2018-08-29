#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
setwd("/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/partial_code/a-basic-truncation-analysis/mean-median-across-all-conditions-analysis")
source('extract_read_data_rhdf5.R')

bc08_fast5_list = list.files(path = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode08_all_reads/0",
                             pattern = '*.fast5', full.names = TRUE, recursive = TRUE)

