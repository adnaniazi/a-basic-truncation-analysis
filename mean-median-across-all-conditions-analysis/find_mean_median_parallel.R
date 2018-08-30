#source('find_mean_median.R')

find_mean_median_parallel <- function(filepath_list, mac_linux){
    library(foreach)
    library(doParallel)
    # Calculate the number of cores
    no_cores <- detectCores() - 10
    # Initiate cluster
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)

    #start time
    strt<-Sys.time()

    #loop
    ls2<-foreach(filepath = filepath_list) %dopar% {
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
            # return zero mean and median if the shift and scale data could not be found in the Fast5 data
            ans <- list(mean_norm_read_data = 0, median_norm_read_data = 0)
        })

    }
    stopCluster(cl)
    return(ls2)
}