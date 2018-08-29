find_mean_median <- function(read_data){
    norm_read_data <- (read_data$raw - read_data$shift)/read_data$scale
    mean_norm_read_data <- mean(norm_read_data)
    median_norm_read_data <- median(norm_read_data)
    return(list(mean_norm_read_data=mean_norm_read_data, 
                median_norm_read_data=median_norm_read_data))
}
