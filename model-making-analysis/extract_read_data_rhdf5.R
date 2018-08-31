extract_read_data_rhdf5 <- function(read_path){
    # extract raw data
    f5_tree <- rhdf5::h5ls(read_path)
    raw_reads <- f5_tree[(which(f5_tree == "/Raw/Reads") + 1), 1]
    raw_data <- rhdf5::h5read(read_path, raw_reads)$Signal
    ra <- rhdf5::h5readAttributes(read_path, raw_reads)
    ugk <- rhdf5::h5readAttributes(read_path, 'UniqueGlobalKey/channel_id/')

    tombo_normalization_key <- rhdf5::h5readAttributes(read_path, 'Analyses/RawGenomeCorrected_000/BaseCalled_template')
    scale <- tombo_normalization_key$scale
    shift <- tombo_normalization_key$shift

    read_id <- ra$read_id
    read_number <- ra$read_number
    start_time <- ra$start_time

    channel_number <- ugk$channel_number
    sampling_rate <- ugk$sampling_rate
    digitisation <- ugk$digitisation
    offset <- ugk$offset
    range <- ugk$range

    # pA normalization
    pa_norm_raw_data <- (raw_data + offset) * (range/digitisation)

    # extract events data
    tmp_path <- rhdf5::h5ls(read_path)[which(rhdf5::h5ls(read_path) == "/Analyses/Basecall_1D_000/BaseCalled_template")[1], 1]
    event_data <- rhdf5::h5read(read_path, tmp_path)$Events

    # tombo normalization
    tombo_norm_raw_data <- (raw_data - shift)/scale

    # time spent in each state as a function of number of samples
    model_state_length <- aggregate(event_data['length'], by = event_data['model_state'], sum)
    model_state_length$length <- as.character(model_state_length$length)
    # extract data into relevant columns
    event_data <- mutate(event_data, rdata1 = raw_data[event_data$start])
    event_data <- mutate(event_data, rdata2 = raw_data[event_data$start+1])
    event_data <- mutate(event_data, rdata3 = raw_data[event_data$start+2])
    event_data <- mutate(event_data, rdata4 = raw_data[event_data$start+3])
    event_data <- mutate(event_data, rdata5 = raw_data[event_data$start+4])

    event_data <- mutate(event_data, tdata1 = tombo_norm_raw_data[event_data$start])
    event_data <- mutate(event_data, tdata2 = tombo_norm_raw_data[event_data$start+1])
    event_data <- mutate(event_data, tdata3 = tombo_norm_raw_data[event_data$start+2])
    event_data <- mutate(event_data, tdata4 = tombo_norm_raw_data[event_data$start+3])
    event_data <- mutate(event_data, tdata5 = tombo_norm_raw_data[event_data$start+4])

    event_data <- mutate(event_data, pdata1 = pa_norm_raw_data[event_data$start])
    event_data <- mutate(event_data, pdata2 = pa_norm_raw_data[event_data$start+1])
    event_data <- mutate(event_data, pdata3 = pa_norm_raw_data[event_data$start+2])
    event_data <- mutate(event_data, pdata4 = pa_norm_raw_data[event_data$start+3])
    event_data <- mutate(event_data, pdata5 = pa_norm_raw_data[event_data$start+4])

    # combine raw, pa, and tombo normalized individual columns into respective merged columns
    # the merged data is in text form with each entry separated by _
    event_data <- event_data %>% unite(pdata, c(pdata1, pdata2, pdata3, pdata4, pdata5))
    event_data <- event_data %>% unite(rdata, c(rdata1, rdata2, rdata3, rdata4, rdata5))
    event_data <- event_data %>% unite(tdata, c(tdata1, tdata2, tdata3, tdata4, tdata5))

    # remove unnecessary columns
    event_data <- event_data %>% select(model_state, rdata, tdata, pdata)
    # aggregate event_data based on model state

    p <- function(..., sep='_') {
        paste(..., sep=sep, collapse=sep)
    }

    aggregated_event_rdata <- aggregate(event_data['rdata'], by=event_data['model_state'], p)
    aggregated_event_pdata <- aggregate(event_data['pdata'], by=event_data['model_state'], p)
    aggregated_event_tdata <- aggregate(event_data['tdata'], by=event_data['model_state'], p)

    # make a single dataframe
    final_df <- left_join(aggregated_event_rdata, aggregated_event_pdata, by='model_state')
    final_df <- left_join(final_df, aggregated_event_tdata, by='model_state')
    final_df <- left_join(final_df, model_state_length, by='model_state')

    rhdf5::H5close()
    return(final_df)
}
