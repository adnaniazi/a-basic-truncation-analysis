import h5py

def get_fast5_read_data(fast5_file_path):
    """
    Read fast5 file attributes to find the start time of the read in seconds
    """
    # Make fast5 file object
    f = h5py.File(fast5_file_path)

    # Extract the read number from file name
    index_start = fast5_file_path.find('read_')
    index_end = fast5_file_path.find('_ch_')
    read_number = fast5_file_path[index_start+5:index_end]

    # Extract the read start time and duration in terms of samples
    group_name = '/Raw/Reads/Read_' + read_number
    group = f[group_name]
    start_time_in_samples = group.attrs.__getitem__('start_time')
    duration = group.attrs.__getitem__('duration')
    read_id = group.attrs.__getitem__('read_id')

    # Extract the sampling rate
    group_name = '/UniqueGlobalKey/channel_id'
    group = f[group_name]
    sampling_rate = group.attrs.__getitem__('sampling_rate')
    channel_id = group.attrs.__getitem__('channel_number')

    # Extract sequence length
    group_name = '/Analyses/Basecall_1D_000/Summary/basecall_1d_template'
    group = f[group_name]
    sequence_length = group.attrs.__getitem__('sequence_length')

    # Compute read start time in seconds
    read_start_time = start_time_in_samples/sampling_rate
    read_end_time = (start_time_in_samples + duration)/sampling_rate

    return read_start_time, read_end_time, sampling_rate, read_id, read_number, channel_id, sequence_length

