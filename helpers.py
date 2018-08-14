import os
import numpy as np
import itertools
from math import sqrt
import json
from ont_fast5_api import fast5_file as f5
from scipy import stats
from sys import platform
import pandas as pd
from random import shuffle
import shelve
import h5py
import matplotlib
if platform == "linux" or platform == "linux2":
    matplotlib.use('Agg')
elif platform == "darwin":
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import multiprocessing


def find_all_fast5s(directory):
    """ Recursively find Fast5 files within directory.
    Source: https://github.com/rrwick/Fast5-to-Fastq/blob/master/fast5_to_fastq.py

    :param directory: folder in which to recursively find fast5 files
    :return fast5s: list of path of fast5 file found recursively in the directory
    """
    fast5s = []
    for dir_name, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith('.fast5'):
                fast5s.append(os.path.join(dir_name, filename))
    return fast5s


def get_read_start_time(fast5_file_path):
    """
    Read fast5 file attributes to find the start time of the read in seconds
    """
    # Make fast5 file object
    f = h5py.File(fast5_file_path)

    # Extract the read number from file name
    index_start = fast5_file_path.find('read_')
    index_end = fast5_file_path.find('_ch_')
    read_id = fast5_file_path[index_start+5:index_end]

    # Extract the read start time and duration in terms of samples
    group_name = '/Raw/Reads/Read_' + read_id
    group = f[group_name]
    start_time_in_samples = group.attrs.__getitem__('start_time')
    duration = group.attrs.__getitem__('duration')

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

    return read_start_time, read_end_time, sampling_rate, read_id, channel_id, sequence_length


def make_channel_start_end_dict(file_paths_list):
    channel_start_end_dict = dict()
    channel_read_id_dict = dict()

    # Do for every file in the list of files provided
    count = 0
    num_files = len(file_paths_list)
    for fast5_file_path in file_paths_list:
        count += 1
        read_start_time, read_end_time, _, read_id, channel_id, sequence_length = get_read_start_time(fast5_file_path)

        key_start_end = 'Ch-' + str(int(channel_id))
        value_start_end = [read_start_time, read_end_time]
        value_read_id = read_id

        try:
            channel_start_end_dict[key_start_end].append(value_start_end)
            channel_read_id_dict[key_start_end].append(value_read_id)

        except:
            channel_start_end_dict.update({key_start_end:[value_start_end]})
            channel_read_id_dict.update({key_start_end:[value_read_id]})

        print("Processing file {} of {}".format(count, num_files))

    print(channel_start_end_dict)
    print(channel_read_id_dict)

    return channel_start_end_dict, channel_read_id_dict


def extract_star_time_and_read_length(fast5_file_path):
    extract_star_time_and_read_length.count +=1
    # Read fast5 file
    f = h5py.File(fast5_file_path)

    # Extract the read number from file name
    index_start = fast5_file_path.find('read_')
    index_end = fast5_file_path.find('_ch_')
    read_id = fast5_file_path[index_start+5:index_end]

    # Extract the read start time and duration in terms of samples
    group_name = '/Raw/Reads/Read_' + read_id
    group = f[group_name]
    start_time_in_samples = group.attrs.__getitem__('start_time')

    group_name = '/UniqueGlobalKey/channel_id'
    group = f[group_name]
    sampling_rate = group.attrs.__getitem__('sampling_rate')

    # Extract sequence length
    group_name = '/Analyses/Basecall_1D_000/Summary/basecall_1d_template'
    group = f[group_name]
    sequence_length = group.attrs.__getitem__('sequence_length')

    read_start_time_in_seconds = start_time_in_samples/sampling_rate
    print(f"Start time: {read_start_time_in_seconds}, Length: {sequence_length}, Read: {extract_star_time_and_read_length.count}")

    return read_start_time_in_seconds, sequence_length


def cummulative_number_of_reads_over_time(file_paths_list, chunksize=200):
    extract_star_time_and_read_length.count = 0

    # get read start time, make a vector
    pool = multiprocessing.Pool()

    result = pool.map(extract_star_time_and_read_length, file_paths_list, chunksize=chunksize)
    result = np.array(result)
    pool.close()
    pool.join()

    # sort the array by time (0th column)
    result = result[result[:, 0].argsort()]
    print('\nSorted Array:', result)

    df = pd.DataFrame(result)

    # Add a read count column
    df['counter'] = range(len(df))

    # Rename all columns
    df.columns = ['read_start_time_in_seconds', 'read_length_in_bp', 'cumm_read_count']

    # Add cummulative read count column
    df['cumm_read_lenth_in_bp'] = df.read_length_in_bp.cumsum()

    # return format ['read_start_time_in_seconds', 'read_length_in_bp', 'cumm_read_count', 'cumm_read_lenth_in_bp']
    return df








    # arrange the vector
    # at every start time, increment the count


def median_normalize(raw_data):
    """ This function median normalizes the data. See Nanoraw paper by Markus Stoiber.

    :param raw_data: raw non pA-normalized data from the Fast5
    :return median_normalized_data: data that has been median normalized
    """

    # Calculate median absolute deviation (MAD)
    read_median = np.median(raw_data)
    #print("Read median:", read_median)

    mad = np.median(abs(raw_data - read_median))
    #print("MAD:", mad)

    median_normalized_data = (raw_data - read_median)/mad

    return median_normalized_data


def pA_normalize(raw_data, range, digitisation, offset):

    """ This function does the pA normalization on raw data from ONT's Fast5 files.
        raw unit = range/digitisation
        pA_normalized_data = raw_unit*(raw_data + offset)

    :param raw_data: raw non pA-normalized data from the Fast5
    :param range: range of ADC (given in Fast5 file)
    :param digitisation: levels of quantization in ADC (given in Fast5 file)
    :param offset: perhaps its the offset of ADC (give in Fast5 file)
    :return pA_normalized_data: data that has been median normalized
    """

    raw_unit = range/digitisation
    pA_normalized_data = raw_unit * (raw_data + offset)

    return pA_normalized_data


def save_filepaths_in_json(file_paths, directory_path='.', file_name='untitled.json'):
    """ This function save the file paths in JSON format

    """
    if directory_path == '.':
        directory_path = os.getcwd()

    with open(os.path.join(directory_path, file_name), 'w') as fp:
        json.dump(file_paths, fp, indent=0)

    print(file_name, 'saved!')


def load_filepaths_from_json(file_name, directory_path='.'):
    """
    This function loads the JSON file containing file paths and returns it as a list.

    """
    if directory_path == '.':
        directory_path = os.getcwd()

    # load json file
    with open(os.path.join(directory_path, file_name), 'r') as fp:
        paths_list = json.load(fp)
    return paths_list