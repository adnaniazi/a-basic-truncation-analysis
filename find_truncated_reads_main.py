import h5py
import pandas as pd
import os
import os.path


def _print_attrs(name, obj):
    """
    Recursively print the f5 tree structure.

    :param name:
    :param obj:
    :return:
    """
    print(name)
    for key, val in obj.attrs.items():
        print("    %s: %s" % (key, val))


def get_fast5_read_data(fast5_file_path):
    """
    Read fast5 file attributes to find the start time of the read in seconds
    """
    # Make fast5 file object
    f = h5py.File(fast5_file_path, 'r')

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

    data = {'read_start_time': read_start_time,
            'read_end_time':  read_end_time,
            'sampling_rate': sampling_rate,
            'read_id': read_id.decode('UTF-8'),
            'read_number': int(read_number),
            'channel_id': int(channel_id),
            'fast5_sequence_length': sequence_length
            }

    # options for traversing the f5 tree
    # f.visititems(_print_attrs)

    return data


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


if __name__ == '__main__':

    TESTMODE = False
    REMOVE_DATAFRAME_FILES = False

    # ------------------------------------------------------ #
    # --------------- FAST5 DATA PROCESSING ---------------- #
    # ------------------------------------------------------ #

    if REMOVE_DATAFRAME_FILES:
        os.system('rm uncl_fast5_data_msgpack_df')
        os.system('rm bc10_fast5_data_msgpack_df')

    # read directories
    bc10_fast5_dir = '/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/' \
                     '20180710_1431_197_lig/backbone_aligned_reads/pass/barcode10_all_reads'
    uncl_fast5_dir = '/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/' \
                     '20180710_1431_197_lig/backbone_aligned_reads/pass/unclassified_all_reads'

    # if fast5 dataframes don't exist, then read the fast5, make dataframes, and save them in the current dir
    if not(os.path.isfile('bc10_fast5_data_msgpack_df')) or not(os.path.isfile('uncl_fast5_data_msgpack_df')):

        # get fast5 files list
        bc10_read_paths = find_all_fast5s(bc10_fast5_dir)
        uncl_read_paths = find_all_fast5s(uncl_fast5_dir)

        # read fast5 files
        bc10_df = pd.DataFrame(columns=['read_start_time',
                                        'read_end_time',
                                        'sampling_rate',
                                        'read_id',
                                        'read_number',
                                        'channel_id',
                                        'fast5_sequence_length'])
        uncl_df = bc10_df

        for i, filepath in enumerate(bc10_read_paths):
            print("Processing file {0} of {1} in Barcode 10 class".format(i, len(bc10_read_paths)))
            read_data = get_fast5_read_data(filepath)
            read_data['class'] = 'barcode10'
            bc10_df = bc10_df.append(read_data, ignore_index=True)
            if TESTMODE:
                if i == 100:
                    break

        bc10_df['class'] = bc10_df['class'].astype('category')
        bc10_df.to_msgpack('bc10_fast5_data_msgpack_df')

        for i, filepath in enumerate(uncl_read_paths):
            print("Processing file {0} of {1} in unclassified class".format(i, len(uncl_read_paths)))
            read_data = get_fast5_read_data(filepath)
            read_data['class'] = 'unclassified'
            uncl_df = uncl_df.append(read_data, ignore_index=True)
            if TESTMODE:
                if i == 100:
                    break

        uncl_df['class'] = uncl_df['class'].astype('category')
        uncl_df.to_msgpack('uncl_fast5_data_msgpack_df')

    # ------------------------------------------------------ #
    # -------------- SAMFILE DATA PROCESSING --------------- #
    # ------------------------------------------------------ #
    BC10_SAM_FILE_PATH = os.path.join('/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/alignment_to_plasmid/barcode10_ecorv_backbone_only/aln.sam')
    UNCL_SAM_FILE_PATH = os.path.join('/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/MinION/20180710_1431_197_lig/alignment_to_plasmid/unclassified_ecorv_backbone_only/aln.sam')
    SAM_TESTMODE = False
    REMOVE_SAM_DATAFRAME_FILES = False

    if REMOVE_SAM_DATAFRAME_FILES:
        os.system('rm uncl_sam_data_msgpack_df')
        os.system('rm bc10_sam_data_msgpack_df')

    if not(os.path.isfile('bc10_sam_data_msgpack_df')) or not(os.path.isfile('uncl_sam_data_msgpack_df')):
        with open(BC10_SAM_FILE_PATH) as f:
            # First two lines contain the SAM header, so ignore them
            next(f)
            next(f)
            bc10_sam_data = {}
            bc10_sam_df = pd.DataFrame(columns=['strand',
                                                'read_id',
                                                'alignment_start_position',
                                                'mapping_quality',
                                                'sam_data_class',
                                                'alignment_length'])
            for i, line in enumerate(f):
                print('Processing record {0} in bacorde10 alignment SAM file'.format(i))
                sp = line.split()
                # filter forward and reverse mapped strands
                if int(sp[1]) == 0:             # forward strand or + strand
                    bc10_sam_data['strand'] = 'positive'
                elif int(sp[1]) == 16:          # reverse strand or - strand
                    bc10_sam_data['strand'] = 'negative'
                bc10_sam_data['read_id'] = sp[0]
                bc10_sam_data['alignment_start_position'] = sp[3]
                bc10_sam_data['mapping_quality'] = sp[4]
                bc10_sam_data['sam_data_class'] = 'barcode10'
                bc10_sam_data['alignment_length'] = len(sp[9])
                bc10_sam_df = bc10_sam_df.append(bc10_sam_data, ignore_index=True)

                if SAM_TESTMODE:
                    if i == 100:
                        break

        bc10_sam_df['sam_data_class'] = bc10_sam_df['sam_data_class'].astype('category')
        bc10_sam_df.to_msgpack('bc10_sam_data_msgpack_df')

        with open(UNCL_SAM_FILE_PATH) as f:
            # First two lines contain the SAM header, so ignore them
            next(f)
            next(f)
            uncl_sam_data = {}
            uncl_sam_df = pd.DataFrame(columns=['strand',
                                                'read_id',
                                                'alignment_start_position',
                                                'mapping_quality',
                                                'sam_data_class',
                                                'alignment_length'])
            for i, line in enumerate(f):
                print('Processing record {0} in unclassified alignment SAM file'.format(i))
                sp = line.split()
                # filter forward and reverse mapped strands
                if int(sp[1]) == 0:             # forward strand or + strand
                    uncl_sam_data['strand'] = 'positive'
                elif int(sp[1]) == 16:          # reverse strand or - strand
                    uncl_sam_data['strand'] = 'negative'
                uncl_sam_data['read_id'] = sp[0]
                uncl_sam_data['alignment_start_position'] = sp[3]
                uncl_sam_data['mapping_quality'] = sp[4]
                uncl_sam_data['sam_data_class'] = 'unclassified'
                uncl_sam_data['alignment_length'] = len(sp[9])
                uncl_sam_df = uncl_sam_df.append(uncl_sam_data, ignore_index=True)

                if SAM_TESTMODE:
                    if i == 100:
                        break

        uncl_sam_df['sam_data_class'] = uncl_sam_df['sam_data_class'].astype('category')
        uncl_sam_df.to_msgpack('uncl_sam_data_msgpack_df')

    # ------------------------------------------------------ #
    # ---------------- MERGING DATA FRAMES ----------------- #
    # ------------------------------------------------------ #

    # First load data frames from the saved files
    bc10_fast5_df = pd.read_msgpack('bc10_fast5_data_msgpack_df')
    uncl_fast5_df = pd.read_msgpack('uncl_fast5_data_msgpack_df')
    bc10_sam_df = pd.read_msgpack('bc10_sam_data_msgpack_df')
    uncl_sam_df = pd.read_msgpack('uncl_sam_data_msgpack_df')

    # merge SAM and FAST5 data frames on the basis of read_ids
    bc10_df = pd.merge(bc10_fast5_df, bc10_sam_df, on='read_id')
    uncl_df = pd.merge(uncl_fast5_df, uncl_sam_df, on='read_id')

    # concatenate the two data frames
    data = pd.concat([bc10_df, uncl_df])

    # MAIN ALGORITHM FOR DETECTING TRUNCATED READS
    # 1. divide data into positive and negative strands
    pos_strand_data = data[data.strand == 'positive']
    neg_strand_data = data[data.strand == 'negative']

    # 2. order the data rows by start time
    pos_strand_data = pos_strand_data.sort_values('read_start_time')
    neg_strand_data = neg_strand_data.sort_values('read_start_time')

    # 3. Get channel wise data
    gap_threshold = 4

    channel_wise_dict_pos_strand = {}
    channel_wise_dict_neg_strand = {}

    for ch in range(1, 512+1):
        print('processing channel {0} of 512 channels (+ strands)'.format(ch))
        channel_wise_data_df_pos_strand = pd.DataFrame(columns=['read_start_time', 'read_end_time', 'sampling_rate',
                                                                'read_id', 'read_number', 'channel_id',
                                                                'fast5_sequence_length', 'class', 'strand',
                                                                'alignment_start_position', 'mapping_quality',
                                                                'sam_data_class', 'alignment_length'])

        pos_strand_ch_wise_data = pos_strand_data[pos_strand_data.channel_id == ch].sort_values('read_start_time')
        pos_strand_ch_wise_data_numrows = pos_strand_ch_wise_data.shape[0]

        for i in range(pos_strand_ch_wise_data_numrows-1):
            # check for two things
            # 1. is alignment of the two fragments back-to-back in time
            # 2. is the start position of the end fragment after the end position of the first fragment in terms of alignment
            if (abs(pos_strand_ch_wise_data.iloc[i+1]['read_start_time'] - pos_strand_ch_wise_data.iloc[i]['read_end_time']) < gap_threshold) \
                    and (int(pos_strand_ch_wise_data.iloc[i+1]['alignment_start_position']) > int(pos_strand_ch_wise_data.iloc[i]['alignment_start_position'])+pos_strand_ch_wise_data.iloc[i]['alignment_length']):
                channel_wise_data_df_pos_strand = channel_wise_data_df_pos_strand.append(pos_strand_ch_wise_data.iloc[i], ignore_index=False)
                channel_wise_data_df_pos_strand = channel_wise_data_df_pos_strand.append(pos_strand_ch_wise_data.iloc[i+1], ignore_index=False)
        # Remove duplicates
        channel_wise_data_df_pos_strand = channel_wise_data_df_pos_strand.drop_duplicates(keep='first')
        # Sort time wise
        channel_wise_data_df_pos_strand = channel_wise_data_df_pos_strand.sort_values('read_start_time')
        channel_wise_dict_pos_strand[str(ch)] = channel_wise_data_df_pos_strand


    for ch in range(1, 512+1):
        print('processing channel {0} of 512 channels (- strands)'.format(ch))
        channel_wise_data_df_neg_strand = pd.DataFrame(columns=['read_start_time', 'read_end_time', 'sampling_rate',
                                                                'read_id', 'read_number', 'channel_id',
                                                                'fast5_sequence_length', 'class', 'strand',
                                                                'alignment_start_position', 'mapping_quality',
                                                                'sam_data_class', 'alignment_length'])

        neg_strand_ch_wise_data = neg_strand_data[neg_strand_data.channel_id == ch].sort_values('read_start_time')

        neg_strand_ch_wise_data_numrows = neg_strand_ch_wise_data.shape[0]

        for i in range(neg_strand_ch_wise_data_numrows-1):
            # check for two things
            # 1. is alignment of the two fragments back-to-back in time
            # 2. is the start position of the end fragment after the end position of the first fragment in terms of alignment

            if (abs(neg_strand_ch_wise_data.iloc[i+1]['read_start_time'] - neg_strand_ch_wise_data.iloc[i]['read_end_time']) < gap_threshold) \
                    and (int(neg_strand_ch_wise_data.iloc[i+1]['alignment_start_position'])+neg_strand_ch_wise_data.iloc[i+1]['alignment_length'] < int(neg_strand_ch_wise_data.iloc[i]['alignment_start_position'])):
                channel_wise_data_df_neg_strand = channel_wise_data_df_neg_strand.append(neg_strand_ch_wise_data.iloc[i], ignore_index=False)
                channel_wise_data_df_neg_strand = channel_wise_data_df_neg_strand.append(neg_strand_ch_wise_data.iloc[i+1], ignore_index=False)
        # Remove duplicates
        channel_wise_data_df_neg_strand = channel_wise_data_df_neg_strand.drop_duplicates(keep='first')
        # Sort time wise
        channel_wise_data_df_neg_strand = channel_wise_data_df_neg_strand.sort_values('read_start_time')
        channel_wise_dict_neg_strand[str(ch)] = channel_wise_data_df_neg_strand

    # ------------------------------------------------------ #
    # ---------------- RESULTS/STATISTICS ------------------ #
    # ------------------------------------------------------ #

    # Fragmented reads count
    fragment_reads_count = 0
    fragments_df = pd.DataFrame(columns=['read_start_time', 'read_end_time', 'sampling_rate',
                                         'read_id', 'read_number', 'channel_id',
                                         'fast5_sequence_length', 'class', 'strand',
                                         'alignment_start_position', 'mapping_quality',
                                         'sam_data_class', 'alignment_length'])
    for ch in range(1, 512+1):
        pos_ch_df = channel_wise_dict_pos_strand[str(ch)]
        neg_ch_df = channel_wise_dict_neg_strand[str(ch)]

        if not(pos_ch_df.empty):
            for i in range(pos_ch_df.shape[0]-1):
                fragments_df = fragments_df.append(pos_ch_df.iloc[i], ignore_index=False)
                if (pos_ch_df.iloc[i+1]['read_start_time'] - pos_ch_df.iloc[i]['read_end_time']) > gap_threshold:
                    fragment_reads_count += 1
            fragment_reads_count += 1
            fragments_df = fragments_df.append(pos_ch_df.iloc[i+1], ignore_index=False)

        if not(neg_ch_df.empty):
            for i in range(neg_ch_df.shape[0]-1):
                fragments_df = fragments_df.append(neg_ch_df.iloc[i], ignore_index=False)
                if (neg_ch_df.iloc[i+1]['read_start_time'] - neg_ch_df.iloc[i]['read_end_time']) > gap_threshold:
                    fragment_reads_count += 1
            fragment_reads_count += 1
            fragments_df = fragments_df.append(neg_ch_df.iloc[i+1], ignore_index=False)

    print('Total number of unclassifed reads: {0}'.format(uncl_df.shape[0]))
    print('Total number of barcode10(DMS treated and made abasic) reads: {0}'.format(bc10_df.shape[0]))
    print('Combined total of unclassified and barcode10 reads: {0}'.format(uncl_df.shape[0]+bc10_df.shape[0]))
    print('Total number of longer reads found by combining fragments: {0}'.format(fragment_reads_count))
    print('Number of read fragments coming from unclassified: {0}'.format(fragments_df.sam_data_class.value_counts()['unclassified']))
    print('Number of read fragments coming from barcode10: {0}'.format(fragments_df.sam_data_class.value_counts()['barcode10']))





