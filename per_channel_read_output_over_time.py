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
import helpers as hlp
import collect_read_lists as crl

from bokeh.io import output_file, show
from bokeh.plotting import figure, output_file
from bokeh.layouts import layout
from bokeh.models import BoxAnnotation
import helpers as hlp
from bokeh.models import NumeralTickFormatter

LOAD_FILEPATHS_FROM_TEXT_FILE = False
DATASET_PATH = "/Users/adnaniazi/mnt/kjempetuja/export/valenfs/data/processed_data/" \
               "MinION/20180710_1431_197_lig/backbone_aligned_reads/pass/barcode10"

file_paths = crl.make_read_lists(folder_path=DATASET_PATH,
                                 load_filepaths_from_text_file=LOAD_FILEPATHS_FROM_TEXT_FILE)

channel_start_end_dict, channel_read_id_dict = hlp.make_channel_start_end_dict(file_paths)

y_list = []
x_list = []
channels_list = []
for k,v in channel_start_end_dict.items():
    p = figure(plot_width=1400, plot_height=720)
    channel_number = int(k.replace('Ch-', ''))

    for v_list in v:
        y_list.append([channel_number, channel_number])
        x_list.append(v_list)
    channels_list.append(k)

# filter reads according to the distance


p.multi_line(xs=x_list, ys=y_list, line_color="red")
p.yaxis[0].formatter = NumeralTickFormatter(format="0")
p.xaxis.axis_label = "Start and end of reads  (time in seconds after experiment start)"
p.yaxis.axis_label = "Channel Number (0 - 512)"
output_file('per_channel_read_output_barcode10.html')
show(p)


if __name__ == '__main__':
    pass
