import sys
import glob
import numpy as np
from datetime import datetime
import csv
import pylab


def write_list_to_csv(target_list, path_to_store_file):
    '''
    Write 1D list to csv as single column
    '''

    temp = []
    for tl in target_list:
        temp.append([tl])  # this step creates a list of lists as such [[a],[b]], to pass to writerows function of csv

    with open(path_to_store_file, "w", newline="") as f:
        print(path_to_store_file)
        writer = csv.writer(f)
        writer.writerows(temp)

def read_csv_to_list(file_path):

    '''
    Read single column csv file as a list
    '''

    var_list = []

    with open(file_path,"r") as file:
        reader = csv.reader(file)

        for row in reader:
            #row = [float(r) for r in row]
            row = float(row[0])
            var_list.append(row)

    return var_list

def movingaverage(values, window_size, mode):
    '''
    Get moving average of list of values with specified window size
    '''

    import numpy as np

    window = np.ones(int(window_size)) / float(window_size)

    return np.convolve(values, window, mode)  # 'valid' 'full' 'same'

def flatten(list_of_lists):

    '''
    Creates a flat list out of a list of lists
    '''

    flat_list = [item for sublist in list_of_lists for item in sublist]

    return flat_list

def get_line_colors(cmap,num_colors):

    '''
    Get list of colors from a matplotlib colormap
    '''

    cm = pylab.get_cmap(cmap)

    color_list = []
    for i in range(num_colors):
        color = cm(1.*i/num_colors)  # color will now be an RGBA tuple
        color_list.append(color)

    return color_list