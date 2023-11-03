"""
Created on Fri Nov  3 14:40:36 2023

@author: Christopher
"""

import numpy as np
from astropy.io import fits
from tkinter import Tk
from tkinter import filedialog
import os


def get_data():

    Tk().withdraw()
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        filename = (os.path.basename(file_path).replace(".fits", ""))
        main(data, filename, file_path)
        data.close()
        
    return


def main(data, filename, path):
    
    num_channels = data[0].shape[0]
    temp = []
    
    for i in range(num_channels):
        max_x, max_y = find_max_value_coordinates(data[0].data[i])
        temp.append([max_x, max_y])
    
    temp = np.vstack((temp))
    x_mean = np.mean(temp[:,0])
    y_mean = np.mean(temp[:,1])
    
    print('shape = ' + str(data[0].shape))
    print('file = ' + filename)
    print('centre = ' + str(round(x_mean)) + ',' + str(round(y_mean)))
    
    return 


def find_max_value_coordinates(arr):
    max_value = np.nanmax(arr)
    max_x, max_y = np.where(arr == max_value)
    return max_x[0], max_y[0]

get_data()
