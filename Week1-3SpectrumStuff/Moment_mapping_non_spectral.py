"""
Created on Sat Oct  7 19:19:14 2023

@author: Christopher

This is my initial attempt at producing a moment map using spectralcube.

This code is named 'non spectral' because I haven't read it into
spectralcube immediately, only when I want the moment plot, i.e all the data
manipulation wasnt done in spectralcube.

This code seems to work, allowing you to input channel range, but can only
produce a 0th moment plot without error.

Also, I believe the spectral axis (frequency) needs to be changed into
velocity if we want a plot that shows anything.

You have to choose files which have already been converted to kelvin.

"""

import numpy as np
from astropy.io import fits
from spectral_cube import SpectralCube as sc
from tkinter import Tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
import radio_beam as beam
from scipy.constants import c
from astropy import units as u
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties import ufloat
from astropy.wcs import WCS

def get_data():

    Tk().withdraw()  
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        filename = (os.path.basename(file_path).replace(".fits", "")).replace("_"," ").replace("."," ")
        main(data, filename, file_path)
        data.close()

    return

def main(data, filename, path):
    
    header = data[0].header
    ref_freq = header['CRVAL3'] #reference frequency
    incr_freq = header['CDELT3'] #frequency increment

    num_channels = data[0].shape[0]  
    frequencies = ((ref_freq + incr_freq * (np.arange(num_channels))) / 10**9 ) *u.GHz
    freq = ((ref_freq + incr_freq * (np.arange(num_channels))) / 10**9 ) 
        
    val_str = input('Range of channel numbers, e.g., in format 100-120:\n')
    val_list = val_str.split('-')

    if len(val_list) == 2:
        val_1 = int(val_list[0])
        val_2 = int(val_list[1])
    else:
        print("Invalid input format. Please use the format 'start-end'.")  
    
    
    for i in range(data[0].shape[0]): # maybe not necessary, here cos tryna do 1st moment without error
        data[0].data[i] = np.nan_to_num(data[0].data[i], nan=0.0)
        data[0].data[i] = np.where(np.isinf(data[0].data[i]), 0.0, data[0].data[i])
    
    print(data[0].data[0])
    filtered_data = filter_data(data, val_1, val_2)    
    print('Shape = ' + str(filtered_data[0].shape)) 
    
    #wcs = WCS(data[0].header)
    cube_data = sc.read(filtered_data)
    moments(cube_data)
    
    return
    
def filter_data(data, val_1, val_2):

    if val_1 >= 0 and val_2 >= 0 and val_1 <= val_2 and val_2 < len(data[0].data):
        data[0].data = data[0].data[val_1:val_2 + 1]
    
    return data

"""
def moments(cube_data):
     # Check for NaN and zero values in the cube_data
    mask = np.isnan(cube_data.filled_data[:])  # Extract the data and fill NaNs
    mask = np.logical_or(mask, np.isinf(cube_data.filled_data[:]))
    mask = np.logical_or(mask, cube_data.filled_data[:] == 0)

    # Handle problematic values by masking or replacing them
    cube_data = np.nan_to_num(cube_data.filled_data[:], nan=0)
    cube_data = np.ma.masked_array(cube_data, mask=mask)

    # Calculate the moment
    moment_1 = cube_data.moment(order=1)
    
    return moment_1
"""

def moments(cube_data):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    moment0 = cube_data.moment(order=0, axis=0)
    print(moment0)
    
    m = moment0.shape[0]
    p = moment0.shape[1]

    x = np.linspace(0, p - 1, p)
    y = np.linspace(0, m - 1, m)
    X, Y = np.meshgrid(x, y)
    
    contourf = ax.contourf(X, Y, moment0 , cmap='magma',
                           levels = 30) # 'viridis')
    plt.colorbar(contourf)#, label='Brightness Temperature, K')
    
    return moment0

get_data()
    
