"""
Created on Thu Nov 16 10:57:28 2023

@author: Christopher
"""

import numpy as np
from astropy.io import fits
from tkinter import Tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from mpl_toolkits.mplot3d import Axes3D as mpl
from matplotlib.lines import Line2D
from scipy.constants import c

########## params ################
#calibration_freq = 220.70908  # GHz
#transition_freq = 220.593987 # GHz
#obs_freq = 220.5725      # GHz
output_filename = 'C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/new_data/CH3CN_K/'
####################################

def get_data():

    Tk().withdraw()
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        filename = (os.path.basename(file_path)).replace(".fits", "_cal.fits")#.replace(".image.pbcor_line.galactic","")
        #main(data, filename, file_path)
        change_header_info(data, output_filename, file_path, filename)
        data.close()
        
    return



def change_header_info(hdul, output_fits_file, path, name):
    """
    Update the header in the output FITS file based on the input FITS file.
    """
    
    print(name)
    header = hdul[0].header
    #print(header)
    if 'k=0' in name.lower():
        transition_freq = 220.747369
        obs_freq = calculations(hdul)
        VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
        header['RESTFRQ'] = transition_freq * 10**9 
        header['VELREF'] = VELREF
        k_number = 0


    elif 'k=1' in name.lower():
        transition_freq = 220.743012
        obs_freq = calculations(hdul)
        VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
        header['RESTFRQ'] = transition_freq * 10**9 
        header['VELREF'] = VELREF
        k_number = 1


    elif 'k=2' in name.lower():
        transition_freq = 220.73027
        obs_freq = calculations(hdul)
        VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
        header['RESTFRQ'] = transition_freq * 10**9 
        header['VELREF'] = VELREF
        k_number = 2


    elif 'k=3' in name.lower():
        transition_freq = header['RESTFRQ']
        obs_freq = header['RESTFRQ']
        VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
        header['RESTFRQ'] = transition_freq * 10**9 
        VELREF = header['VELREF']
        k_number = 3

    elif 'k=4' in name.lower():
        transition_freq = 220.679114
        obs_freq = calculations(hdul)
        VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
        header['RESTFRQ'] = transition_freq * 10**9 
        header['VELREF'] = VELREF
        k_number = 4


    elif 'k=5' in name.lower():
        transition_freq = 220.64112
        obs_freq = calculations(hdul)
        VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
        header['RESTFRQ'] = transition_freq * 10**9 
        header['VELREF'] = VELREF
        k_number = 5


    elif 'k=6' in name.lower():
        transition_freq = 220.593987
        obs_freq = calculations(hdul)
        VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
        header['RESTFRQ'] = transition_freq * 10**9 
        header['VELREF'] = VELREF
        k_number = 6


    elif 'k=7' in name.lower():
        transition_freq = 220.538635
        obs_freq = calculations(hdul)
        k_number = 7
        VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
        header['RESTFRQ'] = transition_freq * 10**9 
        header['VELREF'] = VELREF


    elif 'k=8' in name.lower():
        transition_freq = 220.474767
        obs_freq = calculations(hdul)
        VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
        header['RESTFRQ'] = transition_freq * 10**9 
        header['VELREF'] = VELREF
        k_number = 8

    else:
        title = 'unknown'
        k_number = 'unknown'
        
    #print(obs_freq)
    #print(transition_freq)
    #print(VELREF)
    
    VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
    print("new header will be", VELREF)
    # Read data and header from the input file
    data = hdul[0].data
    #print("old header = ", header['VELREF'],  "\ninput file =", str(path))
    output_fits_file = output_fits_file + name
    print("output file = file =", str(output_fits_file))
    # Write the original data with the updated header to the output file
    fits.writeto(output_fits_file, data, header, overwrite=True)
    
def calculations(data):
    
    header = data[0].header
    ref_freq = header['CRVAL3'] #reference frequency
    incr_freq = header['CDELT3'] #frequency increment
    
    num_channels = data[0].shape[0]  
    #frequencies = ((ref_freq + incr_freq * (np.arange(num_channels))) / 10**9 ) *u.GHz
    freq = ((ref_freq + incr_freq * (np.arange(num_channels))) / 10**9 ) 
    #beam_size = beam_size* u.sr

    temp = []
    
    for i in range(num_channels):
        
        max_value = np.nanmean(data[0].data[i]), freq[i]
        temp.append(max_value)
    
    temp = np.vstack((temp))
    
    mean = gaus_fitting(temp)
    
    return mean
    
def gaussian(x, amplitude, mean, std_dev):
    
    return amplitude * np.exp(-(x - mean) ** 2 / (2 * std_dev ** 2))
    
def gaus_fitting(data):

    x_bar = np.mean(data[:,1])
    x_std = np.std(data[:,1])
    amp = max(data[:,0])-min(data[:,0])
    minimum = min(data[:,0])
    data[:,1] = data[:,1] - minimum
    
    GUESS = amp, x_bar, x_std
    n=0
    while n < 5:
        try:
            params, covariance = curve_fit(gaussian, data[:,1], data[:,0],
                                           p0=GUESS, absolute_sigma=True,
                                           maxfev=3000)
        except ValueError or RuntimeWarning:
            break
        GUESS = params
        n+=1

    amp, mean, std_dev = params

    return mean

get_data()
