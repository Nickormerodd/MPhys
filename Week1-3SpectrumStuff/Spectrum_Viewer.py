# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 13:58:39 2023

@author: Christopher
"""

import numpy as np
from astropy.io import fits
from spectral_cube import SpectralCube as sc
from spectral_cube import DaskSpectralCube
from tkinter import Tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
import radio_beam as beam
from scipy.constants import c
from astropy import units as u
from scipy.optimize import curve_fit
from matplotlib.ticker import FuncFormatter
import uncertainties.unumpy as unp
from uncertainties import ufloat

def get_data():

    Tk().withdraw()  
    file_paths = filedialog.askopenfilenames()
    n = 0
    temp = []
    for file_path in file_paths:
        velocity_name = 0
        vel = []
        n = n + 1
        data = fits.open(file_path)
        filename = (os.path.basename(file_path).replace(".fits", "")).replace("_"," ").replace("."," ")
        #velocity_name = main(data, filename, file_path, vel, 0, 0)
        #temp.append(velocity_name)
        main(data,filename)
        data.close()
        #print(file_path)
        
    return

def main(data, filename):
    print('Shape = ' + str(data[0].shape))
    num_channels = data[0].shape[0]
    header = data[0].header
    print(data[0].shape)
    ref_freq = header['CRVAL3'] #reference frequency
    incr_freq = header['CDELT3'] #frequency increment
    freq = ((ref_freq + incr_freq * (np.arange(num_channels))) / 10**9 )
    
    temp = []    
    for i in range(num_channels):
        max_value = np.nanmean(data[0].data[i]), freq[i]
        temp.append(max_value)
    
    
    channel_no = np.linspace(0,num_channels-1,num_channels)
    spectral = np.hstack((np.vstack((channel_no)),np.vstack((temp))))
    
    print('Mean = ' + str(np.mean(np.sqrt(spectral[:,1]**2))))
    #print(spectral)
    #spectrum(spectral, filename)
    """
    val_str = input('Range of channel numbers you want, e.g. in format 100-120:\n')
    val_list = val_str.split('-')

    if len(val_list) == 2:
        val_1 = float(val_list[0])
        val_2 = float(val_list[1])
    else:
        print("Invalid input format. Please use the format 'start-end'.")
        
    print('Range = {}-{}'.format(val_1,val_2))
    
    #data = filter_data(spectral,val_1,val_2)
    """
    spectrum(spectral, filename)
    cube = sc.read(data)
    
    spectrum(spectral, filename)
    cube = DaskSpectralCube.read(data)
    stats = cube.statistics()
    rms = stats['rms'].value
    epic = float(rms)
    print(type(epic))
    print(epic)
    #print(type(rms))
    print(f'RMS: {rms}')
    print(f'RMS: {rms}', rms**2)
    
    epic_rms = np.mean(np.sqrt(spectral[:,1]**2))
    
    cube_data = cube.filled_data[:]
    epicals = np.hstack((np.vstack((channel_no)),np.vstack((temp))))
    epicals[:,1] = np.sqrt(epicals[:,1]**2)
    
    print(epicals)
    
    spectrum_2(spectral, epicals, filename, epic_rms)
    #print(type(rms))
    return


def spectrum(data, filename):
    #fig = plt.figure()
    #ax = fig.add_subplot()
    fig, ax = plt.subplots(figsize=(12, 6)) 
    
    ax.plot(data[:,2], data[:,1], linewidth=1)
    ax.set_title(filename)
    ax.set_xlabel('Frequency, GHz')
    ax.set_ylabel('Brightness temperature, K')
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    
    ax2 = ax.twiny()
    num_ticks_to_display = 15
    step_size = len(data[:, 0]) // num_ticks_to_display
    channel_numbers = data[::step_size, 0][::-1] 
    ax2.set_xticks(data[::step_size, 0][::-1])
    ax2.set_xticklabels(channel_numbers.astype(int)[::-1])
    ax2.tick_params(axis='x', labelsize=9)
    ax2.set_xlabel('Channel number')
    ax2.spines['top'].set_color('none')
    
    ax.grid(True, linewidth=0.5, alpha=0.7)
    filename = (filename.replace(" ", "_")).replace("image_pbcor_line_galactic_kelvin", "spectrum.png")
    print(filename)
    #plt.savefig(filename,dpi=1000)

    return

def spectrum_2(data, data_2, filename, rms):
    #fig = plt.figure()
    #ax = fig.add_subplot()
    fig, ax = plt.subplots(figsize=(12, 6)) 
    
    ax.plot(data[:,2], data[:,1], linewidth=1, color='blue', label='Data')
    ax.plot(data_2[:,2], data_2[:,1], linewidth=1, color='red', label = 'RMS Data')
    ax.plot(data_2[:,2],np.linspace(rms,rms,len(data_2[:,1])),label='RMS Noise = {:.3f} K'.format(rms),
            color='purple', linestyle='--', alpha=0.7)
    ax.set_title(filename)
    ax.set_xlabel('Frequency, GHz')
    ax.set_ylabel('Brightness temperature, K')
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    """
    ax2 = ax.twiny()
    num_ticks_to_display = 15
    step_size = len(data[:, 0]) // num_ticks_to_display
    channel_numbers = data[::step_size, 0][::-1] 
    ax2.set_xticks(data[::step_size, 0][::-1])
    ax2.set_xticklabels(channel_numbers.astype(int)[::-1])
    ax2.tick_params(axis='x', labelsize=9)
    ax2.set_xlabel('Channel number')
    ax2.spines['top'].set_color('none')
    """
    ax.grid(True, linewidth=0.5, alpha=0.7)
    filename = (filename.replace(" ", "_")).replace("image_pbcor_line_galactic_kelvin", "spectrum.png")
    print(filename)
    ax.legend()
    filename = filename + '_spectrum.png'
    plt.savefig(filename,dpi=1000)

    return

def filter_data(data, val_1, val_2):

    first_column = data[:, 1]
    mask = (first_column >= val_1) & (first_column <= val_2)
    filtered_data = data[mask]
    
    return filtered_data

get_data()
