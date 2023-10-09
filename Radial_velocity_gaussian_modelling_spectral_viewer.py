"""
Created on Fri Oct  6 19:07:00 2023

@author: Christopher Andrews

This code lets you input n .fits files (it doesn't matter if these are
in kelvin or not), choose channel number ranges of interest, input 
molecule name, and input rest frequency of the molecule.

This code then models the molecular line peak as a gaussian, and then extracts
the mean value of the observed frequency of that peak and calculates
the velocity the source is moving at by comparing with the rest frequency
inputted.

This is repeated m times for n files, however many you choose.

It then prints a final velocity which is calculated from the mean of the 
velocities, with an uncertainty on the mean velocity equal to the 
standard deviation of the calculated velocities.
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
        velocity_name = main(data, filename, file_path, vel, 0, 0)
        temp.append(velocity_name)
        data.close()
        
    final_velocity(np.vstack(temp))
    
    return

def main(data, filename, path, vel, spectral, n):
    print('\nFile = ' + str(filename))
    
    header = data[0].header
    ref_freq = header['CRVAL3'] #reference frequency
    incr_freq = header['CDELT3'] #frequency increment
    
    beam_major = np.deg2rad(header['BMAJ'])  # Major axis of the beam 
    beam_minor = np.deg2rad(header['BMIN'])  # Minor axis of the beam
    beam_size = (np.pi/(4 *np.log(2))) * beam_major * beam_minor
    
    num_channels = data[0].shape[0]  
    frequencies = ((ref_freq + incr_freq * (np.arange(num_channels))) / 10**9 ) *u.GHz
    freq = ((ref_freq + incr_freq * (np.arange(num_channels))) / 10**9 ) 
    beam_size = beam_size* u.sr

    temp = []
    
    for i in range(num_channels):
        if not 'kelvin' in filename.lower() and n==0:
            equiv = u.brightness_temperature(frequencies[i])
            data[0].data[i] = (data[0].data[i]*u.Jy/beam_size).to(u.K,equivalencies=equiv)

        max_value = i, np.nanmax(data[0].data[i]), freq[i]
        temp.append(max_value)
    
    if n==0:
        spectral = np.vstack(temp)
        spectrum(spectral, filename)
    
    if not n==0:
        again = input('Anymore ranges in this file? [y/n]\n')
    
    if n == 0 or again == 'y':
        velocity, name = answers(spectral)
        epic = np.array([velocity, name], dtype=object)
        vel.append(epic)
        main(data,filename,path,vel,spectral,1)
    return vel

def answers(spectral):
    
    val_str = input('Range of channel numbers, e.g., in format 100-120:\n')
    val_list = val_str.split('-')

    if len(val_list) == 2:
        val_1 = int(val_list[0])
        val_2 = int(val_list[1])
    else:
        print("Invalid input format. Please use the format 'start-end'.")
        
    print('Range = {}-{}'.format(val_1,val_2))
    molecule_name = input('Name of Molecule:\n')
    print('Name of molecule = {}'.format(molecule_name))
    line_data = filter_data(spectral, val_1, val_2)
    line_data_gaus, params, covariance = gaus_fitting(line_data)
    #perr = np.sqrt(np.diag(covariance))
    line_graph(line_data, line_data_gaus, params, molecule_name, covariance)
    decision = input('Velocity calculation? [y/n]\n')
    if decision=='y':
        velocity = velocity_func(params[1], molecule_name)
    else:
        velocity = np.nan
    
    return velocity, molecule_name

def spectrum(data, filename):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(data[:,2], data[:,1], linewidth = 1)
    ax.set_title(filename)
    ax.set_xlabel('Frequency, GHz')
    ax.set_ylabel('Brightness temperature, K')
    
    ax2 = ax.twiny()
    data[:, 0] = data[::-1, 0]
    ax2.scatter(data[:,0], data[:,1], s=0, color='green')
    ax2.tick_params(axis='x')
    ax2.set_xlabel('Channel number')
    ax.grid(True, linewidth=0.5, alpha=0.7)

    return

def line_graph(data, gaus, params, name, covariance):
    """
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.errorbar(data[:,2], data[:,1], linewidth = 1, 
                label='{} line'.format(name), fmt='x')
    ax.set_title(name)
    ax.set_xlabel('Frequency, GHz')
    ax.set_ylabel('Brightness temperature, K')
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    
    ax2 = ax.twiny()
    data[:, 0] = data[::-1, 0]
    ax2.scatter(data[:,0], data[:,1], s=0, color='green')
    ax2.tick_params(axis='x', labelsize=9)
    ax2.set_xlabel('Channel number')
    
    amp, mean, std_dev = params
    perr = np.sqrt(np.diag(covariance))
    
    y_mean = np.linspace(min(data[:,1]), max(data[:,1])-np.mean(data[:,1])*0.2, 2)
    x_mean = np.linspace(mean,mean,2)
    y_std_1 = np.linspace(0.5*(max(data[:,1])-min(data[:,1])),
                          0.5*(max(data[:,1])-min(data[:,1])), 2)
    x_std_1 = np.linspace(mean,mean+std_dev,2)
    y_std_2 = np.linspace(0.5*(max(data[:,1])-min(data[:,1])),
                          min(data[:,1]), 2)
    x_std_2 = np.linspace(mean+std_dev,mean+std_dev,2)
    
    mean_formatted = '{:.4f}'.format(mean)
    std_formatted = '{:.4f}'.format(std_dev)
    ax.plot(data[:,2], gaus, color='green')
    ax.plot(x_mean,y_mean,color='grey', linewidth=1, linestyle='--',
            label =r'$\bar{\nu}$ = ' + str(mean_formatted) + ' GHz')

    arrow_props = dict(arrowstyle='->', linewidth=1, color='black')
    ax.annotate('', xy=(x_std_1[0], y_std_1[0]), xytext=(x_std_1[1], y_std_1[1]),
                 arrowprops=arrow_props)
    ax.annotate('', xy=(x_std_1[-1], y_std_1[-1]), xytext=(x_std_1[-2], y_std_1[-2]),
                 arrowprops=arrow_props)
    ax.plot([],[], color='black', label= r'$\sigma$ = ' + str(std_formatted) + ' GHz',
            linewidth=1)
    ax.plot(x_std_2,y_std_2,color='blue',linestyle='--',linewidth=1)
    
    def format_func(x, _):
        return '{:.2f}'.format(x)
    ax.xaxis.set_major_formatter(FuncFormatter(format_func))
    
    print('\n''Mean_freq = ' + str(mean) + r' +/- ' + str(perr[1]) + ' GHz')
    print('Std_dev = ' + str(std_dev) + r' +/- ' + str(perr[2]) + ' GHz')
    print('Amplitude = ' + str(amp) + r' +/- ' + str(perr[0]) + ' K')
    ax.grid(True, linewidth=0.5, alpha=0.7)
    ax.legend(fontsize=9.5)
    return
    
def gaussian(x, amplitude, mean, std_dev):
    
    return amplitude * np.exp(-(x - mean) ** 2 / (2 * std_dev ** 2))
    
def gaus_fitting(data):

    x_bar = np.mean(data[:,2])
    x_std = np.std(data[:,2])
    amp = max(data[:,1])-min(data[:,1])
    minimum = min(data[:,1])
    data[:,1] = data[:,1] - minimum
    
    GUESS = amp, x_bar, x_std
    n=0
    while n < 5:
        try:
            params, covariance = curve_fit(gaussian, data[:,2], data[:,1],
                                           p0=GUESS, absolute_sigma=True,
                                           maxfev=1000)
        except ValueError or RuntimeWarning:
            break
        GUESS = params
        n+=1

    amp, mean, std_dev = params
    data[:,1] = data[:,1] + minimum
    data_curve = gaussian(data[:,2],amp,mean,std_dev) + minimum
    params = (amp+minimum,mean,std_dev)
    return data_curve, params, covariance

def filter_data(data, val_1, val_2):

    first_column = data[:, 0]
    mask = (first_column >= val_1) & (first_column <= val_2)
    filtered_data = data[mask]
    
    return filtered_data

def velocity_func(f_obs, name):
    
    f_rest = float(input('What is the rest frequency in GHz for '+str(name)+'\n'))
    velocity = -c*(f_rest-f_obs)/f_rest
    print('v_radial_'+ str(name) + ' = ' + str(velocity*10**(-3)) + ' km/s')
    
    return velocity

def final_velocity(data):
    print('\n')
    data = no_nan(data)
    data[:,0] = data[:,0]*10**(-3)
    if len(data[:,0]) > 1:
        for line in data:
            print('velocity_'+str(line[1])+' = ' + str(line[0]) + ' km/s')
            
        mean_v = np.mean(data[:,0])
        std_dev_v = np.std(data[:,0])
        print('\nMean_velocity = ' +  str(mean_v) + ' +/- ' + str(std_dev_v) + ' km/s')
    return

def no_nan(data):
    temp = []
    for line in data:
        if not np.any(np.isnan(line[0])):
            temp.append(line)
    return np.vstack(temp)

get_data()
