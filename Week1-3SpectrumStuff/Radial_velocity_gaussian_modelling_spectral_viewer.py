"""
Created on Fri Oct  6 19:07:00 2023

@author: Christopher Andrews

This code lets you input n .fits files (it doesn't matter if these are
in kelvin or not), choose channel number ranges of interest, input 
molecule name, and input rest frequency of the molecule.

This code then models the molecular line peak as a gaussian, and then extracts
the mean value of the observed frequency of that peak and calculates
the velocity the source is moving at by comparing with the rest frequency
inputted. The code also evaluates a realistic uncertainty on the brightness 
temperature, which is calculated by finding the uncertainty required for the
reduced chi-square to be within 0.5-1.5.

This is repeated m times for n files, however many you choose.

It then prints a final velocity which is calculated from the mean of the 
velocities, with an uncertainty on the mean velocity equal to the 
standard deviation of the calculated velocities.

A histogram of the velocities are then plotted on a graph, with a gaussian 
plotted through them to show that the spread of the calculated velocities.

This code saves the figures created, so # that off if dont want figure saved.
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
from scipy.stats import norm

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
    print('Number of channels = ' + str(data[0].shape[0]))
    
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

        max_value = np.nanmean(data[0].data[i]), freq[i]
        temp.append(max_value)
    
    if n==0:
        channel_no = np.linspace(0,num_channels-1,num_channels)
        spectral = np.hstack((np.vstack((channel_no)),np.vstack((temp))))
        #print(spectral)
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
    
    val_str = input('Range of channel numbers, e.g. in format 100-120:\n')
    val_list = val_str.split('-')

    if len(val_list) == 2:
        val_1 = int(val_list[0])
        val_2 = int(val_list[1])
    else:
        print("Invalid input format. Please use the format 'start-end'.")
        
    print('Range = {}-{}'.format(val_1,val_2))
    molecule_name = input('Name of Molecule (in latex):\n')
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

    fig, ax = plt.subplots(figsize=(15, 6)) 
    
    ax.plot(data[:,2], data[:,1], linewidth=1)
    ax.set_title(filename)
    ax.set_xlabel('Frequency, GHz')
    ax.set_ylabel('Brightness temperature, K')
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    
    ax2 = ax.twiny()
    num_ticks_to_display = 6 
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
    #plt.savefig('')

    return

def line_graph(data, gaus, params, name, covariance):
    """
    data[:,0] = channel number
    data[:,1] = Brightness temperature
    data[:,2] = frequency
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    low_uncert, up_uncert, mean_uncert = chi_uncertainty(data, gaus)
    
    ax.errorbar(data[:,2], data[:,1], yerr=float(mean_uncert), linewidth = 1, 
                label='Data', fmt='x')
    ax.set_title(name)
    ax.set_xlabel('Frequency, GHz')
    ax.set_ylabel('Brightness temperature, K')
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)

    
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
    x_gaus = np.linspace(min(data[:,2]), max(data[:,2]), 100)
    y_gaus = gaussian(x_gaus,amp,mean,std_dev)
    ax.plot(x_gaus, y_gaus, color='green')
    ax.plot(x_mean,y_mean,color='grey', linewidth=1, linestyle='--',
            label =r'$\bar{\nu}$ = ' + str(mean_formatted) + ' GHz')

    arrow_props = dict(arrowstyle='->', linewidth=1, color='black')
    ax.annotate('', xy=(x_std_1[0], y_std_1[0]), xytext=(x_std_1[1], y_std_1[1]),
                 arrowprops=arrow_props)
    ax.annotate('', xy=(x_std_1[-1], y_std_1[-1]), xytext=(x_std_1[-2], y_std_1[-2]),
                 arrowprops=arrow_props)
    ax.plot([],[], color='black', label= r'$\sigma_{std}$ = ' + str(std_formatted) + ' GHz',
            linewidth=1)
    ax.plot(x_std_2,y_std_2,color='grey',linestyle='--',linewidth=1)
    
    ax.annotate(r'For 0.5$\leq$$\chi^2_{red}$$\leq$1.5 $\Longrightarrow$ ' + str(low_uncert) + '$\leq$$\sigma_{data}$$\leq$' + str(up_uncert),
                xy=(0.5, 0.5), xytext=(0, -145), 
                xycoords='axes fraction', textcoords='offset points', fontsize=11,
                ha='center', va='top')
    
    ax.annotate(r'$\therefore$ $\sigma_{data}$ = ' + str(mean_uncert) + ' K',
                xy=(0.5, 0.5), xytext=(0, -162), 
                xycoords='axes fraction', textcoords='offset points', fontsize=11,
                ha='center', va='top')
    
    def format_func(x, _):
        return '{:.3f}'.format(x)
    ax.xaxis.set_major_formatter(FuncFormatter(format_func))
    
    print('\n''Mean_freq = ' + str(mean) + r' +/- ' + str(perr[1]) + ' GHz')
    print('Std_dev = ' + str(std_dev) + r' +/- ' + str(perr[2]) + ' GHz')
    print('Amplitude = ' + str(amp) + r' +/- ' + str(perr[0]) + ' K')
    ax.grid(True, linewidth=0.5, alpha=0.7)
    ax.legend(fontsize=9.5)
    
    name = ((name.replace('$','')).replace(' ','')).replace("\\rightarrow","")
    
    plt.savefig('Gaussian_CH3CN_'+str(name) + '.png', dpi=1000, bbox_inches='tight')
    
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
                                           maxfev=2000)
        except ValueError or RuntimeWarning:
            break
        GUESS = params
        n+=1

    amp, mean, std_dev = params
    data[:,1] = data[:,1] + minimum
    data_curve = gaussian(data[:,2],amp,mean,std_dev) + minimum
    #data_curve = np.hstack((np.vstack((data_curve)),data[:,2]))
    params = (amp+minimum,mean,std_dev)
    return data_curve, params, covariance

def filter_data(data, val_1, val_2):

    first_column = data[:, 0]
    mask = (first_column >= val_1) & (first_column <= val_2)
    filtered_data = data[mask]
    
    return filtered_data

def velocity_func(f_obs, name):
    
    f_rest = float(input('What is the rest frequency in GHz for '+str(name)+'\n'))
    velocity = c*(f_rest-f_obs)/f_rest
    print('v_radial_'+ str(name) + ' = ' + str(velocity*10**(-3)) + ' km/s')
    
    return velocity

def final_velocity(data):
    print('\n')
    data = no_nan(data)
    if len(data[:,0]) > 1:
        data[:,0] = data[:,0]*10**(-3)    
        
        for line in data:
            print('velocity_'+str(((line[1]))) +' = ' + str(line[0]) + ' km/s')
            
        mean_v = np.mean(data[:,0])
        std_dev_v = np.std(data[:,0])
        velocity_graph(data[:,0], std_dev_v)
        print('\nMean_velocity = ' +  str(mean_v) + ' +/- ' + str(std_dev_v) + ' km/s')
    return

def no_nan(data):
    temp = []
    for line in data:
        if not np.any(np.isnan(line[0])):
            temp.append(line)
    return np.vstack(temp)

def chi_uncertainty(data, gaus):
    
    upper_uncert = (sum((data[:,1]-gaus)**2) / ((len(data)-1)*0.5) )**(0.5)
    lower_uncert = (sum((data[:,1]-gaus)**2) / ((len(data)-1)*1.5))**(0.5)
    mean_uncert = (lower_uncert + upper_uncert)/2
    
    lower_uncert = '{:.2f}'.format(lower_uncert)
    upper_uncert = '{:.2f}'.format(upper_uncert)
    mean_uncert = '{:.2f}'.format(mean_uncert)
    
    print('\nUncertainty = ' + str(lower_uncert) + ' - ' + str(upper_uncert))
    print('Mean uncertainty = ' + str(mean_uncert) + ' K')
    
    return lower_uncert, upper_uncert, mean_uncert

def velocity_graph(velocities, std):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    range_width = std
    
    # Initialize an empty list to store data points
    data_points = []
    
    # Generate 100 data points for each velocity
    num_points = 100
    for velocity in velocities:
        data_points += list(np.linspace(velocity - range_width, velocity + range_width, num_points))
    
    # Create a histogram
    hist, bins, _ = plt.hist(data_points, bins=30, edgecolor='black', alpha=0.7,
                             color='indigo')
    bin_centers = (bins[1:] + bins[:-1]) / 2
    
    # Perform the curve fit with better initial guesses
    mean_estimate = np.mean(data_points)
    stddev_estimate = np.std(data_points)
    params, _ = curve_fit(gaussian, bin_centers, hist, p0=[max(hist), mean_estimate, stddev_estimate])
    
    # Create the fitted Gaussian curve
    fit_curve = gaussian(bin_centers, *params)
    
    # Plot the histogram and the fitted Gaussian curve
    ax.plot(bin_centers, fit_curve, label='Gaussian Fit', color ='gold')
    ax.plot([],[],label="v = {:.3f} km/s".format(params[1]), alpha =0)
    ax.plot([],[],label="$\sigma$ = {:.3f} km/s".format(params[2]), alpha=0)
    ax.set_title('Weighted Velocity Distribution')
    ax.set_xlabel('Velocity')
    ax.set_ylabel('Frequency')
    ax.grid(True, alpha=0.2)
    ax.legend()
    ax.annotate(r'Mean = ' + str('{:.3f}'.format(np.mean(velocities))) + ' $\pm$ ' + str('{:.3f}'.format(std)) + ' km/s',
            xy=(0.5, 0.5), xytext=(0, -145), 
            xycoords='axes fraction', textcoords='offset points', fontsize=11,
            ha='center', va='top')
    ax.annotate(r'Weighted mean = ' + str('{:.3f}'.format(params[1])) + ' $\pm$ ' + str('{:.3f}'.format(params[2])) + ' km/s',
            xy=(0.5, 0.5), xytext=(0, -160), 
            xycoords='axes fraction', textcoords='offset points', fontsize=11,
            ha='center', va='top')
    
    # Display the Gaussian parameters (Amplitude, Mean, and Standard Deviation)

    print("\nWeighted final velocity = ",
          "{} +/- {} km/s".format(params[1], params[2]))
        
    plt.savefig('Weighted_velocity_distribution.png', dpi=1000, bbox_inches='tight')
    
    return

get_data()
