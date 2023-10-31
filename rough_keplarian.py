"""
Created on Tue Oct 31 15:24:56 2023

@author: Christopher
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

G = 6.67430e-11

def get_data():

    Tk().withdraw()  
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        filename = (os.path.basename(file_path).replace(".fits", ""))#.replace("_"," ").replace("."," ")
        main(data, filename, file_path)
        data.close()   
    
    #four_plot(file_paths)
    
    return


def main(data, filename, path):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    header = data[0].header
    
    ctype1 = header['CTYPE1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    crpix1 = header['CRPIX1']
    
    ctype2 = header['CTYPE2']
    crval2 = header['CRVAL2']
    cdelt2 = header['CDELT2']
    crpix2 = header['CRPIX2']
    calibration_freq = 220.70908 #GHz
    
    k_number = 7
    conversion = (c/1000) * (220.538635 - calibration_freq)/calibration_freq
    data[0].data[0] = data[0].data[0] + conversion
    
    # Create the coordinate arrays
    x = (crval1 + (np.arange(data[0].shape[2]) - crpix1) * cdelt1)
    y = (crval2 + (np.arange(data[0].shape[1]) - crpix2) * cdelt2)
    
    x_au = x * 8.1 * 10**3 * 206264.81 * np.pi/180
    y_au = y * 8.1 * 10**3 * 206264.81 * np.pi/180
    
    y_o = np.mean(y_au)
    x_o = np.mean(x_au)
    """
    temp = []
    for line in data[0].data:
        for i in line:
            if i == np.nan or i == 0:
                pass
            else:
                dist = ((y_au[line] - y_o)**2 + (x_au[i] - x_o)**2)**(1/2)
                temp = temp.append(i)
            
    """
    mom = data[0].data[0]
    #print(mom)
    
    v_0 = np.nanmean(mom)
    print(v_0)
    temp = []  # Create an empty list to store the values
    for i in range(mom.shape[0]):  # Loop through rows
        for j in range(mom.shape[1]):  # Loop through columns
            value = mom[i, j]  # Get the value at the current position
            if not (np.isnan(value).any() or value == 0):
                # Calculate distance using x_au and y_au arrays
                dist = np.sqrt((y_au[i] - y_o)**2 + (x_au[j] - x_o)**2) 
                temp.append((np.abs(dist), np.abs(value-v_0)))  # Store the distance and value
                
                # Sort the list based on distances
    temp.sort(key=lambda x: x[0])

    # Now, temp contains a list of (distance, value) pairs, sorted by distance.

    distances, velocity = zip(*temp)

    # Convert distances to AU if needed
    # For example, if your distances are in parsecs, you can convert to AU using: distances = np.array(distances) * 206264.81 * 8.1 * 10**3
    
    # Create a scatter plot
    ax.scatter(distances, velocity)
    ax.set_xlabel('Distance AU')
    ax.set_ylabel('Velocity, km/s')
    ax.set_title('Velocity-Distance')
    #ax.grid(True)
    #plt.show()
    
    GUESS = None
    n=0
    while n < 5:
        try:
            params, covariance = curve_fit(keplarian, distances, velocity,
                                           p0=GUESS, absolute_sigma=True,
                                           maxfev=2000)
        except ValueError or RuntimeWarning:
            break
        GUESS = params
        n+=1
        
    a = params[0]
    print(a)
    M_sol = a**2 / (G * 1.988*10**30)
    #M_sol_pig = '{:.3f}'.format(M_sol)
    
    
    #print(vel)
    
    ax.plot(dist, keplarian(distances,a), label = r'$M_{odot}$= ' + str('{:.3f}'.format(M_sol)))
    
    ax.legend()
    
    return

def keplarian(dist, a):
    
    v = a * dist**(-0.5)
    
    return v

get_data()
