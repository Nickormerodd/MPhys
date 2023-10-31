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
#import uncertainties.unumpy as unp
#from uncertainties import ufloat
from scipy.stats import norm

G = 6.67430**(-11)
AU = 1.496*10**11

def get_data():

    Tk().withdraw()
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        filename = (os.path.basename(file_path).replace(".fits", ""))#.replace("_"," ").replace("."," ")
        main(data, filename, file_path)
        data.close()

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

    mom = data[0].data[0]
    v_0 = np.nanmean(mom)
    print(v_0)
    temp = []  # Create an empty list to store the values
    for i in range(mom.shape[0]):  # Loop through rows
        for j in range(mom.shape[1]):  # Loop through columns
            value = mom[i, j]  # Get the value at the current position
            if not (np.isnan(value).any() or value == 0):
                # Calculate distance using x_au and y_au arrays
                dist = np.sqrt((y_au[i] - y_o)**2 + (x_au[j] - x_o)**2)
                temp.append((np.abs(dist), (value-v_0)))  # Store the distance and value

                # Sort the list based on distances
    temp.sort(key=lambda x: x[0])

    # Now, temp contains a list of (distance, value) pairs, sorted by distance.

    distances, velocity = zip(*temp)

    # Convert distances to AU if needed
    # For example, if your distances are in parsecs, you can convert to AU using: distances = np.array(distances) * 206264.81 * 8.1 * 10**3

    combined = np.hstack((np.vstack((distances)), np.vstack(velocity)))

    pos = []
    neg = []
    for line in combined:
        if line[1] < 0:
            line[1] = np.abs(line[1])
            neg.append(line)
        if line[1] >= 0:
            line[1] = np.abs(line[1])
            pos.append(line)
    pos,neg = np.vstack((pos)), np.vstack((neg))

    # Create a scatter plot
    #ax.scatter(distances, velocity)
    ax.scatter(pos[:,0],pos[:,1], label = '+ve', alpha=0.7, s=15)
    ax.scatter(neg[:,0],neg[:,1], label = '-ve', alpha=0.7, s=15)
    ax.set_xlabel('Distance AU')
    ax.set_ylabel('Velocity, km/s')
    ax.set_title('Velocity-Distance for k=7')

    combined[:,1] = np.abs(combined[:,1])

    new_pos = filtered(pos)
    new_neg = filtered(neg)
    new_comb = filtered(combined)
    a_pos = curve(new_pos)
    a_neg = curve(new_neg)
    a_comb = curve(new_comb)
    M_pos = solar_mass(a_pos)
    M_neg = solar_mass(a_neg)
    M_comb = solar_mass(a_comb)


    ax.plot(new_pos[:,0], keplarian(new_pos[:,0],a_pos),
            label = r'+ve: $M_{\odot}$= ' + str('{:.3f}'.format(M_pos)),
            color = 'green', alpha=0.5)
    ax.plot(new_neg[:,0], keplarian(new_neg[:,0],a_neg),
            label = r'-ve: $M_{\odot}$= ' + str('{:.3f}'.format(M_neg)),
            color = 'red', alpha=0.5)
    ax.plot(new_comb[:,0], keplarian(new_comb[:,0],a_comb),
            label = r'Both: $M_{\odot}$= ' + str('{:.3f}'.format(M_comb)),
            color = 'purple', alpha=0.5)

    ax.legend(fontsize=8, borderaxespad=0.5, frameon=False, edgecolor='black')

    return

def keplarian(distances, a):

    v = a * distances

    return v

def filtered(data):

    piglords = []
    for line in data:
        if not line[0] > 1500:
            piglords.append(line)

    new_data = np.vstack((piglords))

    return new_data

def curve(data):

    GUESS = None
    n=0
    while n < 5:
        try:
            params, covariance = curve_fit(keplarian, data[:,0], data[:,1],
                                           p0=GUESS, absolute_sigma=True,
                                           maxfev=2000)
        except ValueError or RuntimeWarning:
            break
        GUESS = params
        n+=1

    a = (params[0])

    return a

def solar_mass(a):

    Nistance_AU = 1500 * AU
    M_sol = a**2 * Nistance_AU**2/ (G * 1.988*10**30)

    return M_sol

get_data()
