# -*- coding: utf-8 -*-
"""
Created on Thu Nov 9 15:24:03 2023

@author: nickl

This code makes a user input a file (MAKE SURE YOU CHOOSE 'rings_final2.txt').
It takes the second and third columns (the Radius(arcsec) and VROT(KM/S) columns)
and plots them against a v prop r^(-0.5) rotation curve to estimate the mass of
a central body inhibiting Keplerian motion on nearby gas/dust. It also calculates
the chi squared using our own data in rings_final1.txt but since there is no
assosciated uncertainty, this value is very large
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tkinter import Tk, filedialog

########################### use variables ###################################
SaveFigName = 'C:/Users/nickl/linux/Week7/MassEstimates/k=7_bbarolo_curve_fit.png'
SaveFig = True
initial_dir = 'C:/Users/nickl/linux/prog/bbarolo/output/'
min_rad = 0  # minimum radius in arcsec
max_rad = 0.5  # maximum radius in arcsec (0.13"=1000AU)
#############################################################################

############################## global params ##################################
pc = 3.086e+16
G = 6.6743e-11
SunMass = 1.98847e30
scale_factor_x = 8 * 10**3 * np.pi * pc / (3600 * 180)  # turns arcsec into m
scale_factor_y = 1000  # turns km/s into m/s
#############################################################################

def get_data(initial_dir):
    # file picking with a temp window
    root = Tk()
    root.attributes('-topmost', True)  # Set the window to be always on top

    # Prompt the user to select files, starting in the specified directory
    file_path = filedialog.askopenfilename(initialdir=initial_dir, title="Select Folder, then select rings_final2.txt")

    root.destroy() # removes temp window

    file_path1 = file_path.replace('rings_final2.txt', 'rings_final1.txt')

    data_ours = np.genfromtxt(file_path1, usecols=(1, 2), skip_header=1, unpack=True)
    data_bbarolo = np.genfromtxt(file_path, usecols=(1, 2), skip_header=1, unpack=True)


    return data_bbarolo, data_ours

def validate(data):

    rad, vrot = data

    valid_indices = (rad > min_rad) & (rad < max_rad) # validation of data

    rad = rad[valid_indices]
    vrot = vrot[valid_indices]

    return rad, vrot

def model(x, A):
    return A * x**(-0.5)

def chi_squared(data, model, num_params):
    residuals = data - model
    #print(data, model)
    chi_squared = np.sum(residuals**2)
    degrees_of_freedom = len(data) - num_params
    reduced_chi_squared = chi_squared / degrees_of_freedom
    return reduced_chi_squared

def perform_curve_fitting(rad, vrot):

    initial_guess = [np.mean(vrot)] #initial mass

    # Fit the model to the data
    params, covariance = curve_fit(model, rad, vrot, p0=initial_guess, maxfev=10000)

    A = params[0]

    return A

def plot_initial(rad, vrot, A):
    # Create a finer array for plotting the smooth curve
    rad_fine = np.linspace(np.min(rad), np.max(rad), 1000)

    # Create a figure with two subplots, sharing the y-axis
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # First subplot: original data
    ax.plot(rad, vrot, 'o', color='blue', label='Rotation Velocity')
    ax.plot(rad_fine, model(rad_fine, A), '-', color='red', label=f'Fit: ${A:.3f} \cdot r^{{ -0.5}}$')
    ax.set_xlabel('Radius (arcsec)')
    ax.set_ylabel('Rotation Velocity (km/s)')
    ax.set_title('3DBbarolo Rotation Curve')
    ax.grid(alpha=0.3)
    ax.legend()

def scaling_to_si(rad, vrot, A):

    rad_scaled = rad * scale_factor_x
    vrot_scaled = vrot * scale_factor_y
    A_scaled = A * scale_factor_y/ scale_factor_x**(-0.5)

    return rad_scaled, vrot_scaled, A_scaled

def plot_curve(rad, vrot, rad1, vrot1, A, A_scaled):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    rad_fine = np.linspace(np.min(rad), np.max(rad), 1000)

    ax.plot(rad, vrot, 'o', color='indigo', label='Model data', alpha = 0.6)
    ax.plot(rad1, vrot1, 'o', color ='gold', label='Obseverved data')
    ax.plot(rad_fine, model(rad_fine, A_scaled), '-', color='indigo')

    SMass = A_scaled**2 / (G * SunMass)

    ax.plot([],[], label = f'M = {SMass:.2f} M$_\odot$', alpha = 0)
    ax.set_xlabel('Radius (m)')
    ax.set_ylabel('Rotation Velocity (m/s)')
    ax.set_title('SI Scaled Rotation Curve')
    ax.grid(alpha=0.3)
    ax.legend()

    print(f'The value for A in the scaled plot, equal to sqrt(GM), is {A_scaled:.2e}')
    print(f'\nThis gives a value A^2/G = M = {SMass:.3e} Solar Masses')

    if SaveFig:
        plt.savefig(fname=SaveFigName, bbox_inches='tight', dpi=866)
    # Display the plot
    plt.show()

def execute():

    data_bbarolo, data_ours = get_data(initial_dir)

    rad, vrot = validate(data_bbarolo) #returns rad,vrot
    rad1,vrot1 = validate(data_ours)

    # Plot the data and fit
    A = perform_curve_fitting(rad, vrot)
    plot_initial(rad, vrot, A)

    #scaling
    rad_s, vrot_s, A_s = scaling_to_si(rad, vrot, A)
    rad1_s, vrot1_s, A_unused = scaling_to_si(rad1, vrot1, A)

    plot_curve(rad_s, vrot_s, rad1_s, vrot1_s, A, A_s)

    # Calculate the reduced chi-squared value
    data = vrot1_s
    model_data = model(rad_s, A_s)
    chi_squared_red = chi_squared(data, model_data, 1)
    print(f'Reduced Chi-squared = {chi_squared_red:.4f}')

# Execute the code
execute()
