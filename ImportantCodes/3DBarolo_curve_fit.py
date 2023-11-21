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
SaveFigName = 'C:/Users/nickl/linux/Week7/MassEstimates/k=7_bbarolo_free_fit.png'
SaveFig = True
initial_dir = 'C:/Users/nickl/linux/prog/bbarolo/output/' #where is your rings_final2.txt file
min_rad = 0  # minimum radius in arcsec
max_rad = 0.5  # maximum radius in arcsec (0.13"=1000AU) #0.2 gets rid of last point
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

    data_ours = np.genfromtxt(file_path1, usecols=(1, 2, 3), skip_header=1, unpack=True)
    data_bbarolo = np.genfromtxt(file_path, usecols=(1, 2, 3), skip_header=1, unpack=True)


    return data_bbarolo, data_ours

def validate(data):

    rad, vrot, vdisp = data

    valid_indices = (rad > min_rad) & (rad < max_rad) # validation of data

    rad = rad[valid_indices]
    vrot = vrot[valid_indices]
    vdisp = vdisp[valid_indices]

    return rad, vrot, vdisp

def model(x, A):
    return A * x**(-0.5)

def non_keplarian_model(x, A2, n):
    return A2 * x**(n)
"""
def uncertainties():
    PA = 156 #155-158
    PA_err = 2
    inc = 34
    inc_err = 4
    crad = np.pi /180 #converts degrees to radians
    rad = 28.91
    rad_err = 0.24
    ###################
    xp = -x * np.cos()
    yp
    xp_err
    yp_err
"""
def chi_squared(data, model, uncert, num_params):
    #print(data, model)
    residuals = (data - model) / uncert  # Element-wise division
    chi_squared = np.sum(residuals**2)
    degrees_of_freedom = len(data) - num_params
    reduced_chi_squared = chi_squared / degrees_of_freedom
    return reduced_chi_squared

def perform_curve_fitting(rad, vrot, rad1, vrot1, vdisp, vdisp1):

    initial_guess_kep = [np.mean(vrot)]  # Initial guess for Keplerian model
    initial_guess_free = [np.mean(vrot), -0.5]  # Initial guess for non-Keplerian model (A2, n)

    # Fit the Keplerian model to the data
    params, cov = curve_fit(model, rad, vrot, p0=initial_guess_kep,
                                   sigma=vdisp, absolute_sigma=True, maxfev=10000)

    # Fit the non-Keplerian model to the data
    params2, cov2 = curve_fit(non_keplarian_model, rad, vrot, p0=initial_guess_free,
                                     sigma=vdisp1, absolute_sigma=True, maxfev=10000)

    params3, cov3 =curve_fit(model, rad1, vrot1, p0=initial_guess_kep,
                                     sigma=vdisp1, absolute_sigma=True, maxfev=10000)

    A = params[0]
    A_err = np.sqrt(cov[0, 0])
    A1, n = params2
    A1_err, n_err = np.sqrt(np.diag(cov2))[0],np.sqrt(np.diag(cov2))[1]
    A2 = params3[0]
    A2_err = np.sqrt(np.diag(cov3))[0]
    return A, A_err, A1, A1_err, A2, A2_err, n, n_err

def plot_initial(rad, vrot, vdisp, A, n, A2,n_u):
    # Create a finer array for plotting the smooth curve
    rad_fine = np.linspace(np.min(rad), np.max(rad), 1000)

    # Create a figure with two subplots, sharing the y-axis
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # First subplot: original data
    ax.errorbar(rad, vrot, yerr= vdisp/2, fmt='o', color='blue', label='Rotation Velocity')
    ax.plot(rad_fine, model(rad_fine, A), '-', color='red', label=f'Keplarian Fit: ${A:.3f} \cdot r^{{ -0.5}}$')
    #ax.plot(rad_fine, non_keplarian_model(rad_fine, A2, n), '-', color='green', label = f'n = {n:.3f}')
    ax.set_xlabel('Radius (arcsec)')
    ax.set_ylabel('Rotation Velocity (km/s)')
    ax.set_title('3DBbarolo Rotation Curve')
    ax.grid(alpha=0.3)
    ax.legend()

def scaling_to_si(rad, vrot, vdisp, A, A_err, A1, A2, A2_err, n):

    rad_scaled = rad * scale_factor_x
    vrot_scaled = vrot * scale_factor_y
    A_scaled = A * scale_factor_y/ scale_factor_x**(-0.5)
    A_err_s = A_err * scale_factor_y/ scale_factor_x**(-0.5)
    A1_scaled = A1 * scale_factor_y/ scale_factor_x**(n)
    A2_scaled = A2 * scale_factor_y/ scale_factor_x**(-0.5)
    A2_err_s = A2_err * scale_factor_y/ scale_factor_x**(-0.5)
    vdisp_scaled = vdisp* scale_factor_y/2

    return rad_scaled, vrot_scaled, vdisp_scaled, A_scaled, A_err_s, A1_scaled, A2_scaled, A2_err_s

def plot_curve(rad, vrot, vdisp, rad1, vrot1, vdisp1, A_scaled, A1_scaled, A2_scaled, n, n_err):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    rad_fine = np.linspace(np.min(rad), np.max(rad), 1000)
    SMass = A_scaled**2 / (G * SunMass)
    SMass2 = A2_scaled**2 / (G * SunMass)

    ax.errorbar(rad, vrot, yerr= vdisp, fmt= 'o', color='indigo')
    ax.errorbar(rad1, vrot1, yerr= vdisp1, fmt= 'o', color ='gold', label='Observed data')
    #ax.plot(rad1, vrot1, 'o', color ='gold', label='Observed data')
    ax.plot(rad_fine, model(rad_fine, A_scaled), '-', color='indigo', label = f'Model Data\nM = {SMass:.2f} M$_\odot$')
    ax.plot(rad_fine, non_keplarian_model(rad_fine, A1_scaled, n), '-', color='green', label = f'Free Fit\nn = {n:.3f} $\pm$ {n_err:.3f}')

    #SMass = A_scaled**2 / (G * SunMass)

    #ax.plot([],[], label = f'M = {SMass:.2f} M$_\odot$', alpha = 0)
    ax.set_xlabel('Radius (m)')
    ax.set_ylabel('Rotation Velocity (m/s)')
    ax.set_title('SI Scaled Rotation Curve')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

    print(f'The value for A in the scaled plot, equal to sqrt(GM), is {A_scaled:.2e}')

    if SaveFig:
        plt.savefig(fname=SaveFigName, bbox_inches='tight', dpi=866)
    # Display the plot
    plt.show()
    return SMass, SMass2

def execute():

    data_bbarolo, data_ours = get_data(initial_dir)

    rad, vrot, vdisp = validate(data_bbarolo) #returns rad,vrot
    rad1,vrot1, vdisp1 = validate(data_ours)

    #print(vdisp)
    # Plot the data and fit
    A, A_err, A1, A1_err, A2, A2_err, n, n_err = perform_curve_fitting(rad, vrot, rad1, vrot1, vdisp, vdisp1)
    plot_initial(rad, vrot, vdisp, A, n, A2, n_err)

    #scaling
    rad_s, vrot_s, vdisp_s, A_s, A_err_s, A1_unused, A2_unused, A2_err_s_unsused = scaling_to_si(rad, vrot, vdisp, A, A_err, A1, A2, A2_err, n)
    rad1_s, vrot1_s, vdisp1_s, A_unused, A_err_s_unused, A1_s, A2_s, A2_err_s = scaling_to_si(rad1, vrot1, vdisp1, A, A_err, A1, A2, A2_err, n)

    SMass, SMass2 = plot_curve(rad_s, vrot_s, vdisp_s, rad1_s, vrot1_s, vdisp1_s, A_s, A1_s, A2_s, n, n_err)
    SMass_err = ((2*A_s / G)**2 *(A_err_s)**2)**(1/2) / SunMass

    SMass2_err = ((2*A2_s / G)**2 *(A2_err_s)**2)**(1/2) / SunMass
    print(f'{A_err_s/A_s} is the % error in A, for {A_err} and {A_s}')
    print(f'{A2_err_s/A2_s} is the % error in A2, for {A2_err} and {A2_s}')

    ##########################
    # Calculate the reduced chi-squared value
    chi_squared_red = chi_squared(vrot_s, model(rad_s, A_s), vdisp_s, 1)
    chi_squared_red2 = chi_squared(vrot1_s, non_keplarian_model(rad1_s, A2_s, n), vdisp1_s, 1)
    #print(f'{A_s:.3f}, {A2_err_s:.3f}')

    print(f'Reduced Chi-squared of the 3Dbarolo model= {chi_squared_red:.4f}')
    print(f'This gives a value A^2/G = M = {SMass:.3g} +/- {SMass_err:.2f} Solar Masses')

    print(f'\nReduced Chi-squared of our data from free fit= {chi_squared_red2:.4f}')
    print(f'This gives a value A^2/G = M = {SMass2:.3e} +/- {SMass2_err:.2f}Solar Masses')


# Execute the code
execute()
