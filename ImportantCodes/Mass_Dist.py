# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:48:08 2023

@author: nickl
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 9 15:24:03 2023

@author: nickl

This code Calculates a mass distribution using velocity dispersions getting
a high and low limit on the protostar mass
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tkinter import Tk, filedialog

########################### use variables ###################################
SaveFigName = 'C:/Users/nickl/linux/Week7/MassEstimates/MassDist.png'
SaveFig = False
initial_dir = 'C:/Users/nickl/linux/prog/bbarolo/output/' #where is your rings_final2.txt file
min_rad = 0.03  # minimum radius in arcsec
max_rad = 1  # maximum radius in arcsec (0.13"=1000AU) #0.2 gets rid of last point
doi = 1/2 #times all the vdisps by this value
#############################################################################

############################## global params ##################################
pc = 3.086e+16
G = 6.6743e-11
SunMass = 1.98847e30
scale_factor_x = 8 * 10**3 * np.pi * pc / (3600 * 180)  # turns arcsec into m
scale_factor_y = 1000  # turns km/s into m/s
AU = 1.496e+11
#############################################################################

def get_data(initial_dir):
    # file picking with a temp window
    root = Tk()
    root.attributes('-topmost', True)  # Set the window to be always on top

    # Prompt the user to select files, starting in the specified directory
    file_path = filedialog.askopenfilename(initialdir=initial_dir, title="Select Folder, then select rings_final2.txt")

    root.destroy() # removes temp window

    file_path1 = file_path.replace('rings_final2.txt', 'rings_final1.txt')

    data_ours = np.genfromtxt(file_path1, usecols=(1, 2, 3, 14), skip_header=1, unpack=True)
    data_bbarolo = np.genfromtxt(file_path, usecols=(1, 2, 3,14), skip_header=1, unpack=True)


    return data_bbarolo, data_ours

def validate(data):

    rad, vrot, vdisp, verr = data

    valid_indices = (rad > min_rad) & (rad < max_rad) # validation of data

    rad = rad[valid_indices]
    vrot = vrot[valid_indices]
    vdisp = vdisp[valid_indices]
    verr = verr[valid_indices]

    return rad, vrot, vdisp, verr

def model(x, A):
    return A * x**(-0.5)

def non_keplarian_model(x, A2, n):
    return A2 * x**(n)

def perform_curve_fitting(rad, vrot, rad1, vrot1, vdisp, vdisp1, vrot_err,vrot1_err):

    initial_guess = [np.mean(vrot)]
    initial_guess_free = [np.mean(vrot), -0.5]

    # Fit the Keplerian model to the data
    Params, cov = curve_fit(model, rad, vrot, p0=initial_guess, sigma=vrot_err,
                            absolute_sigma = True, maxfev=10000)

    # Fit the non-Keplerian model to the data
    params2, cov2 = curve_fit(model, rad, vrot - vdisp*doi, p0=initial_guess,
                              sigma=vrot_err, absolute_sigma = True,maxfev=10000)

    params3, cov3 =curve_fit(model, rad1, vrot + vdisp*doi, p0=initial_guess,
                             sigma=vrot_err, absolute_sigma = True,maxfev=10000)

    params4, cov4 = curve_fit(non_keplarian_model, rad, vrot, p0=initial_guess_free,
                                     sigma=vrot1_err, absolute_sigma=True, maxfev=10000)

    A_ac = Params[0]
    A_ac_err = np.sqrt(cov[0, 0])
    A_low = params2[0]
    A_low_err = np.sqrt(cov2[0, 0])
    A_high = params3[0]
    A_high_err = np.sqrt(np.diag(cov3))[0]
    A_f,n = params4
    A_f_err, n_err = np.sqrt(np.diag(cov4))[0],np.sqrt(np.diag(cov4))[1]
    return A_ac, A_ac_err, A_low, A_low_err, A_high, A_high_err, A_f, A_f_err, n, n_err

def plot_initial(rad, vrot, vrot_err, rad1, vrot1, vrot1_err, A_ac, A_low, A_high, A_f, A_f_err, n, n_err):
    # Create a finer array for plotting the smooth curve
    rad_fine = np.linspace(np.min(rad), np.max(rad), 1000)

    # Create a figure with two subplots, sharing the y-axis
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # First subplot: original data
    ax.errorbar(rad, vrot, yerr= vrot_err*doi, fmt='o', color='green', label='Model data')
    ax.errorbar(rad1, vrot1, yerr= vrot1_err*doi, fmt='o', color='red', label='Obs data')
    ax.plot(rad_fine, model(rad_fine, A_ac), '--', color='green', label=f'Keplarian Fit: ${A_ac:.3f} \cdot r^{{ -0.5}}$')
    ax.plot(rad_fine, non_keplarian_model(rad_fine, A_f, n), '--', color='red', label=f'Model Fit: ${A_f:.3f} \cdot r^{{n}}$'+f'\nn = ${n:.3f}$')

    ax.set_xlabel('Radius (arcsec)')
    ax.set_ylabel('Rotation Velocity (km/s)')
    ax.set_title('3DBbarolo Rotation Curve')
    ax.grid(alpha=0.3)
    ax.legend()

def scaling_to_si1(rad, vrot, vdisp, vrot_err, A_f, A_f_err, n):

    rad_s = rad * scale_factor_x
    vrot_s = vrot * scale_factor_y
    vdisp_s = vdisp* scale_factor_y
    vrot_err_s = vrot_err * scale_factor_y
    A_f_s = A_f * scale_factor_y/ scale_factor_x**(n)
    A_f_err_s = A_f_err * scale_factor_y/ scale_factor_x**(n)

    return rad_s, vrot_s, vdisp_s, vrot_err_s, A_f_s, A_f_err_s

def scaling_to_si2(data):
    data_s = data * scale_factor_y/ scale_factor_x**(-0.5)
    return data_s

def plot_curve(rad, vrot, vrot_err, rad1, vrot1, vrot1_err, A_ac, A_ac_err, A_low, A_low_err, A_high, A_high_err, A_f, A_f_err, n, n_err):
    plt.figure(figsize=(10,6))
    rad_fine = np.linspace(np.min(rad), np.max(rad), 1000)
    SMass = A_ac**2 / (G * SunMass)
    SMass_low = A_low**2 / (G * SunMass)
    SMass_high = A_high**2 / (G * SunMass)

    SMass_err = ((2*A_ac / G)**2 *(A_ac_err)**2)**(1/2) / SunMass
    SMass_low_err = ((2*A_low / G)**2 *(A_low_err)**2)**(1/2) / SunMass
    SMass_high_err = ((2*A_high / G)**2 *(A_high_err)**2)**(1/2) / SunMass

    plt.errorbar(rad /AU, vrot/1000, yerr = vrot_err/1000,fmt ='o', color='green', label = 'Model data')
    plt.errorbar(rad1 /AU, vrot1/1000, yerr = vrot1_err/1000,fmt ='o', color='indigo', label = 'Observed data')
    plt.plot(rad_fine/AU, non_keplarian_model(rad_fine, A_f, n)/1000, '--', color='indigo', label = f'Free fit n= {n:.2f} +/- {n_err:.2f}')
    plt.plot(rad_fine/AU, model(rad_fine, A_ac)/1000, '--', color='green', label = f'model M={SMass:.2f} +/- {SMass_err:.2f}'+' $M_{\cdot}$')
    plt.plot(rad_fine/AU, model(rad_fine, A_low)/1000, '-', color='red', label = f'low M={SMass_low:.2f} +/- {SMass_low_err:.2f}'+' $M_{\cdot}$')
    plt.plot(rad_fine/AU, model(rad_fine, A_high)/1000, '-', color='blue', label = f'high M={SMass_high:.2f} +/- {SMass_high_err:.2f}'+' $M_{\cdot}$')

    plt.xlabel('Radius (AU)')
    plt.ylabel('Rotation Velocity (km/s)')
    plt.title('SI Scaled Rotation Curve')
    plt.grid(alpha=0.3)
    plt.legend()

    if SaveFig:
        plt.savefig(fname=SaveFigName, bbox_inches='tight', dpi=866)

    plt.show()

def execute():

    data_bbarolo, data_ours = get_data(initial_dir)

    rad, vrot, vdisp, vrot_err = validate(data_bbarolo) #returns rad,vrot
    rad1,vrot1, vdisp1, vrot1_err = validate(data_ours)

    #print(vdisp)
    # Plot the data and fit
    A_ac, A_ac_err, A_low, A_low_err, A_high, A_high_err, A_f, A_f_err, n, n_err= perform_curve_fitting(rad, vrot, rad1, vrot1, vdisp, vdisp1, vrot_err, vrot1_err)
    plot_initial(rad, vrot, vrot_err, rad1, vrot1, vrot1_err, A_ac, A_low, A_high, A_f, A_f_err, n, n_err)

    #scaling
    rad_s, vrot_s, vdisp_s, vrot_err_s, U1, U2 = scaling_to_si1(rad, vrot, vdisp, vrot_err, A_f, A_f_err, n)
    rad1_s, vrot1_s, vdisp1_s, vrot1_err_s, A_f_s, A_f_err_s= scaling_to_si1(rad1, vrot1, vdisp1, vrot1_err, A_f, A_f_err, n)
    A_ac_s= scaling_to_si2(A_ac)
    A_ac_err_s = scaling_to_si2(A_ac_err)
    A_low_s= scaling_to_si2(A_low)
    A_low_err_s = scaling_to_si2(A_low_err)
    A_high_s= scaling_to_si2(A_high)
    A_high_err_s = scaling_to_si2(A_high_err)

    plot_curve(rad_s, vrot_s, vrot_err_s, rad1_s, vrot1_s, vrot1_err_s, A_ac_s,
               A_ac_err_s, A_low_s, A_low_err_s, A_high_s, A_high_err_s,
               A_f_s, A_f_err_s, n, n_err)

# Execute the code
execute()

