# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:48:08 2023

@author: nickl

ITS THE SAME AS NICKBAROLO.PY BUT WITH VELOCITY DISPERSIONS INCLUDED
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
SaveFigName1 = 'C:/Users/nickl/linux/s2/barolo/InitialCurve.png'
SaveFigName2 = 'C:/Users/nickl/linux/s2/barolo/Truedistrobut.png'
SaveFig = False
initial_dir = 'C:/Users/nickl/linux/s2/barolo/' #where is your rings_final2.txt file
min_rad = 0  # minimum radius in arcsec
max_rad = 1  # maximum radius in arcsec (0.13"=1000AU) #0.2 gets rid of last point
doi = 0.3 #times all the vdisps by this value
doiv = 2 #times all v_errs by this value
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

    data_ours = np.genfromtxt(file_path1, usecols=(1, 2, 3, 16), skip_header=1, unpack=True)
    data_bbarolo = np.genfromtxt(file_path, usecols=(1, 2, 3, 16), skip_header=1, unpack=True)


    return data_bbarolo, data_ours

def validate(data):

    rad, vrot, vdisp, vdisperr = data

    valid_indices = (rad > min_rad) & (rad < max_rad) # validation of data

    rad = rad[valid_indices]
    vrot = vrot[valid_indices]
    vdisp = vdisp[valid_indices]
    vdisperr = vdisperr[valid_indices]

    return rad, vrot, vdisp, vdisperr

def keplerian_model(x, A_kep):
    return A_kep * x**(-0.5)

def non_keplerian_model(x, A_nonkep, n):
    return A_nonkep * x**(n)

def perform_curve_fitting(rad_bar, vrot_bar, rad_our, vrot_our, vdisp_bar, vdisp_our, verr_bar,verr_our):

    initial_guess = [np.mean(vrot_bar)]
    initial_guess_free = [np.mean(vrot_bar), -0.5]

    # Fit the Keplerian model to the data
    Params_kep, cov_kep = curve_fit(keplerian_model, rad_bar, vrot_bar, p0=initial_guess, absolute_sigma = True, maxfev=10000)

    # Fit the non-Keplerian model to the data
    params_low, cov_low = curve_fit(keplerian_model, rad_bar, vrot_bar - vdisp_bar*doi, p0=initial_guess, absolute_sigma = True,maxfev=10000)

    params_high, cov_high =curve_fit(keplerian_model, rad_bar, vrot_bar + vdisp_bar*doi, p0=initial_guess, absolute_sigma = True,maxfev=10000)

    params_our, cov_our = curve_fit(non_keplerian_model, rad_bar, vrot_bar, p0=initial_guess_free, absolute_sigma=True, maxfev=10000)


    A_mid = Params_kep[0]
    A_mid_err = np.sqrt(cov_kep[0, 0])
    A_low = params_low[0]
    A_low_err = np.sqrt(cov_low[0, 0])
    A_high = params_high[0]
    A_high_err = np.sqrt(np.diag(cov_high))[0]
    A_free , n_free = params_our
    Aerr_free, nerr_free = np.sqrt(np.diag(cov_our))[0], np.sqrt(np.diag(cov_our))[1]
    print(n_free, nerr_free)
    return A_mid, A_mid_err, A_low, A_low_err, A_high, A_high_err, A_free, Aerr_free, n_free, nerr_free

def plot_initial(rad_bar, vrot_bar, vdisperr_bar, rad_our, vrot_our, vdisperr_our, A_mid, A_low, A_high, A_free, Aerr_free, n_free, nerr_free):
    # Create a finer array for plotting the smooth curve
    rad_fine = np.linspace(np.min(rad_bar), np.max(rad_bar), 1000)

    # Create a figure with two subplots, sharing the y-axis
    plt.figure(figsize=(10,6))

    # First subplot: original data
    plt.errorbar(rad_bar, vrot_bar, yerr= vdisperr_bar*(1/doiv), fmt='o', color='green', label='Model data')
    #plt.errorbar(rad_our, vrot_our, yerr= vdisperr_our*(1/doiv), fmt='o', color='indigo', label='Obs data')
    plt.plot(rad_fine, keplerian_model(rad_fine, A_mid), '-', color='green', label=f'Keplarian Fit: ${A_mid:.3f} \cdot r^{{ -0.5}}$')
    plt.plot(rad_fine, non_keplerian_model(rad_fine, A_free, n_free), '--', color='indigo', label=f'Model Fit: ${A_free:.3f} \cdot r^{{n_free}}$'+f'\nn = ${n_free:.3f} +/- {nerr_free:.3f}$')

    plt.xlabel('Radius (arcsec)', fontsize =12)
    plt.ylabel('Rotation Velocity (km/s)', fontsize = 12)
    plt.title('3DBbarolo Rotation Curve', fontsize = 15)
    plt.grid(alpha=0.3)
    plt.legend(fontsize=12)

    if SaveFig:
        plt.savefig(fname=SaveFigName1, bbox_inches='tight', dpi=800)

    plt.show()

"""
def chi_squared(data, model, uncert, num_params):
    #print(data, model)
    residuals = (data - model) / uncert  # Element-wise division
    chi_squared = np.sum(residuals**2)
    degrees_of_freedom = len(data) - num_params
    reduced_chi_squared = chi_squared / degrees_of_freedom
    return reduced_chi_squared
"""

def scaling_to_si1(rad, vrot, vdisp, vdisperr):

    rad_s = rad * scale_factor_x
    vrot_s = vrot * scale_factor_y
    vdisp_s = vdisp* scale_factor_y
    vdisperr_s = vdisperr * scale_factor_y

    return rad_s, vrot_s, vdisp_s, vdisperr_s,

def scaling_to_si2(A, n):
    A_s = A * scale_factor_y/ scale_factor_x**(n)
    return A_s

def plot_curve(rad, vrot, vrot_err, A_mid, A_mid_err, A_low, A_low_err, A_high, A_high_err, A_free, n_free, nerr_free):
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    rad_fine = np.linspace(np.min(rad), np.max(rad), 1000)

    ## mass calculations
    SMass = A_mid**2 / (G * SunMass)
    SMass_low = A_low**2 / (G * SunMass)
    SMass_high = A_high**2 / (G * SunMass)
    #SMass_err = ((2*A_ac / G)**2 *(A_ac_err)**2)**(1/2) / SunMass
    SMass_low_err = ((2*A_low / G)**2 *(A_low_err)**2)**(1/2) / SunMass
    SMass_high_err = ((2*A_high / G)**2 *(A_high_err)**2)**(1/2) / SunMass

    #### plot
    ax.errorbar(rad /AU, vrot/1000, yerr = vrot_err/1000*(1/doiv),fmt ='o', color='green', label = 'Model data')
    #ax.errorbar(rad1 /AU, vrot1/1000, yerr = vrot1_err/1000*(1/doiv),fmt ='o', color='indigo', label = 'Observed data')
    ax.plot(rad_fine/AU, non_keplerian_model(rad_fine, A_free, n_free)/1000, '--', color='indigo', label = f'Free fit n= {n_free:.2f} +/- {nerr_free:.2f}')
    ax.plot(rad_fine/AU, keplerian_model(rad_fine, A_mid)/1000, '-', color='green') #+/- {SMass_err:.2f}
    ax.plot(rad_fine/AU, keplerian_model(rad_fine, A_low)/1000, '--', color='red')
    ax.plot(rad_fine/AU, keplerian_model(rad_fine, A_high)/1000, '--', color='red')

    Low = SMass_low - SMass_low_err
    High = SMass_high + SMass_high_err
    ax.plot([], [], label=f'${Low:.3f}$ '+'$ M_{\odot}$'+f' < M < ${High:.3f}$'+' $M_{\odot}$', alpha = 0)
    print('Low: ', Low, '\nMid: ', SMass, '\nHigh: ', High)
    #plt.plot([], [], label = r'$M={:.2f}^{{+{:.2f}}}_{{-{:.2f}}}\,M_\odot$'.format(SMass, High, Low), alpha = 0)

    ax.set_xlabel('Radius (AU)', fontsize=12)
    ax.set_ylabel('Rotation Velocity (km/s)', fontsize=12)
    ax.set_title('Rotation Curve CH3CN: $J_{k}=12_{7}-11_{7}$', fontsize=15)
    ax.grid(alpha=0.3)
    ax.legend(fontsize = 12)

    if SaveFig:
        plt.savefig(fname=SaveFigName2, bbox_inches='tight', dpi=400)

    plt.show()



def execute():

    data_bbarolo, data_ours = get_data(initial_dir)

    rad_bar, vrot_bar, vdisp_bar, vdisperr_bar = validate(data_bbarolo) #returns rad,vrot
    rad_our, vrot_our, vdisp_our, vdisperr_our = validate(data_ours)

    #print(vdisp)
    # Plot the data and fit
    A_mid, A_mid_err, A_low, A_low_err, A_high, A_high_err, A_free, Aerr_free, n_free, nerr_free = perform_curve_fitting(rad_bar, vrot_bar, rad_our, vrot_our, vdisp_bar, vdisp_our, vdisperr_bar, vdisperr_our)
    plot_initial(rad_bar, vrot_bar, vdisperr_bar, rad_our, vrot_our, vdisperr_our, A_mid, A_low, A_high, A_free, Aerr_free, n_free, nerr_free)

    #scaling
    n_kep = -0.5
    rad_s, vrot_s, vdisp_s, vdisperr_s = scaling_to_si1(rad_bar, vrot_bar, vdisp_bar, vdisperr_bar)
    A_mid_s= scaling_to_si2(A_mid, n_kep)
    A_mid_err_s = scaling_to_si2(A_mid_err, n_kep)
    A_free_s = scaling_to_si2(A_free, n_free)
    #print(f'{A_mid_s*(1/1000)*AU**-0.5:.3f}, in km/s AU**0.5')
    A_low_s= scaling_to_si2(A_low, n_kep)
    A_low_err_s = scaling_to_si2(A_low_err, n_kep)
    A_high_s= scaling_to_si2(A_high, n_kep)
    A_high_err_s = scaling_to_si2(A_high_err, n_kep)

    #reduced_chi_squared = chi_squared(vrot1_s, model(rad_s, A_ac_s), vrot1_err_s*(1/doiv), 1)
    #print(f'{reduced_chi_squared:.3f}')
    plot_curve(rad_s, vrot_s, vdisperr_s, A_mid_s, A_mid_err_s, A_low_s, A_low_err_s, A_high_s, A_high_err_s,
               A_free_s, n_free, nerr_free)

# Execute the code
execute()
