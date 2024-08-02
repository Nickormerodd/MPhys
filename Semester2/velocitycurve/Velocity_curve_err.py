# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 12:52:59 2023
Edited June 2024

@author: Christopher Andrews

Using a velocity map, this code creates a rotational velocity vs distance
figure and calculates the exponent and mass, with respective uncertainties.

It also calculates a ring-averaged velocity distance. The uncertainties
on each point here are calculated using MAD, which may be an overestimate. 
The code treats this calculated 'uncertainty' on the velocity like a true
uncertainty, but it is really just the distribution of points due to scatter
around the ring. 

The coordinate transformation implementation should be reconsidered,
maybe using an astropy package etc. 

User-defined inputs are at the top.

"""

import numpy as np
from astropy.io import fits
from tkinter import Tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import copy
from statsmodels import robust

G = 6.6743e-11
au = 1.495979e+11
M_o = 1.98847e+30
OUT = 'C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/S2/W4-6/Velocity_curve/'

##############################################################################
#User defined inputs

XO_PIX,YO_PIX = 17,17 # Define central pixel to normalise around

# Radial component to trim the data set to.
TRIM = 1750 # AU
LOWER_TRIM = 250 # AU - allows for better fits 

# Show the rotational figure:
SHOW_ROT_FIG = True

# Show the exponent vs position angle figure:
SHOW_N_PA = True

# Position angle + inclination angle + uncertainties in degrees
POSITION_ANGLE = 170 
INCLINATION = 38
PHI_ERR = 3 # position angle uncertainty
INC_ERR = 1.5 # inclination angle uncertainty

# Radial velocity of source + uncertainty
RADIAL_VELOCITY = 28.9 #km/s
V_RAD_ERR = 0.0707770409 #km/s



##############################################################################

def get_data():

    Tk().withdraw()
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        filename = (os.path.basename(file_path)).replace(".fits", "")#.replace(".image.pbcor_line.galactic","")
        main(data, filename, file_path)
        data.close()
        
    return

def main(data,name,file_path):

    xo_pix,yo_pix = XO_PIX,YO_PIX

    header = data[0].header
    
    # Masking velocities near zero to remove non-real zero rotational velocites
    # which skew the fittings.
    valid_indices = (data[0].data < (np.nanmean(data[0].data - 0.2))) | (data[0].data > (np.nanmean(data[0].data + 0.2)))
    data[0].data[~valid_indices] = np.nan

    data[0].data = data[0].data - np.nanmean(data[0].data)

    amended_data = match_shapes_with_nan(data[0].data)[0]  # This selects the 2D slice.

    x_au, y_au = galactic_coords(amended_data, header)
    
    xo_au = x_au[xo_pix] #central xo point
    yo_au = y_au[yo_pix] #central yo point
    
    x_norm = x_au - xo_au #normalising around the xo point
    y_norm = y_au - yo_au #normalising around the yo point

    main_5(amended_data, x_norm, y_norm, name)

    return

def main_5(data, x_norm,y_norm, name):
    
    incerz = np.linspace(np.deg2rad(34), np.deg2rad(34), 1)
    posizers = np.linspace(np.deg2rad(120), np.deg2rad(190), 500)
    phi_err = np.deg2rad(PHI_ERR)
    i_err = np.deg2rad(INC_ERR)
    
    i_mesh, pa_mesh = np.meshgrid(incerz, posizers, indexing='ij')

    # Initialize an empty list to store the results
    results_array = np.full((len(incerz), len(posizers), 5), np.nan)
    
###############################################################################
    # This is to make the position angle vs exponent plot - not necessarily needed
    # Iterate through the meshgrid
    for i in range(len(incerz)):
        for j in range(len(posizers)):
            inc = i_mesh[i, j]
            phi = pa_mesh[i, j]

    
            r_err, r, theta_err, theta = x_y_uncertainty(x_norm, y_norm, phi, phi_err)
            v_r_err, v_r = velocity_uncertainty(data, 28.9,inc, theta, theta_err, i_err)
            #v_r = rotational_velocity_real(data, x_norm, y_norm)
            r_flat = r.flatten()
            v_r_flat = v_r.flatten()
            v_r_err_flat = v_r_err.flatten()
            rv_array = np.vstack((r_flat,v_r_flat,v_r_err_flat)).T
            
            rv_array = no_nan(rv_array)
            if len(rv_array) == 0:
                continue
            
            new_array = trimmerz(rv_array, 250, 1750)
            if len(new_array) == 0:
                continue

            try:
                best_params, covariance = curve_fit(power_law, new_array[:,0], new_array[:,1],
                                                    sigma = new_array[:,2], maxfev=5000,
                                                    absolute_sigma=True)
                
                A_4, n_4 = best_params
                n_perr = np.sqrt(np.diag(covariance))[1]

            except RuntimeError:
                A_4, n_4 = np.nan, np.nan
    
            results_array[i, j] = [A_4, n_4, inc, phi, n_perr]
    
    # To make the n-exponent vs position angle
    A_closest, closest_n, inclination, position_angle, n_uncert = find_closest_n(results_array)
    
##############################################################################
    
    posas = np.deg2rad(POSITION_ANGLE)
    incatrons = np.deg2rad(INCLINATION)

    # Calculating the all the values
    r_err, r, theta_err, theta = x_y_uncertainty(x_norm, y_norm, position_angle, phi_err)
    v_r_err, v_r = velocity_uncertainty(data, RADIAL_VELOCITY, incatrons, theta, theta_err, i_err)
    
    # To show a view of the 'rotational velocity map' to see if it looks ok
    if SHOW_ROT_FIG:
        fig, ax = plt.subplots()
        cmap = plt.cm.viridis  
        cmap.set_bad(color='black')  
        im = ax.imshow(np.ma.masked_invalid(v_r),cmap=cmap)
        ax.set_title('Rotational Volocity')
        plt.colorbar(im) 
        plt.show()
    
    # Flatten matrices
    r_flat = r.flatten()
    v_r_flat = v_r.flatten()
    r_err_flat = r_err.flatten()
    v_r_err_flat = v_r_err.flatten()
    
    # Combine into one matrix
    result = np.vstack((r_flat,v_r_flat,r_err_flat,v_r_err_flat)).T
    result = no_nan(result) # Remove nans
    result = trimming(result) # Trim data such that radial data is below TRIM
    result = result[result[:, 0].argsort()]
    
    # Result_norm has no bounds applied.
    result_norm = trimmerz(result, LOWER_TRIM, np.inf)
    result_norm = result_norm[result_norm[:, 0].argsort()]
    
    # Curve fitting for whole data
    param, cov = fitting_data(power_law, result_norm, 0)
    params_k, covariance_k = fitting_data(keplarian_fit,result_norm, 0)
    
    # Result_low is trimming 250-inf and has lower-bound velocity: v - v_err
    result_low = trimmerz(result, LOWER_TRIM, np.inf)
    result_low[:,1] = result_low[:,1] - result_low[:,3]
    result_low = result_low[result_low[:, 0].argsort()]
    
    # Curve fitting for result_low
    param_low, cov_low = fitting_data(power_law, result_low, 0)
    params_k_low, covariance_k_low = fitting_data(keplarian_fit,result_low, 0)
    
    # Result_high is trimming 250-inf and has upper-bound velocity: v + v_err
    result_high = trimmerz(result, LOWER_TRIM, np.inf)
    result_high[:,1] = result_high[:,1] + result_high[:,3]
    result_high = result_high[result_high[:, 0].argsort()]
    
    # Curve fitting for result_high
    param_high, cov_high = fitting_data(power_law, result_high, 0)
    params_k_high, covariance_k_high = fitting_data(keplarian_fit,result_high, 0)
    
    # Ring-averaging the rotational velocity
    averaged_array = ring_averaging_with_mad(result)
    
    # Copy array and modify so it is the lower bound
    averaged_array_low = copy.deepcopy(averaged_array)
    averaged_array_low[:,1] = averaged_array_low[:,1] - averaged_array_low[:,3]
    
    # Copy array and modify so it is the upper bound
    averaged_array_high = copy.deepcopy(averaged_array)
    averaged_array_high[:,1] = averaged_array_high[:,1] + averaged_array_high[:,3]
    
    # Curve fitting for the ring-averaged array, lower bound + upper bound
    param_av, cov_av = fitting_data(power_law, trimmerz(averaged_array, LOWER_TRIM, np.inf), 0)
    params_k_av, covariance_k_av = fitting_data(keplarian_fit,trimmerz(averaged_array, LOWER_TRIM, np.inf), 0)
    
    param_av_low, cov_av_low = fitting_data(power_law, trimmerz(averaged_array_low, LOWER_TRIM, np.inf), 0)
    params_k_av_low, covariance_k_av_low = fitting_data(keplarian_fit,trimmerz(averaged_array_low, LOWER_TRIM, np.inf), 0)
    
    param_av_high, cov_av_high = fitting_data(power_law, trimmerz(averaged_array_high, LOWER_TRIM, np.inf), 0)
    params_k_av_high, covariance_k_av_high = fitting_data(keplarian_fit,trimmerz(averaged_array_high, LOWER_TRIM, np.inf), 0)
    
    
    # Normal plotting
    plotting(result, param, cov, params_k, covariance_k, posas, 
                             INCLINATION, param_low, cov_low, params_k_low, 
                             covariance_k_low, param_high, cov_high, params_k_high, 
                             covariance_k_high, name, 0)
    
    # Ring-averaging plotting
    plotting(averaged_array, param_av, cov_av, params_k_av, covariance_k_av, 
             posas, INCLINATION, param_av_low, cov_av_low, params_k_av_low, 
             covariance_k_av_low, param_av_high, cov_av_high, params_k_av_high, 
             covariance_k_av_high, name, 1)
    
    if SHOW_N_PA:
        n_PA(results_array, position_angle)
    
    return

def ring_averaging_with_mad(matrix, bin_size=100):
    """
    Errors for this are still massive, so might want to use another method as 
    the mass-range is very large.
    
    Perform ring averaging using the Median Absolute Deviation (MAD) to calculate uncertainties.
    
    Parameters:
    - matrix: A 2D numpy array where the first column is the radial component in AU,
              the second column is the rotational velocity component,
              and the third column is the radial error (which we don't use for averaging here),
              and the fourth column is the velocity error component.
    - bin_size: Size of each bin in AU.

    Returns:
    - A new 2D numpy array with each row representing the average radial point,
      the average rotational velocity, a placeholder for radial error, 
      and the MAD of the rotational velocity for each 100 AU segment.
    """
    if matrix.shape[1] > 4:
        matrix = matrix[:, :4]

    max_distance = np.max(matrix[:, 0])
    bins = np.arange(0, max_distance + bin_size, bin_size)

    averaged_matrix = []
    mad_matrix = []

    for i in range(1, len(bins)):
        bin_mask = (matrix[:, 0] >= bins[i-1]) & (matrix[:, 0] < bins[i])
        if np.any(bin_mask):
            bin_data = matrix[bin_mask, :]
            
            mean_values = np.mean(bin_data, axis=0)
            mad_values = robust.mad(bin_data[:, 1])  # MAD for the rotational velocity component
            
            mean_values[3] = mad_values  # Update the velocity error with MAD
            
            averaged_matrix.append(mean_values)

    averaged_matrix_np = np.vstack(averaged_matrix) if averaged_matrix else np.empty((0, 4))

    return averaged_matrix_np


##################################################################################################
# This isnt used, but have left it here. The velocity error is the standard deviation
def ring_averaging(matrix):
    """
    May want to modify this so it takes into account weighted average or 
    remove outliers from the average etc.
    
    Perform averaging of the radial velocity component, rotational velocity component,
    and a third data component over every 100 AU, excluding any further columns.
    Add standard deviation for each point.

    Parameters:
    - matrix: A 2D numpy array where the first column is the radial component in AU,
              the second column is the rotational velocity component,
              and the third column is another data component.

    Returns:
    - A new 2D numpy array with each row representing the average radial point,
      the average rotational velocity, and the average of the third data component
      for each 100 AU segment.
    - A new 2D numpy array with each row representing the standard deviation of the radial point,
      the standard deviation of the rotational velocity, and the standard deviation of the third data component
      for each 100 AU segment.
    """
    # Ensure matrix is limited to the first three columns if more exist
    if matrix.shape[1] > 3:
        matrix = matrix[:, :3]

    # Define the bin size for averaging (100 AU)
    bin_size = 100

    # Determine the maximum radial distance to establish bins
    max_distance = np.max(matrix[:, 0])
    bins = np.arange(0, max_distance + bin_size, bin_size)

    # Initialize matrices to store the averaged values and standard deviations
    averaged_matrix = []
    std_matrix = []

    # Loop through each bin and calculate averages and standard deviations
    for i in range(1, len(bins)):
        bin_mask = (matrix[:, 0] >= bins[i-1]) & (matrix[:, 0] < bins[i])
        if np.any(bin_mask):
            mean_values = np.mean(matrix[bin_mask, :], axis=0)
            std_values = np.std(matrix[bin_mask, :], axis=0)
            averaged_matrix.append(mean_values)
            std_matrix.append(std_values)

    # Convert lists to numpy arrays for output
    averaged_matrix_np = np.vstack(averaged_matrix) if averaged_matrix else np.empty((0, 3))
    std_matrix_np = np.vstack(std_matrix) if std_matrix else np.empty((0, 3))
    
    averaged_matrix_nice = np.hstack((averaged_matrix_np, std_matrix_np))

    return averaged_matrix_nice

################################################################################################

def coord_trans(x,y,phi):
    """
    ---------------------------------------------------------------------------
    Have left this in, but should do this coord trans more carefully using 
    astropy packages, and more carefully consider the coord units used etc.
    ---------------------------------------------------------------------------

    """
    inc = np.deg2rad(INCLINATION)
    alpha = phi # +/- np.pi
    
    x_prime = x*np.cos(alpha) - y*np.sin(alpha)
    y_prime = x*np.sin(alpha)/np.cos(inc) + y*np.cos(alpha)/np.cos(inc)

    return x_prime,y_prime

def plotting(array,param,cov,params_k,covariance_k,phi,inc, param_low, cov_low,
             params_k_low, covariance_k_low, param_high, cov_high, params_k_high, 
             covariance_k_high, name, ring):
    
    phi = '{:.0f}'.format(float(phi))
    inc = '{:.0f}'.format(float(inc))

    n = param[1]
    n_err = np.sqrt(np.diag(cov))[1]
    
    # Proportionality constant and mass
    a_scaled = params_k[0] * 10**3 * au**0.5
    a_err_scaled = np.sqrt(np.diag(covariance_k))[0] * 10**3 * au**0.5
    mass = '{:.2f}'.format(a_scaled**2 / (G*M_o))
    mass_err = '{:.2f}'.format((2*a_scaled*a_err_scaled/G) / (M_o))
    
    plt.figure(figsize=(10, 6))

    # Lower bound
    n_low = param_low[1]
    n_err_low = np.sqrt(np.diag(cov_low))[1]
    a_scaled_low = params_k_low[0] * 10**3 * au**0.5
    a_err_scaled_low = np.sqrt(np.diag(covariance_k_low))[0] * 10**3 * au**0.5
    mass_low = '{:.2f}'.format(a_scaled_low**2 / (G*M_o))
    mass_err_low = '{:.2f}'.format((2*a_scaled_low*a_err_scaled_low/G) / (M_o))

    # Fitted data
    fitted_data = power_law(np.linspace(150,1350,40), param[0], param[1])
    fitted_data = np.hstack((np.vstack((np.linspace(150,1350,40))), np.vstack((fitted_data))))
    
    # Keplerian fit
    keplarian = keplarian_fit(np.linspace(150,1350,40), params_k[0])
    keplarian = np.hstack((np.vstack((np.linspace(150,1350,40))), np.vstack((keplarian))))
    print('Keplerian A: = ' + str(params_k[0]))
    
    # Upper bound
    n_high = param_high[1]
    n_err_high = np.sqrt(np.diag(cov_high))[1]
    a_scaled_high = params_k_high[0] * 10**3 * au**0.5
    a_err_scaled_high = np.sqrt(np.diag(covariance_k_high))[0] * 10**3 * au**0.5
    mass_high = '{:.1f}'.format(a_scaled_high**2 / (G*M_o))
    mass_err_high = '{:.3f}'.format((2*a_scaled_high*a_err_scaled_high/G) / (M_o))
    
    m_low = float(mass_low)
    m_high = float(mass_high)
    print('Lower, Upper MASS: ',m_low,m_high)
    print('Normal bound MASS = '+ str(mass) + ' +/- ' + str(mass_err) + '\n')
    
    # Calculating reduced-chi
    chi_red = chi_calc(array, np.vstack((power_law(array[:,0], param[0], param[1]))))

    if ring == 0:
        plt.title('CH$_{3}$CN Velocity-Distance Curve', size=15)
    else:
        plt.title('CH$_{3}$CN Ring-Averaged Velocity-Distance Curve', size=15)

        
    plt.ylabel('Velocity, km/s', size=13)
    plt.xlabel('Distance, AU', size=13)

    plt.errorbar(array[:,0], array[:,1], yerr = array[:,3]*1, fmt='x',
                markersize=6, alpha=0.5, elinewidth=1.5, capsize=3, label = 'Velocity, $\sigma_{v}$')
    
    plt.plot([],[],label =f'i = {INCLINATION}'+'$^{o}$,'+f' $\phi$ = {POSITION_ANGLE}'+'$^{o}$', alpha=0)

    plt.plot(fitted_data[:,0], fitted_data[:,1], label ='n = ' + 
            str('{:.2f}'.format(n)) +' $\pm$ ' +str('{:.2f}'.format(n_err)))
    
    plt.plot(keplarian[:,0], keplarian[:,1], '--', label = 'Keplerian fit:' )

    plt.plot([],[], label = str('{:.1f}'.format(m_low))+ '$M_{\odot}$ < M < '+str(mass_high) + '$M_{\odot}$', alpha=0)

    plt.plot([],[],label='$\chi^{2}_{red}$ = '+str('{:.2g}'.format(chi_red)))
    plt.legend(loc = 'upper right',borderaxespad=0.5, frameon=True, fontsize=13)
    plt.grid(True)

    plt.show()
    return m_low, mass_high

def n_PA(results_array, position_angle):
    
    results_array = np.array([results_array])
    
    flattened_results = results_array.reshape(-1, results_array.shape[-1])
    
    valid = (flattened_results[:,1] > -1 ) & (flattened_results[:,1] < 1)
    flattened_results = flattened_results[valid]
    
    n_values = flattened_results[:, 1]
    
    PA_values = np.rad2deg(flattened_results[:, 3])

    # Extract the corresponding y-value
    min_epic = np.rad2deg(position_angle)
    y_vals = np.linspace(-0.6,max(flattened_results[:, 1]),2)
    x_vals = np.linspace(min_epic,min_epic,2)

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.scatter(PA_values, n_values, color='blue')
    plt.plot(x_vals,y_vals,'--',color='red', alpha=0.5, label='$\phi$ = '+str(round(min_epic)))

    plt.title('Power Law vs. Position Angle, K=7')
    plt.xlabel('Position Angle (PA)')
    plt.ylabel('Power Law Exponent (n)')
    plt.grid(True)
    plt.legend(loc = 'upper right',borderaxespad=0.5, frameon=True, fontsize=12)
    plt.show()

def x_y_uncertainty(x, y, phi, phi_err):
    """
    x_flat = x.flatten()
    y_flat = y.flatten()
    x_mesh, y_mesh = np.meshgrd(x_flat,y_flat)
    """
    inc = np.deg2rad(INCLINATION)
    xp = x*np.cos(phi) + y*np.sin(phi)#*np.cos(inc)
    yp = -x*np.sin(phi) + y*np.cos(phi)#*np.cos(inc)
    r = np.sqrt(xp**2 + yp**2)
    
    xp_err = np.sqrt( (((x*np.sin(phi)) + (y*np.cos(phi)))*phi_err)**2 ) 
    yp_err = np.sqrt( (((x*np.cos(phi)) + (y*np.sin(phi)))*phi_err)**2 ) 
    
    xp_flat = xp.flatten()
    yp_flat = yp.flatten()
    xp_err_flat = xp_err.flatten()
    yp_err_flat = yp_err.flatten()
    
    xp_err_mesh, yp_err_mesh = np.meshgrid(xp_err_flat,yp_err_flat)
    xp_mesh, yp_mesh = np.meshgrid(xp_flat,yp_flat,indexing='ij')
    
    theta = np.arctan2(yp_mesh, xp_mesh)

    theta = (theta + 2 * np.pi) % (2 * np.pi)
    
    allowed_ranges = [(290, 360), (0, 70), (110, 250)]
    
    
    theta = filter_theta(theta, allowed_ranges)
    
    r = np.sqrt(xp_mesh**2 + yp_mesh**2)
    r_err = np.sqrt( (2*xp_mesh/r * xp_err_mesh)**2 + (2*yp_mesh/r * yp_err_mesh)**2) 
    
    theta_err = np.sqrt( (yp_mesh/r * xp_err)**2 + (xp_mesh/r * yp_err_mesh)**2) 
    
    return r_err, r, theta_err, theta

def velocity_uncertainty(v_obs, v_rad, i, theta, theta_err, i_err):
    
    v_err = np.sqrt( ((V_RAD_ERR/v_rad)*v_obs)**2 + V_RAD_ERR**2 ) 
    
    allowed_ranges = [(290, 360), (0, 70), (110, 250)]
    theta = filter_theta(theta, allowed_ranges)
    
    v_r = np.abs(v_obs / (np.sin(i) * np.cos(theta)))
    v_r_err = np.sqrt( (1/(np.sin(i)**2)*(((v_err**2 + (v_r*np.cos(i)*i_err)**2)) )))
                                #+ ((v_r*np.tan(theta)*theta_err)/(np.cos(theta)))**2 )))
    
    return v_r_err, v_r

def galactic_coords(data, header):
    """
    Using the header file to convert x,y into AU from gal coords.

    """
    
    ctype1 = header['CTYPE1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    crpix1 = header['CRPIX1']
    
    ctype2 = header['CTYPE2']
    crval2 = header['CRVAL2']
    cdelt2 = header['CDELT2']
    crpix2 = header['CRPIX2']
    
    #print(data.shape)
    # Create the coordinate arrays
    x_gal = (crval1 + (np.arange(data.shape[1]) - crpix1) * cdelt1)
    y_gal = (crval2 + (np.arange(data.shape[0]) - crpix2) * cdelt2)
    
    x_au = x_gal * 8.1 * 10**3 * 206264.81 * np.pi/180
    y_au = y_gal * 8.1 * 10**3 * 206264.81 * np.pi/180
    
    return x_au,y_au

def match_shapes_with_nan(data):
    """
    Match the shapes of the 2D array with NaNs to apply the coord trans.

    """
    # Extract the shape of the 2D data array (ignoring the first singleton dimension)
    if data.ndim == 3:
        n, m = data.shape[1], data.shape[2]
    elif data.ndim == 2:
        n, m = data.shape[0], data.shape[1]
    max_dim = max(n, m)

    # Create a new data array filled with NaNs of the maximum dimension
    new_data = np.full((1, max_dim, max_dim), np.nan)
    
    # Determine the padding amounts
    pad_n = (max_dim - n) // 2
    pad_m = (max_dim - m) // 2

    end_n = pad_n + n
    end_m = pad_m + m

    # Ensure that the original data is correctly copied into the padded array
    if data.ndim == 3:
        new_data[0, pad_n:end_n, pad_m:end_m] = data[0, :, :]
    elif data.ndim == 2:
        new_data[0, pad_n:end_n, pad_m:end_m] = data[:, :]
    return new_data

def filter_theta(theta, allowed_ranges):
    """
    Replace values of theta that are outside the allowed ranges with np.nan.

    Parameters:
    - theta: A 2D array of angles in radians.
    - allowed_ranges: A list of tuples, each containing the lower and upper bounds of the allowed ranges in degrees.

    Returns:
    A 2D array with the same shape as theta, where angles outside the allowed ranges have been replaced with np.nan.
    """
    # Convert the ranges from degrees to radians
    allowed_ranges_rad = [(np.deg2rad(low), np.deg2rad(high)) for low, high in allowed_ranges]

    # Create a mask for values within the allowed ranges
    mask = np.zeros_like(theta, dtype=bool)
    for low, high in allowed_ranges_rad:
        mask |= ((theta >= low) & (theta <= high))
    
    # Apply the mask to the theta array, replacing out-of-range values with np.nan
    filtered_theta = np.where(mask, theta, np.nan)
    
    return filtered_theta

def no_nan(data):
    temp = []
    for line in data:
        if not np.any(np.isnan(line)):
            temp.append(line)
    if len(temp) == 0:
        return temp
    temp = np.vstack(temp)
    std = np.std(temp[:,1])
    mean = np.mean(temp[:,1])
    
    temperz = []
    for lineo in temp:
        if lineo[1] < (mean + 5*std): # Removes points 5 standard dev away
            temperz.append(lineo)
            
    return np.vstack((temperz))

def trimming(array):
    """
    To trim the radial data-set to user-defined TRIM
    """
    temp = []
    for line in array:
        if not line[0] > TRIM:
            temp.append(line)
    
    if not len(temp) ==0: # Remove r=0 to avoid undefined velocity
        temp=np.vstack(temp)
    return temp

def keplarian_fit(r,a):
    
    return a * r**-0.5

def fitting_data(name, data, minimum):
    """
    Curve fitting

    """
    x = data[:,0] - minimum
    #print(x,data[:,1])
    params, covariance = curve_fit(name, x, data[:,1], #sigma = data[:,3],
                            absolute_sigma=True, maxfev=5000)
    return params, covariance

def power_law(r, A, n):
    return A * r ** n

def find_closest_n(results_array, target_n=-0.5):
    min_diff = float('inf')
    closest_result = None

    for row in results_array:
        for A_value, n_value, inclination, position_angle, n_uncert in row:
            diff = abs(n_value - target_n)

            if diff < min_diff:
                min_diff = diff
                closest_result = (A_value, n_value, inclination, position_angle, n_uncert)

    return closest_result

def trimmerz(array, trim_1, trim_2):
    """
    Can trim data sets to only return data within two bounds.

    """
    #trim_1 lower bound
    #trim_2 upper bound
    temp = []
    for line in array:
        if not (line[0] < trim_1) | (line[0] > trim_2) | (line[1] < 0.00):
            temp.append(line)
    
    if not len(temp) ==0:
        temp=np.vstack(temp)
    
    
    return temp

def chi_calc(data,fit):
    
    pigs = range(0,len(data))
    chi = 0
    for i in pigs:
        if not data[i,3] == 0: 
            chi += ( (data[i,1]-fit[i,0] )**2 / (data[i,3])**2)
        else:
            chi+=0
    
    chi_red = chi / (len(data[:,1]) - 1)
    
    return chi_red


get_data()
