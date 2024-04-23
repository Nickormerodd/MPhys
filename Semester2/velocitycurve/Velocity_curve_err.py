# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 12:52:59 2023

@author: Christopher
"""

import numpy as np
from astropy.io import fits
from tkinter import Tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from mpl_toolkits.mplot3d import Axes3D as mpl
from matplotlib.lines import Line2D
from scipy.constants import c
#import uncertainties.unumpy as unp
#from uncertainties import ufloat

G = 6.6743e-11
au = 1.495979e+11
M_o = 1.98847e+30
OUT = 'C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/S2/W4-6/Velocity_curve/'

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

    xo_pix,yo_pix = 48,49
    xo_pix,yo_pix = 26,33
    xo_pix,yo_pix = 17,22
    xo_pix,yo_pix = 17,17
    #xo_pix,yo_pix = 24,23
    header = data[0].header
    
    valid_indices = (data[0].data < (np.nanmean(data[0].data - 0.2))) | (data[0].data > (np.nanmean(data[0].data + 0.2)))
    data[0].data[~valid_indices] = np.nan
    #print(np.nanmean(data[0].data))     
    print(data[0].shape)
                            
    data[0].data = data[0].data - np.nanmean(data[0].data)

    amended_data = match_shapes_with_nan(data[0].data)[0]  # This selects the 2D slice.
    #v_obs_flattened = amended_data.flatten()  # Flatten the 2D array to a 1D array

    x_au, y_au = galactic_coords(amended_data, header)
    
    xo_au = x_au[xo_pix] #central xo point
    yo_au = y_au[yo_pix] #central yo point
    
    x_norm = x_au - xo_au #normalising around the xo point
    y_norm = y_au - yo_au #normalising around the yo point
    """
    for i in range(inc.shape[0]):
        x_trans, y_trans = coordinate_transform(x_norm, y_norm, Pos_ang[i])
        #x_trans_8, y_trans_8 = coordinate_transform(x_norm, y_norm, PA_8)

        rv_array = rotational_velocity(amended_data, x_trans, y_trans, inc[i])
        #rv_array_8 = rotational_velocity(amended_data, x_trans_8, y_trans_8, i_8)
    
        main_2(rv_array, name, inc[i], Pos_ang[i])
    """
    #main_3(amended_data, x_norm, y_norm, name)
    main_5(amended_data, x_norm, y_norm, name)
    """
    ##################################### add values here
    position_angle = np.deg2rad(170)
    incatrons = np.deg2rad(38)
    x_trans,y_trans = coord_trans(x_norm, y_norm, position_angle)
    print(x_trans)
    v_rot_matrix, combined= rotational_velocity(amended_data, x_trans, y_trans, np.deg2rad(38), 1)
    # combined is r and v_rot
    print(combined)
    combined = combined[combined[:, 0].argsort()]
    #new_plotting(combined)
    
    """
    ######################################
    return

def main_3(amended_data, x_norm, y_norm, name):
    
    phi = np.deg2rad(156.06060606060606)
    phi_err = np.deg2rad(1)
    i = np.deg2rad(34)
    i_err = np.deg2rad(2)
    
    
    r_err, r, theta_err, theta = x_y_uncertainty(x_norm, y_norm, phi, phi_err)
    v_r_err, v_r = velocity_uncertainty(amended_data, 28.673155,i, theta, theta_err, i_err)
    
    r_flat = r.flatten()
    v_r_flat = v_r.flatten()
    r_err_flat = r_err.flatten()
    v_r_err_flat = v_r_err.flatten()
    
    
    result = np.vstack((r_flat,v_r_flat,r_err_flat,v_r_err_flat)).T
    result = no_nan(result)
    result = trimming(result)
    #result = result[result[:, 0].argsort()]
    #main_4(result,phi,i)
    #print(result)
    return

def main_4(result, phi, inc):
    
    param, cov = fitting_data(power_law, result, 0)
    params_k, covariance_k = fitting_data(keplarian_fit,result, 0)
    
    n = param[1]
    n_err = np.sqrt(np.diag(cov))[1]
    
    #a = params_k[0]
    a_scaled = params_k[0] * 10**3 * au**0.5
    a_err_scaled = np.sqrt(np.diag(covariance_k))[0] * 10**3 * au**0.5
    mass = '{:.3f}'.format(a_scaled**2 / (G*M_o))
    mass_err = '{:.3f}'.format(a_err_scaled**2 / (G*M_o))
    
    i_deg = '{:.1f}'.format(np.rad2deg(inc))
    phi_deg = '{:.1f}'.format(np.rad2deg(phi))
    
    result_2 = trimmerz(result, 250, np.inf)
    result_2 = result_2[result_2[:, 0].argsort()]
    
    param_2, cov_2 = fitting_data(power_law, result_2, 0)
    params_k_2, covariance_k_2 = fitting_data(keplarian_fit,result_2, 0)
    
    #result_3 = trimmerz(result)
    
    plotting(result, param, cov, params_k, covariance_k, phi_deg, i_deg, 
             param_2, cov_2,params_k_2, covariance_k_2)
    
    return

def main_5(data, x_norm,y_norm, name):
    
    incerz = np.linspace(np.deg2rad(34), np.deg2rad(34), 1)
    posizers = np.linspace(np.deg2rad(120), np.deg2rad(190), 500)
    phi_err = np.deg2rad(1.2)
    i_err = np.deg2rad(1.2)
    
    i_mesh, pa_mesh = np.meshgrid(incerz, posizers, indexing='ij')

    # Initialize an empty list to store the results
    results_array = np.full((len(incerz), len(posizers), 5), np.nan)

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
            #print(rv_array)
            #rv_array[:,1] = rv_array[:,1] + rv_array[:,2]
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
    
    A_closest, closest_n, inclination, position_angle, n_uncert = find_closest_n(results_array)
    inca = np.rad2deg(inclination)
    #posa = 157
    posas = np.deg2rad(170)
    
    
    #print(f"Closest n: {closest_n} +/- {n_uncert}, Inclination: {inca}, Position Angle: {posas}\n")
    
    incatrons = np.deg2rad(38)
    """
    
    ##################################### add values here
    position_angle = np.deg2rad(170)
    incatrons = np.deg2rad(38)
    x_trans,y_trans = coord_trans(x_norm, y_norm, position_angle)
    v_rot_matrix, combined= rotational_velocity(data[0].data, x_trans, y_trans, np.deg2rad(38), 1)
    # combined is r and v_rot
    new_plotting(combined)
    
    
    ######################################
    
    
    
    """
    r_err, r, theta_err, theta = x_y_uncertainty(x_norm, y_norm, position_angle, phi_err)
    v_r_err, v_r = velocity_uncertainty(data, 28.9, incatrons, theta, theta_err, i_err)
    
    #v_r = v_r[v_r>25] = np.nan
    
    fig, ax = plt.subplots()
    # Use np.ma.masked_invalid to hide the invalid values in the plot
    cmap = plt.cm.viridis  # Choose a colormap that suits your data
    cmap.set_bad(color='black')  # Set color for masked elements
    im = ax.imshow(np.ma.masked_invalid(v_r),cmap=cmap)
    ax.set_title('Rotational Velocity')
    plt.colorbar(im)  # Show the color bar
    plt.show()
    
    r_flat = r.flatten()
    v_r_flat = v_r.flatten()
    r_err_flat = r_err.flatten()
    v_r_err_flat = v_r_err.flatten()
    
    result = np.vstack((r_flat,v_r_flat,r_err_flat,v_r_err_flat)).T
    result = no_nan(result)
    result = trimming(result)
    result = result[result[:, 0].argsort()]
    
    param, cov = fitting_data(power_law, result, 0)
    params_k, covariance_k = fitting_data(keplarian_fit,result, 0)
    
    result_2 = trimmerz(result, 250, np.inf)
    result_2 = result_2[result_2[:, 0].argsort()]
    
    param_2, cov_2 = fitting_data(power_law, result_2, 0)
    params_k_2, covariance_k_2 = fitting_data(keplarian_fit,result_2, 0)
    
    result_3 = trimmerz(result, 250, np.inf)
    result_3[:,1] = result_3[:,1] + result_3[:,3]
    result_3 = result_3[result_3[:, 0].argsort()]
    param_3, cov_3 = fitting_data(power_law, result_3, 0)
    params_k_3, covariance_k_3 = fitting_data(keplarian_fit,result_3, 0)
    
    averaged_array = ring_averaging(result)
    #averaged_array[:,1] = averaged_array[:,1] + averaged_array[:,3]
    #new_columns = np.zeros((len(averaged_array), 2))  # Create a new array with two columns filled with zeros
    #averaged_array = np.hstack((averaged_array, new_columns))
    
    param_av, cov_av = fitting_data(power_law, trimmerz(averaged_array, 10, np.inf), 0)
    params_k_av, covariance_k_av = fitting_data(keplarian_fit,averaged_array, 0)
    
    m_low, m_high = plotting(result, param, cov, params_k, covariance_k, posas, np.rad2deg(incatrons), param_2, cov_2, params_k_2, covariance_k_2, param_3, cov_3, params_k_3, covariance_k_3, name, 0, 0)
    plotting(averaged_array, param, cov, params_k, covariance_k, posas, np.rad2deg(incatrons), param_av, cov_2, params_k_av, covariance_k_2, param_av, cov_av, params_k_av, covariance_k_av, name, m_low, m_high)
    
    n_PA(results_array, position_angle)
    
    return


def ring_averaging(matrix):
    """
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


def coord_trans(x,y,phi):
    #phi = np.deg2rad(phi)
    inc = np.deg2rad(34)
    alpha = np.pi - phi
    
    x_prime = x*np.cos(alpha) + y*np.sin(alpha)
    y_prime = -x*np.sin(alpha)/np.cos(inc) + y*np.cos(alpha)/np.cos(inc)
    """
    xp = x*np.cos(phi) + y*np.sin(phi)
    yp = -x*np.sin(phi) + y*np.cos(phi)
    """
    return x_prime,y_prime

def plotting(array,param,cov,params_k,covariance_k,phi,inc, param_2, cov_2, params_k_2, covariance_k_2, param_3, cov_3, params_k_3, covariance_k_3, name, m_lower, m_higher):
    
    phi = '{:.0f}'.format(float(phi))
    inc = '{:.0f}'.format(float(inc))
    #print(array)
    
    n = param[1]
    n_err = np.sqrt(np.diag(cov))[1]
    
    if 'k=0' in name.lower():
        k_number = 0

    elif 'k=1' in name.lower():
        k_number = 1

    elif 'k=2' in name.lower():
        k_number = 2

    elif 'k=3' in name.lower():
        k_number = 3

    elif 'k=4' in name.lower():
        k_number = 4

    elif 'k=5' in name.lower():
        k_number = 5

    elif 'k=6' in name.lower():
        k_number = 6

    elif 'k=7' in name.lower():
        k_number = 7

    elif 'k=8' in name.lower():
        k_number = 8

    else:
        title = 'unknown'
        k_number = 'unknown'
    
    #a = params_k[0]
    a_scaled = params_k[0] * 10**3 * au**0.5
    a_err_scaled = np.sqrt(np.diag(covariance_k))[0] * 10**3 * au**0.5
    mass = '{:.2f}'.format(a_scaled**2 / (G*M_o))
    mass_err = '{:.2f}'.format((2*a_scaled*a_err_scaled/G) / (M_o))
    
    plt.figure(figsize=(10, 6))
    
    fitted_data = power_law(array[:,0], param[0], param[1])
    fitted_data = np.hstack((np.vstack((array[:,0])), np.vstack((fitted_data))))
    
    keplarian = keplarian_fit(array[:,0], params_k[0])
    keplarian = np.hstack((np.vstack((array[:,0])), np.vstack((keplarian))))
    
    n_2 = param_2[1]
    n_err_2 = np.sqrt(np.diag(cov_2))[1]
    
    #a = params_k[0]
    a_scaled_2 = params_k_2[0] * 10**3 * au**0.5
    a_err_scaled_2 = np.sqrt(np.diag(covariance_k_2))[0] * 10**3 * au**0.5
    mass_2 = '{:.2f}'.format(a_scaled_2**2 / (G*M_o))
    mass_err_2 = '{:.2f}'.format((2*a_scaled_2*a_err_scaled_2/G) / (M_o))

    #fitted_data_2 = power_law(trimmerz(array, 100, np.inf)[:,0], param_2[0], param_2[1])
    fitted_data_2 = power_law(np.linspace(150,1350,40), param_2[0], param_2[1])
    #fitted_data_2 = np.hstack((np.vstack((trimmerz(array, 100, np.inf)[:,0])), np.vstack((fitted_data_2))))
    fitted_data_2 = np.hstack((np.vstack((np.linspace(150,1350,40))), np.vstack((fitted_data_2))))
    
    #keplarian_2 = keplarian_fit(trimmerz(array, 100, np.inf)[:,0], params_k_2[0])
    keplarian_2 = keplarian_fit(np.linspace(150,1350,40), params_k_2[0])
    #keplarian_2 = np.hstack((np.vstack((trimmerz(array, 100, np.inf)[:,0])), np.vstack((keplarian_2))))
    keplarian_2 = np.hstack((np.vstack((np.linspace(150,1350,40))), np.vstack((keplarian_2))))
    
    n_3 = param_3[1]
    n_err_3 = np.sqrt(np.diag(cov_3))[1]
    a_scaled_3 = params_k_3[0] * 10**3 * au**0.5
    a_err_scaled_3 = np.sqrt(np.diag(covariance_k_3))[0] * 10**3 * au**0.5
    mass_3 = '{:.1f}'.format(a_scaled_3**2 / (G*M_o))
    mass_err_3 = '{:.3f}'.format((2*a_scaled_3*a_err_scaled_3/G) / (M_o))
    m_low = float(mass_2) - float(mass_err_2)
    m_high = float(mass_3) + float(mass_err_3)
    print(m_low,m_high)
    #m_low_f = '{:.3f}'.format(m_low)

    fitted_data_3= power_law(trimmerz(array, 100, np.inf)[:,0], param_3[0], param_3[1])
    fitted_data_3 = np.hstack((np.vstack((trimmerz(array, 100, np.inf)[:,0])), np.vstack((fitted_data_3))))
    
    keplarian_3 = keplarian_fit(trimmerz(array, 200, np.inf)[:,0], params_k_3[0])
    keplarian_3 = np.hstack((np.vstack((trimmerz(array, 200, np.inf)[:,0])), np.vstack((keplarian_3))))
    print('Keplerian A: = ' + str(params_k_3[0]))
    
    title = 'CH3CN_K=' + str(k_number) +' Velocity_Curve.pdf'
    output_path = OUT + title
    
    if not 'v_off' in name.lower():
        if m_lower == 0:
            plt.title('CH$_{3}$CN: K=' + str(k_number) +' Velocity-Distance Curve')
        else:
            plt.title('CH$_{3}$CN: K=' + str(k_number) +' Ring-Averaged Velocity-Distance Curve')
    else:
        if m_lower == 0:
            plt.title('CH$_{3}$CN Velocity-Distance Curve')
        else:
            plt.title('CH$_{3}^{\:13}$CN Ring-Averaged Velocity-Distance Curve')
            #np.savetxt('C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/S2/W4-6/Velocity_curve/CH3CN_ring-av.txt', array)
        
    plt.ylabel('Velocity, km/s')
    plt.xlabel('Distance, AU')
    
    chi_red = chi_calc(array,fitted_data)

    plt.errorbar(array[:,0], array[:,1]*1, yerr = array[:,3]*1/30, fmt='x',
                markersize=6, alpha=0.5, elinewidth=1.5, capsize=3, label = 'Velocity Dispersion, $\sigma_{v}$')
    
    plt.plot([],[],label =f'i = {inc}'+'$^{o}$,'+f' $\phi$ = 170'+'$^{o}$', alpha=0)
    """
    plt.plot(fitted_data[:,0], fitted_data[:,1], label ='n = ' + 
            str('{:.2f}'.format(n)) +'$\pm$' +str('{:.2f}'.format(n_err)))
    
    plt.plot(keplarian[:,0], keplarian[:,1], label = 'Keplarian fit:\nM = ' +
             str(mass) + ' $\pm$ '+str(mass_err) + ' $M_{\odot}$', alpha=0)
    """
    #plt.plot([],[],alpha=0, label=' ')
    plt.plot(fitted_data_2[:,0], fitted_data_2[:,1]*1, label ='n = ' + 
            str('{:.2f}'.format(n_2)) +' $\pm$ ' +str('{:.2f}'.format(n_err_2*2.2)))
    
    plt.plot(keplarian_2[:,0], keplarian_2[:,1]*1, '--', label = 'Keplerian fit:' )#+
           #  str(mass_2) + ' $\pm$ '+str(mass_err_2)+ ' $M_{\odot}$', alpha=0)
    #plt.plot(fitted_data_3[:,0], fitted_data_3[:,1], label ='n = ' + 
    #       str('{:.3f}'.format(n_3)) +'$\pm$' +str('{:.2f}'.format(n_err_3)))
    
    #plt.plot(keplarian_3[:,0], keplarian_3[:,1], label = 'Keplarian fit\nM = ' +
    #         str(mass_3) + ' $\pm$ '+str(mass_err_3)+ ' $M_{\odot}$')
    #print('Angular Velocity = ' + str(keplarian[:,1]/keplarian[:,0]))
    if m_lower == 0:
        plt.plot([],[], label = str('{:.1f}'.format(m_low))+ '$M_{\odot}$ < M < '+str(mass_3) + '$M_{\odot}$', alpha=0)
    else:
        plt.plot([],[], label = str('{:.1f}'.format(m_lower))+ '$M_{\odot}$ < M < '+str(m_higher) + '$M_{\odot}$', alpha=0)
        #plt.plot([],[], label = str('5.6$M_{\odot}$ < M < 14.9$M_{\odot}$'), alpha=0)
        
    plt.plot([],[],label='$\chi^{2}_{red}$ = ' + str('{:.2f}'.format(chi_red)))
    plt.legend(loc = 'upper right',borderaxespad=0.5, frameon=True, fontsize=12)
    plt.grid(True)
    #plt.savefig(output_path,bbox_inches='tight')
    #np.savetxt('Keplarian_data_k=7.csv',keplarian,delimiter='')
    if not m_lower == 0:
        #plt.savefig('Velocity-Distance_RingAv_XCLASS.pdf',bbox_inches='tight')
        title2 = 'CH3CN_K=' + str(k_number) +'_RingAv_Velocity_Curve.pdf'
        plt.savefig(title2,bbox_inches='tight')
        
    plt.show()
    return m_low, mass_3

def n_PA(results_array, position_angle):
    
    results_array = np.array([results_array])
    
    flattened_results = results_array.reshape(-1, results_array.shape[-1])
    #print(flattened_results)
    
    valid = (flattened_results[:,1] > -1 ) & (flattened_results[:,1] < 1)
    flattened_results = flattened_results[valid]
    
    n_values = flattened_results[:, 1]
    
    PA_values = np.rad2deg(flattened_results[:, 3])
    #min_x_index = np.argmin(flattened_results[:, 0])

    # Extract the corresponding y-value
    #min_PA = flattened_results[min_x_index, 3]
    min_epic = np.rad2deg(position_angle)
    y_vals = np.linspace(-0.6,max(flattened_results[:, 1]),2)
    x_vals = np.linspace(min_epic,min_epic,2)
   # minerz_pa = flattened_results[min]
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.scatter(PA_values, n_values, color='blue')
    plt.plot(x_vals,y_vals,'--',color='red', alpha=0.5, label='$\phi$ = '+str(round(min_epic)))
    #plt.plot(PA_values, n_values)
    plt.title('Power Law vs. Position Angle, K=7')
    plt.xlabel('Position Angle (PA)')
    plt.ylabel('Power Law Exponent (n)')
    plt.grid(True)
    plt.legend(loc = 'upper right',borderaxespad=0.5, frameon=True, fontsize=12)
    #plt.savefig('Power_law_vs_PA_k=7.png', bbox_inches='tight',dpi=500)
    plt.show()
    #plt.savefig('Power_law_vs_PA_k=7.png', bbox_inches='tight',dpi=500)

def x_y_uncertainty(x, y, phi, phi_err):
    """
    x_flat = x.flatten()
    y_flat = y.flatten()
    x_mesh, y_mesh = np.meshgrd(x_flat,y_flat)
    """
    inc = np.deg2rad(34)
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
    ########################
    theta = (theta + 2 * np.pi) % (2 * np.pi)
    
    allowed_ranges = [(290, 360), (0, 70), (110, 250)]
    
    
    theta = filter_theta(theta, allowed_ranges)
    
    r = np.sqrt(xp_mesh**2 + yp_mesh**2)
    r_err = np.sqrt( (2*xp_mesh/r * xp_err_mesh)**2 + (2*yp_mesh/r * yp_err_mesh)**2) 
    
    theta_err = np.sqrt( (yp_mesh/r * xp_err)**2 + (xp_mesh/r * yp_err_mesh)**2) 
    
    return r_err, r, theta_err, theta

def velocity_uncertainty(v_obs, v_rad, i, theta, theta_err, i_err):
    """
    ################
    allowed_ranges = [(290, 360), (0, 70), (110, 250)]
    
    mask = np.zeros_like(theta, dtype=bool)
    
    # Set mask to True for angles within the allowed ranges
    for low, high in allowed_ranges:
        mask |= (theta >= low) & (theta <= high)
    
    v_obs = np.where(mask, v_obs, np.nan)
    
    # Calculate the rotational velocity, avoiding division by zero
    sin_inclination = np.sin(i)
    cos_theta = np.cos(theta)
    
    # Mask out the values where cos(theta) is zero (semi-minor axis)
    cos_theta_masked = np.where(mask, cos_theta, np.nan)
    
    # Compute the rotational velocity
    v_rot = v_obs / (sin_inclination * cos_theta_masked)
    
    # Replace nan and inf with zero after the calculation
    #v_rot = np.nan_to_num(v_rot, nan=0.0, posinf=0.0, neginf=0.0)
    
    ##################
    """
    
    #print(np.nanmax(theta))
    #v_err = np.sqrt( ((0.236/v_rad)*v_obs)**2 + 0.236**2 )
    v = (v_obs - v_rad)/(np.sin(i)*np.cos(theta))
    v_err = np.sqrt( ((0.0707770409/v_rad)*v_obs)**2 + 0.0707770409**2 ) 
    #v = (v_obs - v_rad)/(np.sin(i)*np.cos(theta))
    
    v_r = np.abs(v_obs / (np.sin(i) * np.cos(theta)))
    v_r_err = np.sqrt( (1/(np.sin(i)**2)*(((v_err**2 + (v*np.cos(i)*i_err)**2)) ))) #+ ((v*np.tan(theta)*theta_err)/(np.cos(theta)))**2))) 
    
    
    #v_r = np.abs(v_obs / (np.sin(i) * np.cos(theta)))

    
    return v_r_err, v_r

def galactic_coords(data, header):
    
    #header = data[0].header
    
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
    #print('Std = ' + str(std))
    #print('Mean = ' + str(mean)) 
    
    temperz = []
    for lineo in temp:
        if lineo[1] < (mean + 7*std):
            temperz.append(lineo)
            
    return np.vstack((temperz))

def trimming(array):
    temp = []
    for line in array:
        if not line[0] > 1750:
            temp.append(line)
    
    if not len(temp) ==0:
        temp=np.vstack(temp)
    return temp

def keplarian_fit(r,a):
    
    return a * r**-0.5

def fitting_data(name, data, minimum):
    
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
            chi += ( (data[i,1]-fit[i,1] )**2 / (data[i,3] / 18)**2)
        else:
            chi+=0
    
    chi_red = chi / (len(data[:,1]) - 1)
    
    return chi_red

def rotational_velocity_real(v_obs, x_trans, y_trans, inclination_deg=34, A_kep_SI=0):
    # Generate meshgrid from the flattened transposed coordinates
    x_mesh, y_mesh = np.meshgrid(x_trans, y_trans)
    
    # Calculate the radius r from the meshgrid
    r = np.sqrt(x_mesh**2 + y_mesh**2)
    
    # Scale the observed velocity
    v_obs_scaled = v_obs# - 29  # Assuming 29 is the systemic velocity to be subtracted
    
    # Calculate theta in radians from the meshgrid and convert to degrees
    theta = np.rad2deg(np.arctan2(y_mesh, x_mesh))
    
    # Normalize theta to be within the range [0, 360)
    theta = (theta + 360) % 360
    
    # Define the allowed ranges in degrees, excluding the semi-minor axis
    allowed_ranges = [(0, 70), (115, 250), (295, 360)]
    
    # Initialize mask to include all angles
    mask = np.zeros_like(theta, dtype=bool)
    
    # Set mask to True for angles within the allowed ranges
    for low, high in allowed_ranges:
        mask |= (theta >= low) & (theta <= high)
    
    # Apply the mask to v_obs, setting out-of-range values to np.nan
    v_obs_masked = np.where(mask, v_obs_scaled, np.nan)
    
    # Calculate the Keplerian velocity for regions where the velocity is zero or masked
    v_keplerian = A_kep_SI / np.sqrt(r)
    v_obs_masked = np.where((v_obs_masked == 0) | np.isnan(v_obs_masked), v_keplerian, v_obs_masked)
    
    # Calculate the rotational velocity, avoiding division by zero
    sin_inclination = np.sin(np.deg2rad(inclination_deg))
    cos_theta = np.cos(np.deg2rad(theta))
    
    # Mask out the values where cos(theta) is zero (semi-minor axis)
    cos_theta_masked = np.where(mask, cos_theta, np.nan)
    
    # Compute the rotational velocity
    v_rot = v_obs_masked / (sin_inclination * cos_theta_masked)
    
    # Replace nan and inf with zero after the calculation
    v_rot = np.nan_to_num(v_rot, nan=0.0, posinf=0.0, neginf=0.0)
    """
    # Plot the rotational velocity field
    fig, ax = plt.subplots()
    im = ax.imshow(v_rot, extent=[x_mesh.min(), x_mesh.max(), y_mesh.min(), y_mesh.max()], cmap='viridis')
    ax.set_title('Rotational Velocity')
    plt.colorbar(im)
    plt.show()
    """
    return v_rot

def rotational_velocity(v_obs_scaled, x_trans, y_trans, inclination_deg, numba):
    # Generate meshgrid from the flattened transposed coordinates
    x_mesh, y_mesh = np.meshgrid(x_trans, y_trans)
    
    # Calculate the radius r from the meshgrid
    r = np.sqrt(x_mesh**2 + y_mesh**2)
    
    # Scale the observed velocity
     # Assuming 29 is the systemic velocity to be subtracted
    
    # Calculate theta in radians from the meshgrid and convert to degrees
    theta = np.rad2deg(np.arctan2(y_mesh, x_mesh))
    
    # Normalize theta to be within the range [0, 360)
    theta = (theta + 360) % 360
    
    # Define the allowed ranges in degrees, excluding the semi-minor axis
    allowed_ranges = [(0, 75), (110, 255), (285, 360)]
    
    # Initialize mask to include all angles
    mask = np.zeros_like(theta, dtype=bool)
    
    # Set mask to True for angles within the allowed ranges
    for low, high in allowed_ranges:
        mask |= (theta >= low) & (theta <= high)
    
    # Apply the mask to v_obs, setting out-of-range values to np.nan
    v_obs_masked = np.where(mask, v_obs_scaled, np.nan)
    
    # Calculate the Keplerian velocity for regions where the velocity is zero or masked
    #v_keplerian = A_kep_SI / np.sqrt(r)
    #v_obs_masked = np.where((v_obs_masked == 0) | np.isnan(v_obs_masked), v_keplerian, v_obs_masked)
    
    # Calculate the rotational velocity, avoiding division by zero
    sin_inclination = np.sin((inclination_deg))
    cos_theta = np.cos(np.deg2rad(theta))
    
    # Mask out the values where cos(theta) is zero (semi-minor axis)
    cos_theta_masked = np.where(mask, cos_theta, np.nan)
    
    # Compute the rotational velocity
    v_rot = v_obs_masked / (sin_inclination * cos_theta_masked)
    
    # Replace nan and inf with zero after the calculation
    v_rot_matrix = np.nan_to_num(v_rot, nan=0.0, posinf=0.0, neginf=0.0)
    
    if numba == 1:
        # Plot the rotational velocity field
        fig, ax = plt.subplots()
        im = ax.imshow(v_rot, extent=[x_mesh.min(), x_mesh.max(), y_mesh.min(), y_mesh.max()], cmap='viridis')
        ax.set_title('Rotational Velocity')
        plt.colorbar(im)
        plt.show()
    #print(np.nanmax(v_rot))
    combined =  np.vstack((r.flatten(),np.abs(v_rot_matrix.flatten()))).T
    
    return v_rot_matrix, combined

get_data()
