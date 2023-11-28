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
    header = data[0].header
    
    valid_indices = (data[0].data < (np.nanmean(data[0].data - 0.4))) | (data[0].data > (np.nanmean(data[0].data + 0.4)))
    data[0].data[~valid_indices] = np.nan
    #print(np.nanmean(data[0].data))     
                              
                            
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
    main_3(amended_data, x_norm, y_norm, name)
    main_5(amended_data, x_norm, y_norm)
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

def main_5(data, x_norm,y_norm):
    
    incerz = np.linspace(np.deg2rad(32), np.deg2rad(35), 5)
    posizers = np.linspace(np.deg2rad(120), np.deg2rad(200), 400)
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
            r_flat = r.flatten()
            v_r_flat = v_r.flatten()
            v_r_err_flat = v_r_err.flatten()
            rv_array = np.vstack((r_flat,v_r_flat,v_r_err_flat)).T
            
            rv_array = no_nan(rv_array)
            if len(rv_array) == 0:
                continue
            #print(rv_array)
            #rv_array[:,1] = rv_array[:,1] + rv_array[:,2]
            new_array = trimmerz(trimming(rv_array), 250, np.inf)
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
    posas = np.rad2deg(position_angle)
    print(f"Closest n: {closest_n} +/- {n_uncert}, Inclination: {inca}, Position Angle: {posas}\n")
    
    incatrons = np.deg2rad(34)
    
    r_err, r, theta_err, theta = x_y_uncertainty(x_norm, y_norm, position_angle, phi_err)
    v_r_err, v_r = velocity_uncertainty(data, 28.9, incatrons, theta, theta_err, i_err)
    
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
    
    plotting(result, param, cov, params_k, covariance_k, posas, np.rad2deg(incatrons), param_2, cov_2, params_k_2, covariance_k_2, param_3, cov_3, params_k_3, covariance_k_3)
    
    n_PA(results_array, position_angle)
    return

def plotting(array,param,cov,params_k,covariance_k,phi,inc, param_2, cov_2, params_k_2, covariance_k_2, param_3, cov_3, params_k_3, covariance_k_3):
    
    phi = '{:.0f}'.format(float(phi))
    inc = '{:.0f}'.format(float(inc))
    
    n = param[1]
    n_err = np.sqrt(np.diag(cov))[1]
    
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

    fitted_data_2 = power_law(trimmerz(array, 200, np.inf)[:,0], param_2[0], param_2[1])
    fitted_data_2 = np.hstack((np.vstack((trimmerz(array, 200, np.inf)[:,0])), np.vstack((fitted_data_2))))
    
    keplarian_2 = keplarian_fit(trimmerz(array, 200, np.inf)[:,0], params_k_2[0])
    keplarian_2 = np.hstack((np.vstack((trimmerz(array, 200, np.inf)[:,0])), np.vstack((keplarian_2))))
    
    n_3 = param_3[1]
    n_err_3 = np.sqrt(np.diag(cov_3))[1]
    a_scaled_3 = params_k_3[0] * 10**3 * au**0.5
    a_err_scaled_3 = np.sqrt(np.diag(covariance_k_3))[0] * 10**3 * au**0.5
    mass_3 = '{:.1f}'.format(a_scaled_3**2 / (G*M_o))
    mass_err_3 = '{:.3f}'.format((2*a_scaled_3*a_err_scaled_3/G) / (M_o))
    m_low = float(mass_2) - float(mass_err_2)
    #m_low_f = '{:.3f}'.format(m_low)

    fitted_data_3= power_law(trimmerz(array, 200, np.inf)[:,0], param_3[0], param_3[1])
    fitted_data_3 = np.hstack((np.vstack((trimmerz(array, 200, np.inf)[:,0])), np.vstack((fitted_data_3))))
    
    keplarian_3 = keplarian_fit(trimmerz(array, 200, np.inf)[:,0], params_k_3[0])
    keplarian_3 = np.hstack((np.vstack((trimmerz(array, 200, np.inf)[:,0])), np.vstack((keplarian_3))))
    
    plt.title('CH$_{3}$CN: K=7 Velocity-Distance Curve')
    plt.ylabel('Velocity, km/s')
    plt.xlabel('Distance, AU')

    plt.errorbar(array[:,0], array[:,1], yerr = array[:,3]*0, xerr = array[:,2]*0, fmt='x',
                markersize=6, alpha=0.5, elinewidth=1.5)
    
    plt.plot([],[],label =f'i = {inc}'+'$^{o}$,'+f' $\phi$ = {phi}'+'$^{o}$', alpha=0)
    """
    plt.plot(fitted_data[:,0], fitted_data[:,1], label ='n = ' + 
            str('{:.2f}'.format(n)) +'$\pm$' +str('{:.2f}'.format(n_err)))
    
    plt.plot(keplarian[:,0], keplarian[:,1], label = 'Keplarian fit:\nM = ' +
             str(mass) + ' $\pm$ '+str(mass_err) + ' $M_{\odot}$', alpha=0)
    """
    #plt.plot([],[],alpha=0, label=' ')
    plt.plot(fitted_data_2[:,0], fitted_data_2[:,1], label ='n = ' + 
            str('{:.2f}'.format(n_2)) +' $\pm$ ' +str('{:.2f}'.format(n_err_2)))
    
    plt.plot(keplarian_2[:,0], keplarian_2[:,1], '--', label = 'Keplarian fit:' )#+
           #  str(mass_2) + ' $\pm$ '+str(mass_err_2)+ ' $M_{\odot}$', alpha=0)
    #plt.plot(fitted_data_3[:,0], fitted_data_3[:,1], label ='n = ' + 
    #       str('{:.3f}'.format(n_3)) +'$\pm$' +str('{:.2f}'.format(n_err_3)))
    
    #plt.plot(keplarian_3[:,0], keplarian_3[:,1], label = 'Keplarian fit\nM = ' +
    #         str(mass_3) + ' $\pm$ '+str(mass_err_3)+ ' $M_{\odot}$')
    plt.plot([],[], label = str('{:.1f}'.format(m_low))+ '$M_{\odot}$ < M < '+str(mass_3) + '$M_{\odot}$', alpha=0)
    plt.legend(loc = 'upper right',borderaxespad=0.5, frameon=True, fontsize=12)
    plt.grid(True)
    plt.savefig('Velocity-Distance_curve_K=7.png',bbox_inches='tight',dpi=500)
    plt.show()
    #plt.savefig('Velocity-Distance_curve_K=7.png',bbox_inches='tight',dpi=500)
    return

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
    plt.savefig('Power_law_vs_PA_k=7.png', bbox_inches='tight',dpi=500)
    plt.show()
    #plt.savefig('Power_law_vs_PA_k=7.png', bbox_inches='tight',dpi=500)

def x_y_uncertainty(x, y, phi, phi_err):
    """
    x_flat = x.flatten()
    y_flat = y.flatten()
    x_mesh, y_mesh = np.meshgrd(x_flat,y_flat)
    """
    xp = x*np.cos(phi) + y*np.sin(phi)
    yp = -x*np.sin(phi) + y*np.cos(phi)
    r = np.sqrt(xp**2 + yp**2)
    
    xp_err = np.sqrt( (((x*np.sin(phi)) + (y*np.cos(phi)))*phi_err)**2 ) 
    yp_err = np.sqrt( (((x*np.cos(phi)) + (y*np.sin(phi)))*phi_err)**2 ) 
    
    xp_flat = xp.flatten()
    yp_flat = yp.flatten()
    xp_err_flat = xp_err.flatten()
    yp_err_flat = yp_err.flatten()
    
    xp_err_mesh, yp_err_mesh = np.meshgrid(xp_err_flat,yp_err_flat,indexing='ij')
    xp_mesh, yp_mesh = np.meshgrid(xp_flat,yp_flat,indexing='ij')
    
    theta = np.arctan2(yp_mesh, xp_mesh)
    allowed_ranges = [(290, 360), (0, 70), (110, 250)]
    theta = filter_theta(theta, allowed_ranges)
    
    r = np.sqrt(xp_mesh**2 + yp_mesh**2)
    r_err = np.sqrt( (2*xp_mesh/r * xp_err_mesh)**2 + (2*yp_mesh/r * yp_err_mesh)**2) 
    
    theta_err = np.sqrt( (yp_mesh/r * xp_err)**2 + (xp_mesh/r * yp_err_mesh)**2) 
    
    return r_err, r, theta_err, theta

def velocity_uncertainty(v_obs, v_rad,i, theta, theta_err, i_err):
    
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
    n, m = data.shape[1], data.shape[2]
    max_dim = max(n, m)

    # Create a new data array filled with NaNs of the maximum dimension
    new_data = np.full((1, max_dim, max_dim), np.nan)
    
    # Determine the padding amounts
    pad_n = (max_dim - n) // 2
    pad_m = (max_dim - m) // 2

    end_n = pad_n + n
    end_m = pad_m + m

    # Ensure that the original data is correctly copied into the padded array
    new_data[0, pad_n:end_n, pad_m:end_m] = data[0, :, :]
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
        if lineo[1] < (mean + 10*std):
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
    params, covariance = curve_fit(name, x, data[:,1], sigma = data[:,3],
                            absolute_sigma=True, maxfev=3000)
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
        if not (line[0] < trim_1) | (line[0] > trim_2):
            temp.append(line)
    
    if not len(temp) ==0:
        temp=np.vstack(temp)
    
    
    return temp

get_data()

