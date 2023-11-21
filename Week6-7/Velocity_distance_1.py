"""
Created on Thu Nov 16 16:06:58 2023

@author: Christopher
"""

"""
Created on Thu Nov  9 18:45:26 2023

@author: Christopher

Input a moment 1 map fits files and it will produce a rot velcity-distance
plot, fit a keplarian curve to the data, and extract the mass of the
protostar. It will fit two curves, one being to the whole data (purple),
and the other being to the data at which it peaks, slightly trimmed data 
(orange curve), and the mass is estimated. Mass from orange curve is more 
uncertain than that from purple curve. However, orange curve will hopefullly
fit to the data more accurately, representing keplarian motion of the disc
to a better degree.

You can ammend the inclination angle and position angles at the top, and add 
more values to it. The code iterates through them. 

#########
This code has main_2 hashtagged out
########


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

#ina = 30
#PA = 150

#inc = np.linspace(np.deg2rad(ina-15),np.deg2rad(ina+15),7)
#Pos_ang = np.linspace(np.deg2rad(PA-30),np.deg2rad(PA+30),7)

inc = np.array([np.deg2rad(35), np.deg2rad(89), np.deg2rad(40), np.deg2rad(30)])#, np.deg2rad(28)])
Pos_ang = np.array([np.deg2rad(165), np.deg2rad(160), np.deg2rad(160), np.deg2rad(156)])#, np.deg2rad(160)])

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
    
    valid_indices = (data[0].data < (np.nanmean(data[0].data - 0.2))) | (data[0].data > (np.nanmean(data[0].data + 0.2)))
    data[0].data[~valid_indices] = np.nan
    #print(np.max(data[0].data))     
                              
                            
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
    return

def main_2(rv_array, name, i, PA):
    
    rv_array = no_nan(rv_array)
    
    rv_array = rv_array[rv_array[:, 0].argsort()]
    #max_y_value = find_max_value_coordinates(no_nan(rv_array))
    #filtered_data_2 = filter_data(rv_array, rv_array[max_y_value,0], 4000)
    
    filtered_data = filter_data(rv_array, 200, 3000)
    params, cov = fitting_data(filtered_data,0)
    fitted_y_data = keplarian_fit(filtered_data[:,0], params[0])

    a_scaled = params[0] * 10**3 * au**0.5
    mass = a_scaled**2 / (G*M_o)
    fitted_data = np.hstack((np.vstack((filtered_data[:,0])), np.vstack((fitted_y_data))))
    #fitted_data = fitted_data[fitted_data[:, 0].argsort()]
    
    max_y_value = find_max_value_coordinates(no_nan(rv_array)) 
    
    filter_temp = filter_data(rv_array, 200, 1750)
    filter_temp = filter_data_2(filter_temp, np.nanmean(filter_temp[:,1]), 1000)
    weighted_mean_x = weighted_mean_calc(filter_temp)
    value = decision(weighted_mean_x,rv_array[max_y_value,0])
    #print(max_y_value)
    #filtered_data_2 = filter_data(rv_array, rv_array[max_y_value,0], 4000)
    filtered_data_2 = filter_data(rv_array, value, 3000)
    #filtered_data_3 = filter_data(rv_array, rv_array[max_y_value,0], 6000)
    
    #filtered_data_3 = filter_data_2(filtered_data_2, params[0])
    minimum =  np.nanmin(filtered_data_2[:,0]) - 70
    
    #print(minimum)
    params_2, cov_2 = fitting_data(filtered_data_2, minimum)
    fitted_y_data_2 = keplarian_fit(filtered_data_2[:,0] - minimum, params_2[0]) #+ minimum
    a_scaled_2 = params_2[0] * 10**3 * au**0.5
    
    #filtered_data_3 = filter_data_2(filtered_data_2, params_2[0])
    
    mass_2 = a_scaled_2**2 / (G*M_o)
    fitted_data_2 = np.hstack((np.vstack((filtered_data_2[:,0])), np.vstack((fitted_y_data_2))))
    """
    params_3, cov_3 = fitting_data(filtered_data_3, minimum)
    fitted_y_data_3 = keplarian_fit(filtered_data_3[:,0] - minimum, params_3[0]) #+ minimum
    a_scaled_3 = params_3[0] * 10**3 * au**0.5
    
    fitted_data_3 = np.hstack((np.vstack((filtered_data_3[:,0])), np.vstack((fitted_y_data_3))))
    """
    masserz = mass_(fitted_data_2)
    #masserz_3 = mass_(fitted_data_3)
    
    params_3, cov_3 = fitting_data(fitted_data_2, 0)
    a_scaled_3 = params_3[0] * 10**3 * au**0.5
    
    #filtered_data_3 = filter_data_2(filtered_data_2, params_2[0])
    
    fitted_y_data_3 = keplarian_fit(filtered_data_2[:,0], params_3[0])
    fitted_data_3 = np.hstack((np.vstack((filtered_data_2[:,0])), np.vstack((fitted_y_data_3))))
    
    mass_3 = a_scaled_3**2 / (G*M_o)
    #print(mass_3)
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.scatter(rv_array[:,0], rv_array[:,1],s=3)
    ax.plot(fitted_data_3[:,0], fitted_data_3[:,1], color='green')
    ax.plot(fitted_data_2[:,0], fitted_data_2[:,1], color='pink')
    """
    #plotting_2(rv_array, i, PA, fitted_data, mass, name, fitted_data_2,
    #         masserz)
    
    # Convert data to a NumPy array if it's not already
    
    #main_3(rv_array)
    
    new_data = trimming(rv_array)

    # Split the data into velocity and distance
    velocity = new_data[:, 1]
    distance = new_data[:, 0]

    # Perform the curve fitting
    best_params , _ = curve_fit(power_law, distance, velocity)

    # Extract A and n
    A_4, n_4 = best_params
    n_perr = np.sqrt(np.diag(_))[1]
    
    inc_deg = np.rad2deg(i)
    PA_deg = np.rad2deg(PA)
    
    print(f'Name = {name} \nPower = {n_4} +/- {n_perr}\nfor {inc_deg},{PA_deg}\n')
    
    return

def main_3(data, x_norm, y_norm, name):
    
    # Define the ranges for inclination and position angles
    incerz = np.linspace(np.deg2rad(32), np.deg2rad(40), 40)
    posizers = np.linspace(np.deg2rad(150), np.deg2rad(160), 100)
    
    # Create a meshgrid for iteration
    i_mesh, pa_mesh = np.meshgrid(incerz, posizers, indexing='ij')

    # Initialize an empty list to store the results
    results_array = np.full((len(incerz), len(posizers), 5), np.nan)

    # Iterate through the meshgrid
    for i in range(len(incerz)):
        for j in range(len(posizers)):
            inclination = i_mesh[i, j]
            position_angle = pa_mesh[i, j]

            # Transform coordinates
            x_trans, y_trans = coordinate_transform(x_norm, y_norm, position_angle)

            # Calculate rotational velocity
            rv_array = rotational_velocity(data, x_trans, y_trans, inclination)
            #print(rv_array)
            rv_array = no_nan(rv_array)
            if len(rv_array) == 0:
                continue
            #print(rv_array)
            new_array = trimming(rv_array)
            if len(new_array) == 0:
                continue
            #print(new_array)
            # Fit the power law and extract parameters
            try:
                best_params, covariance = curve_fit(power_law, new_array[:,0], new_array[:,1],
                                                    maxfev=5000)
                
                A_4, n_4 = best_params
                n_perr = np.sqrt(np.diag(covariance))[1]

            except RuntimeError:
                A_4, n_4 = np.nan, np.nan


            # Append the results
            results_array[i, j] = [A_4, n_4, inclination, position_angle, n_perr]

    # Convert the results to a 2D numpy array
    #results_array = np.array(results).reshape(len(incerz), len(posizers), 3)
    #print(results_array)
    A_closest, closest_n, inclination, position_angle, n_uncert = find_closest_n(results_array)
    inca = np.rad2deg(inclination)
    posa = np.rad2deg(position_angle)
    
    print(f"Closest n: {closest_n} +/- {n_uncert}, Inclination: {inca}, Position Angle: {posa}\n")
    
    x_tranerz, y_tranerz = coordinate_transform(x_norm, y_norm, position_angle)

            # Calculate rotational velocity

    final_array = rotational_velocity(data, x_tranerz, y_tranerz, inclination)
    final_array = trimming(final_array)
    model_plot(final_array, float(closest_n), float(inclination), float(position_angle),
               float(A_closest), name, float(n_uncert))
    plot_n_contours(results_array, incerz, posizers, name)
    n_PA(results_array, name)
    #print(results_array)
        
    return

def model_plot(array, n, i, PA, A, name, n_uncert):
    plt.figure(figsize=(10, 6))
    
    i_deg = '{:.1f}'.format(np.rad2deg(i))
    PA_deg = '{:.1f}'.format(np.rad2deg(PA))
    
    array = no_nan(array)
    array = array[array[:, 0].argsort()]
    
    fitted_data = power_law(array[:,0], A, n)
    fitted_data = np.hstack((np.vstack((array[:,0])), np.vstack((fitted_data))))
    
    param_kep, _ = fitting_data(array, 0)
    keplarian = keplarian_fit(array[:,0], param_kep[0])
    keplarian = np.hstack((np.vstack((array[:,0])), np.vstack((keplarian))))
    
    a_scaled = param_kep[0] * 10**3 * au**0.5
    mass = '{:.3f}'.format(a_scaled**2 / (G*M_o))
    
    plt.title(f'Model Velocity-Distance for {name}')
    plt.ylabel('Velocity, km/s')
    plt.xlabel('Distance, AU')

    plt.errorbar(array[:,0], array[:,1], fmt='x',
                markersize=6, alpha=0.5)
    
    plt.plot(fitted_data[:,0], fitted_data[:,1], label ='n = ' + 
            str('{:.3f}'.format(n)) +'$\pm$' +str('{:.2f}'.format(n_uncert))+ f'\ni = {i_deg}\n$\phi$ = {PA_deg}')
    
    plt.plot(keplarian[:,0], keplarian[:,1], label = 'Keplarian fit\nM = ' +
             str(mass) + ' $M_{\odot}$')
    
    plt.legend(loc = 'upper right',borderaxespad=0.5, frameon=True)
    plt.grid(True)
    plt.show()
    return

def plot_n_contours(results_array, incerz, posizers, name):
    # Extract n values from results_array
    n_values = results_array[:,:,1]  # Assuming n values are at index 1

    # Create a meshgrid for plotting
    i_mesh, pa_mesh = np.meshgrid(np.rad2deg(incerz), np.rad2deg(posizers), indexing='ij')

    # Create the contour plot
    plt.figure(figsize=(10, 6))
    cp = plt.contourf(i_mesh, pa_mesh, n_values, levels=20, cmap='viridis')
    plt.colorbar(cp, label='n value')
    plt.title('Contour Plot of n Values')
    plt.xlabel('Inclination Angle')
    plt.ylabel('Position Angle')

    plt.show()
    #plt.savefig('Contour_model.png', dpi=1000)

def n_PA(results_array, name):
    
    results_array = np.array([results_array])
    
    flattened_results = results_array.reshape(-1, results_array.shape[-1])
    #print(flattened_results)
    
    valid = (flattened_results[:,1] > -1 ) & (flattened_results[:,1] < 1)
    flattened_results = flattened_results[valid]
    
    n_values = flattened_results[:, 1]
    
    PA_values = np.rad2deg(flattened_results[:, 3])
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.scatter(PA_values, n_values, color='blue')
    #plt.plot(PA_values, n_values)
    plt.title(f'n vs. PA for {name}')
    plt.xlabel('Position Angle (PA)')
    plt.ylabel('Power Law Exponent (n)')
    plt.grid(True)
    plt.show()

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

def coordinate_transform(x,y,PA):
    
    x_trans = x*np.cos(PA)+y*np.sin(PA)
    y_trans = -x*np.sin(PA)+y*np.cos(PA)

    return x_trans,y_trans


def plotting(array, i, PA, fitted_data, mass, name, fitted_data_2, 
              mass_2):
    """
    r in column 0
    v in column 1
    r uncert in column 2
    v uncert in column 3
    
    """
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    i_deg = np.rad2deg(i)
    PA_deg = np.rad2deg(PA)
    
    ax.set_title('Velocity-Distance for '+ str(name), size=8)
    ax.set_ylabel('Velocity, km/s', size=8)
    ax.set_xlabel('Distance, AU', size=8)

    ax.errorbar(array[:,0], array[:,1], fmt='x',
                markersize=2, alpha=0.5)
    ax.plot(fitted_data[:,0], fitted_data[:,1], label ='M =' + 
            str('{:.3f}'.format(mass)) + ' M$_{\odot}$', color='purple')
    
    ax.plot(fitted_data_2[:,0], fitted_data_2[:,1], label ='M =' + 
            str('{:.3f}'.format(mass_2)) + ' M$_{\odot}$', color='orange')
    
    #ax.plot(fitted_data_3[:,0], fitted_data_3[:,1], label ='M =' + 
    #        str('{:.3f}'.format(mass_3)) + ' M$_{\odot}$', color='green')
    ax.plot([],[], label=f'i = {i_deg}$^o$\n$\phi$ = {PA_deg}$^o$',alpha=0)
    ax.legend(loc = 'upper right', fontsize=8, borderaxespad=0.5, frameon=False)
    
    
    return

def plotting_2(array, i, PA, fitted_data, mass, name, fitted_data_2, 
              mass_2):
    temp = []
    for line in array:
        if not line[0] > 3000:
            temp.append(line)
            
    array=np.vstack(temp)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    i_deg = np.rad2deg(i)
    PA_deg = np.rad2deg(PA)
    
    ax.set_title('Velocity-Distance for '+ str(name), size=8)
    ax.set_ylabel('Velocity, km/s', size=8)
    ax.set_xlabel('Distance, AU', size=8)

    ax.errorbar(array[:,0], array[:,1], fmt='x',
                markersize=2, alpha=0.5)
    ax.plot(fitted_data[:,0], fitted_data[:,1], label ='M =' + 
            str('{:.3f}'.format(mass)) + ' M$_{\odot}$', color='purple')
    
    ax.plot(fitted_data_2[:,0], fitted_data_2[:,1], label ='M =' + 
            str('{:.3f}'.format(mass_2)) + ' M$_{\odot}$', color='orange')
    
    #ax.plot(fitted_data_3[:,0], fitted_data_3[:,1], label ='M =' + 
    #        str('{:.3f}'.format(mass_3)) + ' M$_{\odot}$', color='green')
    ax.plot([],[], label=f'i = {i_deg}$^o$\n$\phi$ = {PA_deg}$^o$',alpha=0)
    ax.legend(loc = 'upper right', fontsize=8, borderaxespad=0.5, frameon=False,)
    
    return


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

def rotational_velocity(v_obs, x_trans, y_trans, i):

    x_flat = x_trans.flatten()
    y_flat = y_trans.flatten()
    
    x_mesh, y_mesh = np.meshgrid(x_flat, y_flat, indexing='ij')

    theta = np.arctan2(y_mesh, x_mesh)
    
    allowed_ranges = [(290, 360), (0, 70), (110, 250)]
    theta = filter_theta(theta, allowed_ranges)
    
    r = np.sqrt(x_mesh**2 + y_mesh**2)

    v_rot = v_obs / (np.sin(i) * np.cos(theta)) 
    v_rot = np.abs(v_rot)

    r_flat = r.flatten()
    v_rot_flat = v_rot.flatten()
    
    # Stack flattened arrays side by side to create the result
    result = np.vstack((r_flat, v_rot_flat)).T

    return result

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

def filter_data(data, val_1, val_2):

    first_column = data[:, 0]
    mask = (first_column >= val_1) & (first_column <= val_2)
    filtered_data = data[mask]
    
    return filtered_data

def filter_data_2(data, val_1, val_2):
    
    first_column = data[:, 1]
    mask = (first_column >= val_1) & (first_column <= val_2)
    filtered_data = data[mask]
    
    return filtered_data

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
        if lineo[1] < (mean + 30*std):
            temperz.append(lineo)
            
    return np.vstack((temperz))

def keplarian_fit(r,a):
    
    return a * r**-0.5

def fitting_data(data, minimum):
    
    x = data[:,0] - minimum
    #print(x,data[:,1])
    params, covariance = curve_fit(keplarian_fit, x, data[:,1],
                            absolute_sigma=True, maxfev=3000)
    return params, covariance

def find_max_value_coordinates(arr):
    max_value = np.nanmax(arr[:,1])
    max_y = np.where(arr == max_value)
    return max_y[0]

def find_mean_value_coordinates(arr):
    mean_value = np.nanmean(arr[:,1])
    mean_y = np.where(arr == mean_value)
    return mean_y[0]

def mass_(data):
    
    a = data[:,1] * data[:,0]**0.5
    a_scaled = a * 10**3 * au**0.5
    #print(a_scaled)
    midpoint = len(a_scaled) // 6
    first_half = a_scaled[:midpoint]
    a_epic = np.mean(first_half)
    masserz = a_epic**2 / (G*M_o)
    #print('masserz = ' + str(masserz))
    
    return masserz

def weighted_mean_calc(data):
    """
    hist, bins = np.histogram(data, bins=1000)
    # Calculate the weighted mean
    weighted_mean = np.average(bins[:-1], weights=hist)
    print("Weighted Mean:", weighted_mean)
    """
    weighted_mean = np.sum(data[:,0] * data[:,1]) / np.sum(data[:,1])
    #print("Weighted Mean:", weighted_mean)
    return weighted_mean

def decision(mean,maximum):
    
    if mean > maximum:
        value = maximum
    else:
        value = mean
    
    return value

def power_law(r, A, n):
    return A * r ** n

def trimming(array):
    temp = []
    for line in array:
        if not line[0] > 1750:
            temp.append(line)
    
    if not len(temp) ==0:
        temp=np.vstack(temp)
    return temp

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


get_data() #this gets data
