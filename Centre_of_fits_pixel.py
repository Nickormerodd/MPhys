"""
Created on Sat Nov  4 22:37:22 2023

@author: Christopher
"""


import numpy as np
from astropy.io import fits
from tkinter import Tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def get_data():

    Tk().withdraw()
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        filename = (os.path.basename(file_path).replace(".fits", ""))
        main(data, filename, file_path)
        #main_2(data)
        data.close()
        
    return

def main(data, filename, path):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    num_channels = data[0].shape[0]
    shape = (data[0].shape[1], data[0].shape[2])
    temp = np.zeros(shape)

    for i in range(num_channels):
        temp = temp + data[0].data[i]

    # Clean the data by replacing NaNs with zeros
    temp = temp/num_channels
    temp[np.isnan(temp)] = 0

    m, p = temp.shape
    X, Y = np.meshgrid(range(p), range(m))
    xdata = np.column_stack((X.ravel(), Y.ravel()))

    def func(xdata, A, mean_x, mean_y, sigma):
        x, y = xdata[:, 0], xdata[:, 1]
        return A * np.exp(-((x - mean_x) ** 2 + (y - mean_y) ** 2) / (2 * sigma ** 2))

    # Initial guess for the parameters
    initial_guess = (temp.max(), np.argmax(temp) % temp.shape[1], np.argmax(temp) // temp.shape[1], 1.0)

    popt, _ = curve_fit(func, xdata, temp.ravel(), p0=initial_guess)

    # Extract the (x, y) coordinates of the Gaussian fit's maximum
    new_max_x, new_max_y = int(popt[1]), int(popt[2])

    max_x,max_y = find_max_value_coordinates(temp)

    print('Position of maximum value (x, y) = ' + str(max_x) + ', ' + str(max_y) + ' pixels')
    print('Position of maximum value from Gaussian fit (x, y) = ' + str(new_max_x) + ', ' + str(new_max_y) + ' pixels')

    # Generate the fitted surface using the parameters obtained
    fitted_surface = func(xdata, *popt).reshape(m, p)

    surf = ax.plot_surface(X, Y, fitted_surface, cmap='viridis', linewidth=0, antialiased=False)

    # Customize the plot (labels, title, etc.) as needed
    ax.set_title('Location of Protostar - Gaussian Fit')

    # Add a colorbar
    #fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='K')

    plt.show()




def find_max_value_coordinates(arr):
    max_value = np.nanmax(arr)
    max_x, max_y = np.where(arr == max_value)
    return max_x[0], max_y[0]

def main_2(data):
    
    num_channels = data[0].shape[0]
    shape = (data[0].shape[1], data[0].shape[2])
    temp = np.zeros(shape)
    
    for i in range(num_channels):
        temp = temp + data[0].data[i]
        #temp.append([max_x, max_y])
    temp = temp/num_channels
    max_x,max_y = find_max_value_coordinates(temp)
    
    print('epic')
    print(max_x,max_y)
    return

def gaussian(x, A, mean, sigma):
    return A * np.exp(-(x - mean)**2 / (2 * sigma**2))


get_data()
