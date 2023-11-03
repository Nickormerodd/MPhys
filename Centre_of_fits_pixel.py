"""
Created on Fri Nov  3 14:40:36 2023

@author: Christopher
"""

import numpy as np
from astropy.io import fits
from tkinter import Tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt


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
    
    temp = temp / num_channels
    max_x,max_y = find_max_value_coordinates(temp)
    print('epic')
    print('centre (x,y) = ' +str(max_x) + ',' + str(max_y) + ' pixels')
    
    m, p = temp.shape  
    X, Y = np.meshgrid(range(p), range(m))

    surf = ax.plot_surface(X, Y, temp, cmap='magma', linewidth=0, antialiased=False)

    # Customize the plot (labels, title, etc.) as needed
    #ax.set_xlabel('X')
    
    #ax.set_ylabel('Y')
    #ax.set_zlabel('Z')
    ax.set_title('Location of Protostar')
    
    # Add a colorbar
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='K')
    
    return 


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

get_data()
