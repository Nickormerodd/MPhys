"""
Created on Tue Oct 10 15:57:51 2023

@author: Christopher
"""

import numpy as np
from astropy.io import fits
from spectral_cube import SpectralCube as sc
from tkinter import Tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
import radio_beam as beam
from scipy.constants import c
from astropy import units as u
from scipy.optimize import curve_fit
from matplotlib.ticker import FuncFormatter
import uncertainties.unumpy as unp
from uncertainties import ufloat


def get_data():

    Tk().withdraw()  
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        filename = (os.path.basename(file_path).replace(".fits", "")).replace("_"," ").replace("."," ")
        main(data, filename, file_path)
        data.close()    
    return

def main2(data, name, path):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    m = data[0].shape[1]
    p = data[0].shape[2]
    
    print(data[0].shape)

    x = np.linspace(0, p - 1, p)
    y = np.linspace(0, m - 1, m)
    X, Y = np.meshgrid(x, y)
    
    median = np.nanmedian(data[0].data[0])
    print(median)
    
    ref_velocity = 0  # You can adjust this value as needed
    data[0].data[0] = data[0].data[0] + ref_velocity
    
    # Create a colormap with blue for values less than median and red for values greater than median
    cmap = plt.cm.get_cmap('bwr')  # You can choose any colormap you prefer
    
    contourf = ax.contourf(X, Y, data[0].data[0], cmap=cmap, levels=30, 
                           vmin=data[0].data[0].min(), 
                           vmax=data[0].data[0].max())
    plt.colorbar(contourf, label='km/s')
    
    #contour = ax.contour(X, Y, data[0].data[0], cmap='magma', levels=10, alpha=0.01)
    ax.set_title('')
    
    return




def main(data, name, path):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    m = data[0].shape[1]
    p = data[0].shape[2]
    
    print(data[0].shape)

    x = np.linspace(0, p - 1, p)
    y = np.linspace(0, m - 1, m)
    X, Y = np.meshgrid(x, y)
    
    median = np.nanmedian(data[0].data[0])
    print(median)
    
    ref_velocity =  30.02346263833746
    data[0].data[0] = data[0].data[0] + ref_velocity #- np.nanmean(data[0].data[0])
    
    contourf = ax.contourf(X, Y, data[0].data[0] , cmap='bwr',
                           levels = 30) # 'viridis')
    plt.colorbar(contourf, label='km/s')
    #contour = ax.contour(X,Y,data[0].data[0], cmap = 'black', levels = 10,
     #                    alpha = 0.1)
    ax.set_title('$CH_{3}CN-J_{K}=12_{2}\\rightarrow11_{2}$')
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    plt.savefig('SPW25_ch914.png', dpi=1000)
    
    
    return

get_data()
