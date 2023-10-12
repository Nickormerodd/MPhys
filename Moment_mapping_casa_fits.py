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
from scipy.constants import c
import os

def get_data():

    Tk().withdraw()  
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        filename = (os.path.basename(file_path).replace(".fits", ""))#.replace("_"," ").replace("."," ")
        main(data, filename, file_path)
        data.close()    
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
    
    name = name.replace("_image","") + ".png"
    print(name)
    
    calibration_freq = 220.70908 #GHz
    
    if '914' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{2}\\rightarrow11_{2}$'
        conversion = (c/1000) * (220.73027 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion
    elif '958' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{3}\\rightarrow11_{3}$'
    elif '1020' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{4}\\rightarrow11_{4}$'
        conversion = (c/1000) * (220.679114 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion
    elif '1096' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{5}\\rightarrow11_{5}$'
        conversion = (c/1000) * (220.64112 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion
    else:
        title = 'unknown'
    
    contourf = ax.contourf(X, Y, data[0].data[0] , cmap='bwr',
                           levels = 50) # 'viridis')
    plt.colorbar(contourf, label='km/s')
    ax.set_title(title)
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    plt.savefig(name, dpi=1000)
    
    return

get_data()
