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
    
    #four_plot(file_paths)
    
    return
                                     

def main(data, name, path):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    m = data[0].shape[1]
    p = data[0].shape[2]

    #print(data[0].shape)
    header = data[0].header
    #print(header)
    #    Extract relevant header values
    ctype1 = header['CTYPE1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    crpix1 = header['CRPIX1']
    
    ctype2 = header['CTYPE2']
    crval2 = header['CRVAL2']
    cdelt2 = header['CDELT2']
    crpix2 = header['CRPIX2']
    
    # Create the coordinate arrays
    x = (crval1 + (np.arange(data[0].shape[2]) - crpix1) * cdelt1)
    y = (crval2 + (np.arange(data[0].shape[1]) - crpix2) * cdelt2)

    median = np.nanmedian(data[0].data[0])
    print(median)

    name = name.replace("_image","") + "_gal" +".png"
    print(name)

    calibration_freq = 220.70908 #GHz

    if '875' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{0}\\rightarrow11_{0}$'
        conversion = (c/1000) * (220.747369 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion
    elif '890' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{1}\\rightarrow11_{1}$'
        conversion = (c/1000) * (220.743012 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion
    elif '914' in name.lower():
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
    elif '1190' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{6}\\rightarrow11_{6}$'
        conversion = (c/1000) * (220.593987 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion
    elif '1300' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{7}\\rightarrow11_{7}$'
        conversion = (c/1000) * (220.538635 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion
    elif '1430' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{8}\\rightarrow11_{8}$'
        conversion = (c/1000) * (220.474767 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion
    else:
        title = 'unknown'

    contourf = ax.contourf(x, y, data[0].data[0] , cmap='bwr',
                           levels = 50) # 'viridis')
    plt.colorbar(contourf, label='km/s')
    ax.set_title(title)

    #ax.set_xticks([])
    #ax.set_yticks([])
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    
    ax.set_xlabel(ctype1)
    ax.set_ylabel(ctype2)
    #plt.title('Contour Plot')
    ax.invert_xaxis()  # Invert the x-axis to match the standard Galactic coordinate system
    #ax.invert_yaxis()  # Invert the y-axis
    ax.grid(alpha=0.2)
    
    plt.savefig(name, dpi=1000, bbox_inches = 'tight')

    return
