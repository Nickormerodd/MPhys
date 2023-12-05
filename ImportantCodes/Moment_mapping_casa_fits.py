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
#from matplotlib.patches import FancyArrrowPatch as fap

def get_data():

    Tk().withdraw()
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        filename = (os.path.basename(file_path).replace(".fits", ""))#.replace("_"," ").replace("."," ")
        main(data, filename, file_path)
        data.close()
        print(file_path)

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

    calibration_freq = 220.70908 #GHz

    if '875' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{0}\\rightarrow11_{0}$'
        k_number = 0
        conversion = (c/1000) * (220.747369 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '890' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{1}\\rightarrow11_{1}$'
        k_number = 1
        conversion = (c/1000) * (220.743012 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '914' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{2}\\rightarrow11_{2}$'
        k_number = 2
        conversion = (c/1000) * (220.73027 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '958' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{3}\\rightarrow11_{3}$'
        k_number = 3

    elif '1020' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{4}\\rightarrow11_{4}$'
        k_number = 4
        conversion = (c/1000) * (220.679114 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '1096' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{5}\\rightarrow11_{5}$'
        k_number = 5
        conversion = (c/1000) * (220.64112 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '1190' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{6}\\rightarrow11_{6}$'
        k_number = 6
        conversion = (c/1000) * (220.593987 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '1300' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{7}\\rightarrow11_{7}$'
        k_number = 7
        conversion = (c/1000) * (220.538635 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '1430' in name.lower():
        title = '$CH_{3}CN-J_{K}=12_{8}\\rightarrow11_{8}$'
        k_number = 8
        conversion = (c/1000) * (220.474767 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '_CH3_13CN_K=1' in name.upper():
        title = '$CH_{3}(13)CN-J_{K}=12_{1}\\rightarrow11_{1}$'
        k_number = 1
        conversion = (c/1000) * (220.63384 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '_CH3_13CN_K=2' in name.upper():
        title = '$CH_{3}(13)CN-J_{K}=12_{2}\\rightarrow11_{2}$'
        k_number = 2
        conversion = (c/1000) * (220.62113 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '_CH3_13CN_K=3' in name.upper():
        title = '$CH_{3}(13)CN-J_{K}=12_{3}\\rightarrow11_{3}$'
        k_number = 3
        conversion = (c/1000) * (220.6000 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '_CH3_13CN_K=4' in name.upper():
        title = '$CH_{3}(13)CN-J_{K}=12_{4}\\rightarrow11_{4}$'
        k_number = 4
        conversion = (c/1000) * (220.57032 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    elif '_CH3_13CN_K=5' in name.upper():
        title = '$CH_{3}(13)CN-J_{K}=12_{5}\\rightarrow11_{5}$'
        k_number = 5
        conversion = (c/1000) * (220.53224 - calibration_freq)/calibration_freq
        data[0].data[0] = data[0].data[0] + conversion

    else:
        title = 'unknown'
        k_number = 'unknown'

    if 'CH3_13CN' in name.upper():
        name = "CH3(13)CN_k=" + str(k_number) + "_gal" +".png"
    else:
        name = "CH3CN_k=" + str(k_number) + "_gal" +".png"
    print(name)

    x_au = x * 8.1 * 10**3 * 206264.81 * np.pi/180
    y_au = y * 8.1 * 10**3 * 206264.81 * np.pi/180

    x_norm = np.abs(1000 / (8.1 * 10**3 * 206264.81 * np.pi/180))
    y_norm = np.abs(1000 / (8.1 * 10**3 * 206264.81 * np.pi/180))

    #print(len(y))
    #print(np.mean(y))

    x_point = np.min(x) + +(np.mean(x)/(len(x)*160))
    y_point = np.min(y)+(np.mean(y)/(len(y)*150))
    x_values = np.linspace(x_point, x_point+x_norm,2)
    y_values = np.linspace(y_point,y_point,2)
    ax.plot(x_values,y_values,color='black')#, label='1000 AU')
    ax.plot([],[],label='1000 AU', color='black', alpha=0)

    #y_au = 1000
    y_norm = 1000 / (8.1 * 10**3 * 206264.81 * np.pi/180)
    x_center = x_point + x_norm / 2  # Center of the x-direction line
    y_point = np.min(y) - y_norm / 2  # Center of the y-direction line
    x_values_y = np.linspace(x_center, x_center, 2)
    y_values_y = np.linspace(y_point, y_point + y_norm, 2)
    ax.plot(x_values_y, y_values_y, color='black')

    arrow_props = dict(arrowstyle='->', linewidth=1, color='black')
    contourf = ax.contourf(x, y, data[0].data[0] , cmap='bwr',
                           levels = 50) # 'viridis')
    plt.colorbar(contourf, label='km/s')
    ax.set_title(title)
    #centre = 48,46
    print(data[0].shape)
    #ax.scatter(x[48],y[46],label='centre')


    ax.set_xlabel("GLON deg")
    ax.set_ylabel("GLAT deg")

    ax.invert_xaxis()  # Invert the x-axis to match the standard Galactic coordinate system
    #ax.invert_yaxis()  # Invert the y-axis
    ax.grid(alpha=0.2)
    ax.legend(loc='lower right', fontsize=7.5, borderaxespad=0.5, frameon=False, edgecolor='black')

    plt.savefig(name, dpi=1000, bbox_inches='tight')

    return


get_data()
