# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 13:59:24 2024

@author: Christopher
"""

import numpy
import matplotlib.pyplot as plt
from scipy.stats import iqr
from astropy.io import fits
#from matplotlib import gridspec
import numpy as np
#import networkx as nx

from tkinter import Tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.patches import Ellipse
from matplotlib.ticker import FormatStrFormatter
from astropy.wcs import WCS
from matplotlib.lines import Line2D
#from reproject import reproject_interp
#from mpl_toolkits.mplot3d import Axes3D as mpl
#from matplotlib.lines import Line2D
#from scipy.constants import c
#import uncertainties.unumpy as unp
#from uncertainties import ufloat

G = 6.6743e-11
au = 1.495979e+11
M_o = 1.98847e+30
OUT = 'C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/S2/W7-9/linewidth/'
CONTIN = 'C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/S2/new_data/continuum_Xreg.fits'
TEMP = 'C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/S2/W4-6/XCLASS/CH3-13CN_6it/T_rot_CH3-13CN_test_8.fits'
RED = 'C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/S2/new_data/red.fits'

def get_data():

    Tk().withdraw()
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        #data = fits.open(file_path)
        #filename = (os.path.basename(file_path)).replace(".fits", "")#.replace(".image.pbcor_line.galactic","")
        #main(data)
        #data.close()
        print(file_path)
        
        return file_path
    
    
def main():
    
    #temp_hdul = fits.open(get_data())
    temp_hdul = fits.open(TEMP)
    continuum = fits.open(CONTIN)
    red = fits.open(RED)
    continuum_data_2d = continuum[0].data[0, 0, :, :]


    #print('shape = ', temp_hdul[0].shape)
    print('contin shape = ', continuum_data_2d.shape)

    #temp_hdul[0].data = match_shapes_with_nan(temp_hdul[0].data)
    
    line_map, av_line_map = line_eq(temp_hdul[0].data)
    
    plotting(line_map, temp_hdul, red, continuum, 0)
    plotting(av_line_map, temp_hdul, red, continuum, 1)
    
    return


def plotting(line_map, temp_hdul, red, continuum, n):
    fig = plt.figure(figsize=(9,6))
    
    #continuum[0].data[0][0] = continuum[0].data[0, 0, :, :]
    red_2d = red[0].data[0,:,:]
    continuum_2d = continuum[0].data[0, 0, :, :]
    x_m, y_m, x_gal, y_gal = galactic_coords(temp_hdul[0], temp_hdul[0].header)
    x_m_c, y_m_c, x_gal_c, y_gal_c = galactic_coords(continuum[0], continuum[0].header)
    x_m_r, y_m_r, x_gal_r, y_gal_r = galactic_coords(red[0], red[0].header)

    xx, yy = np.meshgrid(x_gal, y_gal)
    xx_r, yy_r = np.meshgrid(x_gal_r, y_gal_r)
    xx_c, yy_c = np.meshgrid(x_gal_c, y_gal_c)
    
    mask = (0.028 > red_2d) | (red_2d > 0.0448)
    red_2d[mask] = np.nan
    
    pig = plt.contourf(xx, yy, line_map, cmap='bwr', levels=100)
    plt.contour(xx_c,yy_c,continuum_2d, levels = 25, alpha = 0.3, cmap='copper',
                linewidths = 1)
    red_cont = plt.contour(xx_r, yy_r, red_2d, levels=8, alpha=0.5, cmap='Reds',
                       linewidths=2, label='Red-shifted $^{13}$CO outflow')
    cbar = plt.colorbar(pig, label='km/s')
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.xlabel("GLON", size=12)
    plt.ylabel("GLAT", size=12)
    
    if n == 0:
        plt.title('Thermal Linewidth ($\sigma_{th,m}$) for CH$_{3}^{\: 13}$CN', size=14)
        name = OUT + 'Thermal_linewidth_contour_CH3-13-CN.pdf'
        name2 = OUT + 'Thermal_linewidth_CH3-13-CN.fits'
    elif n == 1:
        plt.title('Thermal Linewidth, $\sigma_{th,m_{av}}}}$', size=14)
        name = OUT + 'Average_Thermal_contour_linewidth.pdf'
        name2 = OUT + 'Average_Thermal_linewidth.fits'
    
    # Get beam parameters from header
    bmaj = temp_hdul[0].header['BMAJ'] * 3600  # Convert from degrees to arcseconds
    bmin = temp_hdul[0].header['BMIN'] * 3600
    bpa = temp_hdul[0].header['BPA']
    
    # Add beam ellipse
    beam = Ellipse(xy=(0.1, 0.1), width=2*bmin, height=2*bmaj, angle=bpa, transform=plt.gca().transAxes,
                   fill=True, color='black', alpha=0.5)
    plt.gca().add_patch(beam)
    plt.gca().text(0.1, 0.1,'', color='white', ha='center', va='center', transform=plt.gca().transAxes)
    
    num_ticks = 6

    # Set the tick labels on both axes
    x_ticks = np.linspace(x_gal.min(), x_gal.max(), num_ticks)
    y_ticks = np.linspace(y_gal.min(), y_gal.max(), num_ticks)
    
    # Format tick labels nicely
    x_tick_labels = [f"{tick:.5f}" for tick in x_ticks]
    y_tick_labels = [f"{tick:.5f}" for tick in y_ticks]

    # Set the ticks and labels
    plt.xticks(ticks=x_ticks, labels=x_tick_labels)
    plt.yticks(ticks=y_ticks, labels=y_tick_labels)
    plt.gca().invert_xaxis()
    
    continuum_legend = Line2D([0], [0], color='gray', lw=1, label='Continuum contours')
    red_outflow_legend = Line2D([0], [0], color='red', lw=2, label='Red-shifted $^{13}$CO outflow')
    
    plt.legend(handles=[continuum_legend, red_outflow_legend], fontsize=11, loc='lower right')   

    
    plt.savefig(name, bbox_inches='tight')
    plt.show()
    
    print(name2)
    
    header = temp_hdul[0].header
    
    hdu = fits.PrimaryHDU(line_map, header=header)

    # Create an HDU list and append the Primary HDU
    hdul = fits.HDUList([hdu])
    
    

    #hdul.writeto(name2, overwrite=True)
    return

def line_eq(temp_map):
    
    m = 42
    m_av = 2.33
    
    therm_linewidth = 288 * ( (42)**(-1/2) ) * ( (temp_map/10)**(1/2) ) * 10**(-3) #km/s
    av_linewidth = 288 * ( (2.33)**(-1/2) ) * ( (temp_map/10)**(1/2) ) * 10**(-3) #km/s
    
    return therm_linewidth, av_linewidth

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
    if data.data.ndim == 4:
        x_gal = (crval1 + (np.arange(data.shape[3]) - crpix1) * cdelt1)
        y_gal = (crval2 + (np.arange(data.shape[2]) - crpix2) * cdelt2)        
        
    elif data.data.ndim ==2:
        x_gal = (crval1 + (np.arange(data.shape[1]) - crpix1) * cdelt1)
        y_gal = (crval2 + (np.arange(data.shape[0]) - crpix2) * cdelt2)
        
    elif data.data.ndim ==3:
        x_gal = (crval1 + (np.arange(data.shape[2]) - crpix1) * cdelt1)
        y_gal = (crval2 + (np.arange(data.shape[1]) - crpix2) * cdelt2)

    x_m = x_gal * 8.1 * 10**3 * 206264.81 * np.pi/180 * au
    y_m = y_gal * 8.1 * 10**3 * 206264.81 * np.pi/180 * au

    return x_m,y_m, x_gal, y_gal

def match_shapes_with_nan(data):
    # Extract the shape of the 2D data array (ignoring the first singleton dimension)
    if data.ndim == 4:
        n, m = data.shape[2], data.shape[3]
        max_dim = max(n, m)
        # Create a new data array filled with NaNs of the maximum dimension
        new_data = np.full((1, 1, max_dim, max_dim), np.nan)
        
        # Determine the padding amounts
        pad_n = (max_dim - n) // 2
        pad_m = (max_dim - m) // 2

        end_n = pad_n + n
        end_m = pad_m + m

        # Ensure that the original data is correctly copied into the padded array
        new_data[0, 0, pad_n:end_n, pad_m:end_m] = data[0, 0, :, :]
        
    elif data.ndim ==2:
        n, m = data.shape[0], data.shape[1]
                
        max_dim = max(n, m)
        
        # Create a new data array filled with NaNs of the maximum dimension
        new_data = np.full((max_dim, max_dim), np.nan)
        
        # Determine the padding amounts
        pad_n = (max_dim - n) // 2
        pad_m = (max_dim - m) // 2
        
        end_n = pad_n + n
        end_m = pad_m + m
        
        # Ensure that the original data is correctly copied into the padded array
        new_data[pad_n:end_n, pad_m:end_m] = data[:, :]
        

    return new_data

def coord_trans(x,y):
    phi = np.deg2rad(170)
    inc = np.deg2rad(38)
    #z=
    alpha = phi + np.pi
    
    x_prime = x*np.cos(alpha) + y*np.sin(alpha)
    y_prime = x*np.sin(alpha)/np.cos(inc) + y*np.cos(alpha)/np.cos(inc)

    return x_prime,y_prime

main()
    
    
