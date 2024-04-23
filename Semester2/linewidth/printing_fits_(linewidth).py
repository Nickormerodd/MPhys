# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 14:19:47 2024

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
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
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
RED = 'C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/S2/new_data/red.fits'
BLUE = 'C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/S2/new_data/blue.fits'
LINE = 'C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/S2/W4-6/XCLASS/CH3CN_10it/BestResult___parameter__V_width_Gauss___molecule__CH3CN_v=0____component__1___Iter__1.fits'

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
    
    hdul = fits.open(get_data())
    #hdul = fits.open(LINE)
    print('shape = ', hdul[0].shape)
    #hdul[0].data = match_shapes_with_nan(hdul[0].data)
    continuum = fits.open(CONTIN)
    red = fits.open(RED)
    blue = fits.open(BLUE)
    continuum_data_2d = continuum[0].data[0, 0, :, :]
    print(red[0].shape)

    
    plotting(hdul[0].data, continuum, blue, red, hdul)
    
    return

def plotting(line_map, continuum, blue, red, hdul):
    
    fig = plt.figure(figsize=(9,6))
    #print(hdul[0].header)
    
    continuum_2d = continuum[0].data[0, 0, :, :]
    red_2d = red[0].data[0,:,:]
    blue_2d = blue[0].data[0,:,:]
    x_m, y_m, x_gal, y_gal = galactic_coords(hdul[0], hdul[0].header)
    x_m_c, y_m_c, x_gal_c, y_gal_c = galactic_coords(continuum[0], continuum[0].header)
    x_m_r, y_m_r, x_gal_r, y_gal_r = galactic_coords(red[0], red[0].header)
    x_m_b, y_m_b, x_gal_b, y_gal_b = galactic_coords(blue[0], blue[0].header)
    xo_pix, yo_pix = 23, 24

    xx, yy = np.meshgrid(x_gal, y_gal)
    xx_c, yy_c = np.meshgrid(x_gal_c, y_gal_c)
    xx_r, yy_r = np.meshgrid(x_gal_r, y_gal_r)
    xx_b, yy_b = np.meshgrid(x_gal_b, y_gal_b)
    
    
    mask = (0.028 > red_2d) | (red_2d > 0.0448)
    red_2d[mask] = np.nan
    
    
    pig = plt.contourf(xx,yy,line_map, cmap='jet', levels=100)
    plt.contour(xx_c,yy_c,continuum_2d, levels = 25, alpha = 0.5, cmap='Greys',
                linewidths = 1, label = 'Continuum contours')
    #red_cont = plt.contour(xx_r, yy_r, red_2d, levels=8, alpha=0.5, cmap='Reds',
    #                   linewidths=2, label='Red-shifted $^{13}$CO outflow')
    #blue_cont = plt.contour(xx_b,yy_b,blue_2d, levels = 7, alpha = 0.6, cmap='Blues',
     #           linewidths = 2)
    #plt.clabel(red_cont, fontsize=8)
    plt.colorbar(pig, label='Ratio')
    plt.xlabel("GLON", size=12)
    plt.ylabel("GLAT", size=12)

    
    plt.title('CH$_{3}$CN - CH$_{3}^{\: 13}$CN Ratio', size=13)
    name = OUT + 'CH3CN_CH3-13-CN_ratio.png'
    
    # Get beam parameters from header
    bmaj = hdul[0].header['BMAJ'] * 3600  # Convert from degrees to arcseconds
    bmin = hdul[0].header['BMIN'] * 3600
    bpa = hdul[0].header['BPA']
    
    # Add beam ellipse
    beam = Ellipse(xy=(0.1, 0.1), width=2*bmin, height=2*bmaj, angle=bpa, transform=plt.gca().transAxes,
                   fill=True, color='black', alpha=0.5)
    plt.gca().add_patch(beam)
    plt.gca().text(0.2, 0.2,'', color='white', ha='center', va='center', transform=plt.gca().transAxes)
    
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
    
    #plt.legend(handles=[continuum_legend, red_outflow_legend], fontsize=11, loc='lower right')   
    plt.legend(handles=[continuum_legend], fontsize=11, loc='lower right')
    #plt.savefig(name, bbox_inches='tight', dpi =300)
    plt.show()
    
    return


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
        #x_gal = x_gal[::-1]
        #y_gal = y_gal[::-1]
        
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
    
    