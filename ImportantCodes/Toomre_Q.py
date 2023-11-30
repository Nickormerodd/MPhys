"""
Created on Thu Nov 30 14:08:05 2023

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
#from scipy.constants import c
from astropy import units as u
from radio_beam import Beam
#import uncertainties.unumpy as unp
#from uncertainties import ufloat
#from astropy.convolution import equivalencies

G = 6.6743e-11
au = 1.495979e+11 #AU to m
M_o = 1.98847e+30
pc = 3.086e+16
#k = 1.3806488e-23
R = 8.3144621
A_kep = 103.70690152932114 # km/s AU**0.5
A_kep_SI = A_kep * 10**3 * au**0.5

#bbarolo_path = 'C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/week_7-9/velocity_distance/bbarolo_data/densprof.txt'

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
    
    print(data[0].shape)
    #print(data[0].data[0])
    
    m = data[0].shape[2] # y shape
    p = data[0].shape[3] # x shape
    
    header = data[0].header
    ref_freq = header['CRVAL3'] #reference frequency
    incr_freq = header['CDELT3'] #frequency increment
    pix_1 = np.deg2rad(header['CDELT1'])
    pix_2 = np.deg2rad(header['CDELT2'] )
    pixel_size = pix_1*pix_2
    pixel_area = 0.5 * pixel_size**2 * (8.1 * 10**3 * pc)**2
    beam_major = np.deg2rad(header['BMAJ'])  # Major axis of the beam 
    beam_minor = np.deg2rad(header['BMIN'])  # Minor axis of the beam
    beam_sizer = (np.pi/(4 *np.log(2))) * beam_major * beam_minor
    print(beam_minor, beam_major, beam_sizer)
    beam_size = beam_sizer * u.sr
    data[0].data[0] = (data[0].data[0]*u.Jy/u.beam).to(u.Jy/u.sr, equivalencies=u.beam_angular_area(beam_size))
    F_v = data[0].data[0] * beam_sizer
    print(str(np.nanmax(F_v)) + ' Jansky')
    temps = temp_calc()
    wave = 1.36
    
    xo_pix,yo_pix = 15,23
    x_m, y_m, x_gal, y_gal = galactic_coords(data[0], header)
    xo_m = x_m[xo_pix] #central xo point
    yo_m = y_m[yo_pix] #central yo point
    x_norm = x_m - xo_m #normalising around the xo point
    y_norm = y_m - yo_m #normalising around the yo point
    
    ang_mom = np.zeros((m, p))
    M = np.zeros((m, p))
    for i in range(m): # i is y
        for j in range(p): #  j is x
            #print(F_v[0,i,j])
            ang_mom[i,j] = ang_mom_calc(x_norm[j], y_norm[i])
            M[i,j] = mass_calc(wave, temps, F_v[0,i,j])
    
    print(str(np.nanmax(M)) + ' max sol mass')
    print(str(np.nansum(M)) + ' sum sol mass')
    
    surf_density = (M * M_o)/pixel_area 
    print(pixel_area)
    #print(np.nanmax(surf_density))
    #print(np.nanmin(surf_density))
    
    print(surf_density.shape)
    print(ang_mom.shape)
    cs = speed_sound(temps)
    Q = Toomre_Q(cs, ang_mom, surf_density)
    print('Q max = ' +str(np.nanmax(Q)))
    plotting(Q,x_gal,y_gal)
    return

def plotting(Q, x_gal, y_gal):
    # Create a contour plot with Q as the contour levels
    plt.contour(x_gal, y_gal, Q, levels = 50)
    
    # Add labels, titles, colorbar, etc. as needed
    plt.xlabel("GLON deg")
    plt.ylabel("GLAT deg")
    #plt.invert_xaxis() 
    plt.title('Contour Map of Q')
    plt.colorbar(label='Q Contour Levels')
    
    # Show the plot
    plt.show()

def mass_calc(wave,T,F_v):
    """
    wavelength : in mm
    T : in K
    k_v : 1.99 cm^2 g^-1
    F_v : From contintuum in Jy
    d : in pc, (8.1 kpc)
    M : solar masses
    """
    d = 8.1*10**3
    k_v = 1.99
    M = 0.12*(np.exp(1.439*(wave*(T/10))**(-1) -1))*((k_v/0.01)**(-1))*(F_v)*((d/100)**2)*((wave)**3)
    
    return M

def temp_calc():
    
    temp = 200
    
    return temp

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
    x_gal = (crval1 + (np.arange(data.shape[3]) - crpix1) * cdelt1)
    y_gal = (crval2 + (np.arange(data.shape[2]) - crpix2) * cdelt2)
    
    x_m = x_gal * 8.1 * 10**3 * 206264.81 * np.pi/180 * au
    y_m = y_gal * 8.1 * 10**3 * 206264.81 * np.pi/180 * au
    
    return x_m,y_m, x_gal, y_gal

def Toomre_Q(cs,ang_mom,E):
    
    Q = cs*ang_mom / (G*E*np.pi)
    
    return Q

def ang_mom_calc(x,y):
    
    r = np.sqrt(x**2 + y**2)
    
    if r == 0:
        ang_mom = np.nan
    else:
        try:
            v_rot = A_kep_SI * r **-0.5
            ang_mom = v_rot / r
    
        except RuntimeWarning:
            ang_mom = np.nan
    
    return ang_mom

def speed_sound(T):
    
    Molar = 2.8 * 10e-3
    gam = 1.4
    cs = np.sqrt(gam*R*T/Molar)
    
    return cs

get_data()
