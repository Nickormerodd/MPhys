# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 12:12:06 2023

@author: nickl
"""

import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
from tkinter import Tk, filedialog
from astropy.io import fits
#from astropy.modeling import models, fitting
#from mpl_toolkits.mplot3d import Axes3D
#from astropy.wcs import WCS
from astropy import units as u
#from radio_beam import Beam

############################## global params ##################################
pc = 3.086e+16
G = 6.6743e-11
SunMass = 1.98847e30
au = 1.496e+11
distance = 8.1*10**3  # Example distance in parsecs
kappa_nu = 1.99  # Dust opacity in cm^2/g
lambda_mm = 1.3  # Wavelength in millimeters
temperature_K = 50  # Temperature in Kelvin
doi = 1 #times the FWHM beams by this much
SaveFig = True
SaveFigName = 'C:/Users/nickl/linux/week10/massbefore.png'
SaveFigName2 = 'C:/Users/nickl/linux/week10/MassSubtraction.png'
initial_dir = 'C:/Users/nickl/linux/data/Fits_Jansky' #where is your rings_final2.txt file
#############################################################################

def get_data(initial_dir):
    # file picking with a temp window
    root = Tk()
    root.attributes('-topmost', True)  # Set the window to be always on top

    # Prompt the user to select files, starting in the specified directory
    file_path = filedialog.askopenfilename(initialdir=initial_dir, title="Select Folder, then select rings_final2.txt")

    root.destroy() # removes temp window

    hdu_list = fits.open(file_path)
    data = hdu_list[0].data
    header = hdu_list[0].header
    hdu_list.close()
    return data, header

def flux_calc(data, header):

    print(data[0].shape) # flux (1), then y, then x

    m = data[0].shape[1] # y shape
    p = data[0].shape[2] # x shape
   # ref_freq = header['CRVAL3'] #reference frequency
    #incr_freq = header['CDELT3'] #frequency increment
    beam_major = np.deg2rad(header['BMAJ'])  # Major axis of the beam
    beam_minor = np.deg2rad(header['BMIN'])  # Minor axis of the beam
    beam_sizer = (np.pi/(4 *np.log(2))) * beam_major * beam_minor
    #print(f'{beam_minor:.3e} {beam_major:.3e} {beam_sizer:.3e} minor, major and size in rad, rad, sr')
    data_quantity = data[0] * u.Jy / u.beam
    data_converted = data_quantity.to(u.Jy / u.sr, equivalencies=u.beam_angular_area(beam_sizer * u.sr))

    F_v = data_converted.value * beam_sizer # in Jy
    print(f'{np.nanmax(F_v):.3f} max flux in Jansky in one pixel')
    """
    xo_pix,yo_pix = 15,23
    x_m, y_m, x_gal, y_gal = galactic_coords(data[0], header)
    xo_m = x_m[xo_pix] #central xo point
    yo_m = y_m[yo_pix] #central yo point
    x_norm = x_m - xo_m #normalising around the xo point
    y_norm = y_m - yo_m #normalising around the yo point
    """
    M = np.zeros((m, p))
    z = np.zeros((m, p))
    for i in range(m): # i is y
        for j in range(p): #  j is x
            M[i, j] = mass_calc(F_v[0, i, j]*u.Jy) #takes 'mass' at each point
            z[i,j] = F_v[0,i,j] #takes Jy at each point

    print(f'{np.nansum(z):.3f} sum of Flux in Jy')
    print(f'{np.nanmax(M):.3f} max sol mass in one pixel')
    print(f'{np.nansum(M):.3f} sum sol mass')
    #print(np.nanmax(surf_density))
    #print(np.nanmin(surf_density))

    return z

def Jy_image_and_plot(data, header, z):

    bmaj = np.deg2rad(header['BMAJ'])  # Beam major axis in radians
    bmin = np.deg2rad(header['BMIN'])  # Beam minor axis in radians
    beam_area_sr = np.pi * (bmaj * bmin) / (4 * np.log(2))  # Beam area in steradians

    pixel_scale_lon = abs(header['CDELT1'])
    pixel_scale_lat = abs(header['CDELT2'])
    beam_major_pixels = header['BMAJ'] / pixel_scale_lon
    beam_minor_pixels = header['BMIN'] / pixel_scale_lat
    fwhm_stddev_maj_y = doi*beam_major_pixels / (2 * np.sqrt(2 * np.log(2)))
    fwhm_stddev_min_x = doi*beam_minor_pixels / (2 * np.sqrt(2 * np.log(2)))
    print(f'BMAJ = {fwhm_stddev_maj_y:.3f}, BMIN = {fwhm_stddev_min_x:.3f}')

    flux_density = data
    flux_density_Jy = (data * u.Jy / u.beam).to(u.Jy / u.sr, equivalencies=u.beam_angular_area(beam_area_sr * u.sr)) # Convert flux from Jy/beam to Jy/sr
    Jy_image = flux_density_Jy * beam_area_sr *u.sr #changinf Jy/sr to Jy

    x_stddev =  fwhm_stddev_min_x # Convert FWHM to standard deviation
    y_stddev = fwhm_stddev_maj_y

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y = np.meshgrid(np.arange(data[0].shape[2]), np.arange(data[0].shape[1]))
    contour = ax.contourf(x, y, z, 50, cmap='magma')
    cbar = fig.colorbar(contour, ax=ax)
    cbar.set_label('Intensity Scale')
    ax.set_title('Flux against position in protostellar source')
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Flux in Jy')  # Optional: Add a label for the z axis
    if SaveFig:
        plt.savefig(fname=SaveFigName, bbox_inches='tight', dpi=400)

    plt.show()
    plt.show()

    z[np.isnan(z)] = 0.0

    amplitude = np.nanmax(z)
    max_index = np.unravel_index(np.argmax(z), z.shape)
    x_max, y_max = max_index[1], max_index[0]

    gaussian = amplitude * np.exp(-((x - x_max) ** 2 / (2 * x_stddev ** 2)) - ((y - y_max) ** 2 / (2 * y_stddev ** 2)))

    image_data_subtracted = flux_density - gaussian #which is in jy/beam
    image_data_subtracted = (image_data_subtracted * u.Jy / u.beam).to(u.Jy / u.sr, equivalencies=u.beam_angular_area(beam_area_sr * u.sr)) #now in Jy/sr
    z_subtracted = image_data_subtracted * beam_area_sr *u.sr/(u.Jy) #now in Jy
    z_subtracted_2d = z_subtracted.reshape(x.shape)
    print(f'{np.nansum(z_subtracted):.3f} sum of Flux in Jy after subtraction')
    Mass_after_sub = mass_calc(z_subtracted*u.Jy)
    print(f'{np.nansum(Mass_after_sub):.3f} sum of sol mass after subtraction')

    #the difference between doing only above zero (3.5 sol mass instead of 7.3)
    z_no0 = z_subtracted[z_subtracted > 0]
    print(f'{np.nansum(z_no0):.3f} sum of Flux22 in Jy after subtraction')
    Mass_after_sub2 = mass_calc(z_no0*u.Jy)
    print(f'{np.nansum(Mass_after_sub2):.3f} sum of sol mass22 after subtraction')
    mask = z_subtracted_2d > 0
    z_no0_2d = np.where(mask, z_subtracted_2d, 0)


    er, err, x_gal, y_gal = galactic_coords(data, header)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    contour1 = ax1.contourf(x, y, z_no0_2d, 50, cmap='viridis') #replace z_subtracted_2d with z_no0_2d
    ax1.plot_surface(x, y, gaussian, color = 'red', alpha = 0.2)
    cbar1 = fig1.colorbar(contour1, ax=ax1)
    cbar1.set_label('Intensity Scale')
    ax1.set_xlabel('X Axis')
    ax1.set_ylabel('Y Axis')
    ax1.set_zlabel('Flux in Jy')
    ax1.set_title('Subtracted Brightness')

    if SaveFig:
        plt.savefig(fname=SaveFigName2, bbox_inches='tight', dpi=400)

    plt.show()

    return

def galactic_coords(data, header):

    #header = data[0].header

    #ctype1 = header['CTYPE1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    crpix1 = header['CRPIX1']

   # ctype2 = header['CTYPE2']
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

def mass_calc(Jy_image):
    mass_sol = 0.12 * (np.exp(1.439 * (lambda_mm)**(-1)*(temperature_K/10)**(-1)) - 1) * (Jy_image/(1*u.Jy)) *(kappa_nu/ 0.01)**(-1) * (distance/ 100)**2 * (lambda_mm)**3
    return mass_sol


def main():
    data, header = get_data(initial_dir)
    z = flux_calc(data, header)
    Jy_image_and_plot(data, header, z)

main()

