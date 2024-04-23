# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 14:08:05 2023

@author: Christopher
path: /mnt/c/users/christopher/onedrive/documents/4th_year_physics/MPhys_Project/s2/w4-6/ToomreQ/Toomre_Q_3.py

"""

import numpy as np
from astropy.io import fits
from tkinter import Tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker
#from mpl_toolkits.mplot3d import Axes3D as mpl
from matplotlib.lines import Line2D
#from scipy.constants import c
from astropy import units as u
from radio_beam import Beam
#import uncertainties.unumpy as unp
#from uncertainties import ufloat
#from astropy.convolution import equivalencies
from astropy.wcs import WCS
from reproject import reproject_interp
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize
#from astropy.coordinates import SkyCoord
from astropy.coordinates import SkyCoord, Galactic

G = 6.6743e-11
au = 1.495979e+11 #AU to m
M_o = 1.98847e+30
pc = 3.086e+16
#k = 1.3806488e-23
R = 8.3144621
#A_kep = 103.70690152932114 # km/s AU**0.5 ##nick only got 54.569 but its similar
#A_kep_SI = A_kep * 10**3 * au**0.5
#A_kep_SI_2 = 54.569 * 10**3 * au**0.5
A_kep_SI = 28544680664.453144
#bbarolo_path = 'C:/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/week_7-9/velocity_distance/bbarolo_data/densprof.txt'

def get_data():
    
    #continuum_path = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/new_data/continuum_small.fits'
    continuum_path = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/new_data/continuum_small.fits'
    #continuum_path = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/new_data/continuum_pacman.fits'
    #file_path = '/mnt/c/users/christopher/onedrive/documents/4th_year_physics/mphys_project/s2/w4-6/xclass/new_ch3_13cn_small/Temp_ch3_13cn_test.fits'
    #file_path = '/mnt/c/users/christopher/onedrive/documents/4th_year_physics/mphys_project/s2/w4-6/xclass/CH3-13CN_6it/T_rot_CH3-13CN_test_8.fits'
    file_path = '/mnt/c/users/christopher/onedrive/documents/4th_year_physics/mphys_project/s2/w4-6/xclass/CH3-13CN_6it/T_rot_small.fits'
    #file_path = '/mnt/c/users/christopher/onedrive/documents/4th_year_physics/mphys_project/s2/w4-6/xclass/CH3-13CN_6it/BestResult___parameter__T_rot___molecule__CH3C-13-N_v=0____component__1___Iter__6.fits'
    #file_path = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/Mphys_Project/Data-28-09-23/T_rot_CH3C-13-N.fits'
    continuum = fits.open(continuum_path)
    #Tk().withdraw()
    #file_paths = filedialog.askopenfilenames()
    #for file_path in file_paths:
    data_temp = fits.open(file_path)
    #continuum = fits.open(continuum_path)
    filename = (os.path.basename(file_path)).replace(".fits", "")#.replace(".image.pbcor_line.galactic","")
    main(continuum, data_temp, filename, file_path)
    data_temp.close()
    #continuum = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/Mphys_Project/Data-28-09-23/Continuum_Toomre_Q.fits'
    continuum.close()
    return

def main(data,data_temp,name,file_path):

    print(data[0].shape)
    print(data_temp[0].shape)
    #print(data[0].data[0])
    # Assuming the relevant data is in the last two dimensions of the continuum image
    continuum_data_2d = data[0].data[0, 0, :, :]
    
    # Extract WCS for the 2D slice of the continuum
    wcs_continuum_2d = WCS(data[0].header, naxis=2)
    
    # Reproject temperature map to match continuum image's WCS and pixel scale
    data_temp, footprint = reproject_interp(data_temp[0], wcs_continuum_2d, shape_out=continuum_data_2d.shape)
    
    # Create a new header for the reprojected temperature map
    new_header = data[0].header.copy()
    new_header['NAXIS'] = 2  # Update the number of axes
    for key in ['NAXIS3', 'NAXIS4', 'CRPIX3', 'CRPIX4', 'CDELT3', 'CDELT4', 'CRVAL3', 'CRVAL4', 'CTYPE3', 'CTYPE4']:
        if key in new_header:
            del new_header[key]

    # Save the reprojected temperature map
    data_temp = fits.PrimaryHDU(data=data_temp, header=new_header)
    #print(data_temp.data)
    
    data[0].data = match_shapes_with_nan(data[0].data)
    #print(data_temp.data)
    data_temp.data = match_shapes_with_nan(data_temp.data)
    print(data[0].shape)
    
    m = data[0].shape[2] # y shape
    p = data[0].shape[3] # x shape

    header = data[0].header
    header_temp = data_temp.header
    #print(header)
    ref_freq = header['CRVAL3'] #reference frequency
    incr_freq = header['CDELT3'] #frequency increment
    pix_1 = np.deg2rad(header['CDELT1'])
    pix_2 = np.deg2rad(header['CDELT2'] )
    pixel_size = pix_1*pix_2
    print(pix_1,pix_2)
    print(pixel_size)
    pixel_area = 0.5 * np.abs(pixel_size**1) * (8.1 * 10**3 * pc)**2
    beam_major = np.deg2rad(header['BMAJ'])  # Major axis of the beam
    beam_minor = np.deg2rad(header['BMIN'])  # Minor axis of the beam
    
    pix_1_temp = np.deg2rad(header_temp['CDELT1'])
    pix_2_temp = np.deg2rad(header_temp['CDELT2'] )
    pixel_size_temp = pix_1_temp*pix_2_temp
    print(pix_1_temp,pix_2_temp)
    print(pixel_size_temp)
    #pixel_area = 0.5 * np.abs(pixel_size**1) * (8.1 * 10**3 * pc)**2
    #beam_major = np.deg2rad(header['BMAJ'])  # Major axis of the beam
    #beam_minor = np.deg2rad(header['BMIN'])  # Minor axis of the beam
    
    beam_major_arcsec = header['BMAJ'] * 3600
    beam_minor_arcsec = header['BMIN'] * 3600
    print(f"FWHM Major Axis: {beam_major_arcsec} arcsec")
    print(f"FWHM Minor Axis: {beam_minor_arcsec} arcsec")
    print(np.sqrt(beam_major_arcsec*beam_minor_arcsec))
    
    beam_sizer = (np.pi/(4 *np.log(2))) * beam_major * beam_minor
    print(f'{beam_minor:.3e} {beam_major:.3e} {beam_sizer:.3e} minor, major and size in rad, rad, sr')
    beam_size = beam_sizer * u.sr
    data[0].data[0] = (data[0].data[0]*u.Jy/u.beam).to(u.Jy/u.sr, equivalencies=u.beam_angular_area(beam_size))
    F_v = data[0].data[0] * beam_sizer #* np.cos(np.deg2rad(38))
    print(f'{np.nanmax(F_v):.3f} max flux in Jansky')
    print(f'{np.nansum(F_v):.3f} sum flux in Jansky')
    print(f'{m*p} No of points')
    #temps = temp_calc()
    wave = 1.36
    
    #data_temp[0].data = data_temp[0].data / np.cos(np.deg2rad(35))

    #xo_pix,yo_pix = find_max_value_coordinates(data[0].data[0][0])
    #xo_pix,yo_pix = 62, 62
    xo_pix,yo_pix = 17,17
    #xo_pix,yo_pix = 24,25
    print('centre:')
    print(xo_pix,yo_pix)
    
    x_m, y_m, x_gal, y_gal = galactic_coords(data[0], header)
    xo_m = x_m[xo_pix] #central xo point
    yo_m = y_m[yo_pix] #central yo point
    x_norm = x_m - xo_m #normalising around the xo point
    y_norm = y_m - yo_m #normalising around the yo point
    x_trans,y_trans = coord_trans(x_norm, y_norm)
    #x_trans, y_trans = coorderz_trans_3(x_gal,y_gal)
    
    y_temp = []
    for y_value in y_trans:
        if -1100 < y_value/au < 1100: #2400 #1200
            y_temp.append(y_value)
        else:
            y_temp.append(np.nan)
    
    y_trans = np.vstack(y_temp)
    
    x_temp = []
    for x_value in x_trans:
        if -1100 < x_value/au < 1100: #2600 #1200
            x_temp.append(x_value)
        else:
            x_temp.append(np.nan)
    
    x_trans = np.vstack(x_temp)
    
    xx, yy = np.meshgrid(x_trans, y_trans)
    ang_vel = np.zeros((m, p))
    M = np.zeros((m, p))
    part1 = np.zeros((m, p))
    part3 = np.zeros((m, p))
    #temps = np.zeros((m, p))
    M_star = np.zeros((m, p))
    
    for i in range(m): # i is y
        for j in range(p): #  j is x
            #print(F_v[0,i,j])
            #temps[i,j] = temp_calc(i,j)
            ang_vel[i,j] = ang_vel_calc(x_trans[j], y_trans[i])
            M[i,j], part1[i,j],part3[i,j] = np.abs(mass_calc(wave, data_temp.data[i,j], F_v[0,i,j], xx[i,j],yy[i,j]))
            #M_star[i,j] = mass_calc()
            
    print(M)
    print(f'{np.nanmax(M):.3f} max sol mass')
    print(f'{np.nanmin(M):.3f} min sol mass')
    print(f'{np.nansum(M):.3f} sum sol mass for r < 2000AU')

    surf_density = (M * M_o)/pixel_area
    print(f'{pixel_area:.3g} pixel area in m^2')
    print('Max surf density = '+ str(np.nanmax(surf_density/M_o)))
    
    ang_vel_flat = ang_vel.flatten()
    r_flat = (np.sqrt(xx**2+yy**2)).flatten()
    ang_dist = np.stack((r_flat,ang_vel_flat)).T
    ang_dist_trimmed = trimming(ang_dist)
    ang_vel_av = np.nanmean(ang_dist_trimmed[:,1])
    ang_vel_med = np.nanmedian(ang_dist_trimmed[:,1])
    
    print(str(ang_vel_av)+' = mean ang vel')
    print(str(ang_vel_med)+' = median ang vel')    
    cs = speed_sound(data_temp.data)
    print('cs max = ' +str(np.nanmax(cs)))
    print('cs min = ' +str(np.nanmin(cs)))
    
    Q = Toomre_Q(cs, ang_vel, surf_density)
    Q_ss = Toomre_Q(cs,ang_vel_av,surf_density)
    
    print('Q max = ' +str(np.nanmax(Q)))
    print('Q min = ' +str(np.nanmin(Q)))
    print('Q_ss max = ' +str(np.nanmax(Q_ss)))
    print('Q_ss min = ' +str(np.nanmin(Q_ss)))
    #plotting(Q,x_gal,y_gal)
    # Create a contour plot with Q as the contour levels
    print(str(np.nanmax(data_temp.data))+', '+str(np.nanmin(data_temp.data)) + ' max,min temperature')
    # Stack Q, x_gal, and y_gal together to form a 2x2 matrix
    #print(xx/au)
    valid_indices = (-1100 < xx/au) & (xx/au < 1100) #2600
    valid_indices_2 = (-1100 < yy/au) & (yy/au < 1100) #2400
    #valid_indices = (-1200 < xx/au) & (xx/au < 1200)
    #valid_indices_2 = (-1200 < yy/au) & (yy/au < 1200)
    valid_indices_combined = valid_indices & valid_indices_2
    xx[~valid_indices] = 0
    yy[~valid_indices_2] = 0
    
    
    #xx = np.ma.masked_where(~valid_indices, xx)
    #yy = np.ma.masked_where(~valid_indices_2, yy)
    print(np.nanmax(xx/au))
    print(np.nanmin(xx/au))
    
    print(np.nanmax(ang_vel),np.nanmin(ang_vel))
    
    #Q = np.where(Q < 100, -1, Q)
    combined_matrix = np.stack(((Q), xx, yy), axis=-1)  
    combined_q_ss = np.stack((Q_ss,xx,yy), axis=-1)
    #print(combined_matrix)
    #plotting(combined_matrix, 0)
    plotting_new(combined_matrix, 0)
    """
    fig, ax = plt.subplots()
    # Use np.ma.masked_invalid to hide the invalid values in the plot
    cmap = plt.cm.viridis  # Choose a colormap that suits your data
    cmap.set_bad(color='black')  # Set color for masked elements
    im = ax.imshow(np.ma.masked_invalid(data_temp.data),cmap=cmap)
    ax.set_title('')
    plt.colorbar(im)  # Show the color bar
    plt.show()
    """
    
    #main_2(data,header,x_norm,y_norm,F_v,data_temp)
    return

def main_2(data,header,x_norm,y_norm,z,data_temp):
    
    x_m, y_m, x_gal, y_gal = galactic_coords(data[0], header)
    #x,y = x*
    pixel_scale_lon = abs(header['CDELT1'])
    pixel_scale_lat = abs(header['CDELT2'])
    beam_major_pixels = header['BMAJ'] / pixel_scale_lon 
    beam_minor_pixels = header['BMIN'] / pixel_scale_lat 
    
    beam_major = np.deg2rad(header['BMAJ']) *8.1 * 10**3 * pc
    beam_minor = np.deg2rad(header['BMIN']) *8.1 * 10**3 * pc
    
    fwhm_stddev_maj_y = beam_major / (2 * np.sqrt(2 * np.log(2))) #*8.1 * 10**3 * pc
    fwhm_stddev_min_x = beam_minor / (2 * np.sqrt(2 * np.log(2))) #*8.1 * 10**3 * pc
    
    amplitude = np.nanmax(z)
    x_stddev =  fwhm_stddev_min_x # Convert FWHM to standard deviation
    y_stddev = fwhm_stddev_maj_y
    #x_max, y_max = 24, 24
    x_mid = x_m[24]
    y_mid = y_m[25]
    #x_mid = x_m[62]
    #y_mid = y_m[62]
    x = x_m - x_mid
    y = y_m - y_mid
    x,y=np.meshgrid(x,y)
    gaussian = amplitude * np.exp(-((x - 0) ** 2 / (2 * x_stddev ** 2)) - ((y - 0) ** 2 / (2 * y_stddev ** 2)))
    print(gaussian.shape)
    print(z.shape)
    z = np.squeeze(z) - np.abs(gaussian)
    
    m = data[0].shape[2] # y shape
    p = data[0].shape[3] # x shape
    ang_vel = np.zeros((m, p))
    M = np.zeros((m, p))
    part1 = np.zeros((m, p))
    part3 = np.zeros((m, p))
    #temps = np.zeros((m, p))
    M_star = np.zeros((m, p))
    wave = 1.36
    
    for i in range(m): # i is y
        for j in range(p): #  j is x
            #print(F_v[0,i,j])
            #temps[i,j] = temp_calc(i,j)
            #ang_vel[i,j] = ang_vel_calc(x_trans[j], y_trans[i])
            M[i,j], part1[i,j],part3[i,j] = mass_calc(wave, data_temp.data[i,j], z[i,j], x[i,j],y[i,j])
    
    
    print(str(np.nansum(np.abs(M)))+' = MASSSSSS')
    
    fig=plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    epicalnations = ax.contourf(x,y,z,cmap='magma')
    fat = ax.plot_surface(x,y,gaussian)
    plt.colorbar(epicalnations)
    plt.show()
    
    return
           
           

def plotting(matrix, numba):
    """
    #epicals = plt.contourf(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], levels = 50)
    no2=20
    colors = plt.cm.jet(np.linspace(0, 1, round(no2/2)))
    
    # Create a new colormap from the list of colors
    cmap = ListedColormap(colors)
    #no2=22
    # Define the boundaries for each color
    bounds = np.linspace(0, no2, round(no2/2))  # 36 intervals, 37 boundaries
    
    # Create a 'Normalize' object which maps the data values to the interval [0, 1]
    norm = plt.Normalize(min(bounds), max(bounds))
    
    fig = plt.figure(figsize=(12,6))
    epic = plt.contourf(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], levels = 60, cmap=cmap, norm=norm)
    contour_levels = np.arange(0, 10, 2)
    contour_levels_2 = np.arange(0, 9, 2)

    # Create contour lines up to level 30
    epic2 = plt.contour(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], levels=contour_levels, cmap='Purples', alpha=0.6)
    #epic3 = plt.contour(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], levels=contour_levels_2, cmap='Greys', alpha=0.6)
    # plt.clabel will now only label contours up to the maximum level specified, which is 30
    plt.clabel(epic2, inline=True, fontsize=8)
    #plt.clabel(epic3, inline=True, fontsize=8, colors='white')
    """
    #numba=0
    
    if numba == 0:
        
        no2 = 30 # Set to your maximum value of interest for the colorbar
        over_value_color = 'red'  # Color to use for values over no2
        
        # Generate a colormap with no2/2 colors from the 'jet' colormap
        colors = plt.cm.jet(np.linspace(0, 1, no2//2))  # Ensure we have a color for each even number up to no2
        colors_2 = plt.cm.jet(np.linspace(1, 1, 2))
        # Create a new colormap from the list of colors
        cmap = ListedColormap(colors)
        cmap_2 = ListedColormap(colors_2)
        
        # Define the boundaries for each color
        bounds = np.linspace(0, no2, cmap.N+1)
        bounds_2 = np.linspace(28, 72, cmap_2.N+1)
        
        # Create a 'Normalize' object which maps the data values to the interval [0, 1]
        norm = Normalize(vmin=0, vmax=no2)
        norm_2 = Normalize(vmin=28, vmax=1000) #72
        norm_3 = Normalize(vmin=30, vmax=1000)
        contour_levels = np.arange(0, 10, 2)
        
        # Create the plot
        fig = plt.figure(figsize=(7.5,5.5)) #7.5,5,5 #6.5,5,5 #0.34
        epic_2 = plt.contourf(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], norm=norm_3, cmap=cmap_2)
        epic = plt.contourf(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], levels=bounds, cmap=cmap, norm=norm)
        
        ######################
        #epic_2 = plt.imshow(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], vmin = 28, cmap=cmap_2)
        #epic = plt.imshow(matrix[:,:,0], cmap='jet', interpolation='nearest')
        #####################
        
        #epic_2 = plt.contourf(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], norm=norm_3, cmap=cmap_2)
        epic2 = plt.contour(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], levels=contour_levels, cmap='Purples', alpha=0.8)
        # Create the colorbar with ticks up to no2
        cbar = plt.colorbar(epic, ticks=np.arange(0, no2, 2))  # Ticks every 2 units
        
        # Set the tick labels manually, add '+' to the last label
        tick_labels = [str(tick) for tick in range(0, no2-2, 2)] + ['28+']
        cbar.set_ticklabels(tick_labels, size=10)
        
        # Optionally, change the color of the last tick label
        cbar.ax.get_yticklabels()[-1].set_color(over_value_color)
        cbar.ax.get_yticklabels()[-1].set_weight('bold')  # Optionally set the last label to bold
        #plt.clabel(epic2, inline=True, fontsize=10)
        """
        fig = plt.figure(figsize=(12,6))
        pig = plt.contourf(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], cmap='jet')
        pig_2 = plt.contour(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], cmap='Greys')
        plt.colorbar(pig)
        plt.clabel(pig_2, inline=True, fontsize=8)
        """
        plt.xlabel("AU", size=12)
        plt.ylabel("AU", size=12)
        plt.title('Toomre Q - Keplerian - CH$_{3}^{\:13}$CN', size=13)
        #plt.savefig('/mnt/c/users/christopher/onedrive/documents/4th_year_physics/MPhys_Project/Week_7-9/Toomre_Q_map_Keplarian_MAIN_lower_newcoord.png',bbox_inches='tight', dpi=800)
        #plt.savefig('/mnt/c/users/christopher/onedrive/documents/4th_year_physics/MPhys_Project/s2/w4-6/ToomreQ/Toomre_Q_CH3-13CN_small.pdf',bbox_inches='tight')
    elif numba == 1:
        fig = plt.figure(figsize=(12,6))
        pig = plt.contourf(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], cmap='jet')
        pig_2 = plt.contour(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], cmap='Greys')
        plt.colorbar(pig)
        plt.clabel(pig_2, inline=True, fontsize=8)
        plt.title('Toomre Q - Sound Speed only - CH$_{3}$CN')
        plt.xlabel("AU")
        plt.ylabel("AU")
        #plt.savefig('/mnt/c/users/christopher/onedrive/documents/4th_year_physics/MPhys_Project/s2/w4-6/ToomreQ/Toomre_Q_CH3-13CN.pdf',bbox_inches='tight')
    plt.show()
    """
    # Add labels, titles, colorbar, etc. as needed
    plt.xlabel("AU")
    #plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(nbins=6)) #plots 6 points along x axis
    plt.ylabel("AU")
    #plt.invert_xaxis()
    plt.title('Contour Map of Q')
    plt.colorbar(epic, label='Q Contour Levels')
    #path = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/Mphys_Project/'
    plt.savefig('Toomre_Q_map.png',bbox_inches='tight')

    # Show the plot
    plt.show()
    """

def plotting_new(matrix, numba):
    
    fig = plt.figure(figsize=(7.5,6.5))
    epic = plt.contourf(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0]*0.34, cmap='jet',
                        vmax=5, levels = 100)
    cbar = plt.colorbar(epic)
    
    plt.xlabel("AU", size=12)
    plt.ylabel("AU", size=12)
    plt.title('Toomre Q - Keplerian - CH$_{3}$CN', size=13)
    plt.show()
    
    return

def mass_calc(wave,T,F_v,x,y):
    """
    wavelength : in mm
    T : in K
    k_v : 1.99 cm^2 g^-1
    F_v : From contintuum in Jy
    d : in pc, (8.1 kpc)
    M : solar masses
    """
    d = 8.1*10**3
    g=100
    k_v = 1.99
    
    r = np.sqrt((x**2+y**2)) / au
    
    part1 = (np.exp(1.439 * (1.36*(T/10))**(-1) ) -1 )
    
    part2 = ((k_v/0.01)**(-1))
    part3 = (F_v)*((d/100)**2)
    part4 = ((wave)**3)
    
    if not r < 50:
        M = 0.12*part1*part2*part3*part4 #* g
    else:
        M = np.nan
        
    
    return M, part1, part3

def temp_calc(i,j):
    """
    r = np.sqrt((x**2+y**2)) / au
    
    if r > 0:
        temp = -0.04*r + 270
    else:
        temp = np.nan
    """
    #temperature = temp_array[i,j]
        
    
    #temp = 100
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
    x_gal = (crval1 + (np.arange(data.shape[3]) - crpix1) * cdelt1)
    y_gal = (crval2 + (np.arange(data.shape[2]) - crpix2) * cdelt2)

    x_m = x_gal * 8.1 * 10**3 * 206264.81 * np.pi/180 * au
    y_m = y_gal * 8.1 * 10**3 * 206264.81 * np.pi/180 * au

    return x_m,y_m, x_gal, y_gal

def Toomre_Q(cs,ang_vel,E):

    Q = cs*ang_vel / (G*E*np.pi)

    return Q

def ang_vel_calc(x,y):

    r = np.sqrt(x**2 + y**2)
    """
    if r ==0:
        ang_vel = np.nan
    else:
    """
    try:
        v_rot = A_kep_SI * r **(-0.5)
        ang_vel = v_rot / r
        #print(v_rot)
    except RuntimeWarning:
        ang_vel = np.nan

    return ang_vel

def speed_sound(T):

    Molar = 2.8 * 10e-3
    gam = 1.4
    cs = np.sqrt(gam*R*T/Molar)

    return cs

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

def find_max_value_coordinates(arr):
    max_value = np.nanmax(arr)
    max_x, max_y = np.where(arr == max_value)
    return max_x[0], max_y[0]

def coord_trans(x,y):
    phi = np.deg2rad(170)
    inc = np.deg2rad(38)
    #z=
    alpha = phi + np.pi
    
    x_prime = x*np.cos(alpha) + y*np.sin(alpha)
    y_prime = x*np.sin(alpha)/np.cos(inc) + y*np.cos(alpha)/np.cos(inc)
    """
    inclination_angle = 157  # Inclination angle in degrees
    position_angle = 34  # Position angle in degrees
S
    # Convert to SkyCoord object
    coords = SkyCoord(x * u.pc, y * u.pc, z * u.pc, frame='galactic', representation_type='cartesian')
    
    # Apply rotation around the z-axis by the position angle
    rotated_coords = coords.transform_to(coords.skyoffset_frame(rotation=position_angle * u.deg))
    
    # Apply inclination by rotating around the y-axis
    inclined_coords = rotated_coords.transform_to(rotated_coords.skyoffset_frame(rotation=inclination_angle * u.deg))
    
    # Extract the primed coordinates
    x_prime = inclined_coords.cartesian.x.value
    y_prime = inclined_coords.cartesian.y.value
    
    
    transformation_matrix = np.array([
    [np.cos(phi), -np.sin(phi) * np.cos(inc)],
    [np.sin(phi), np.cos(phi) * np.cos(inc)]
    ])

    # Apply the transformation matrix to the coordinates
    xy_prime = np.dot(transformation_matrix, np.array([x, y]))
    print(np.vstack((xy_prime)))
    x_prime,y_prime = xy_prime
    """
    return x_prime,y_prime
   # return x,y
    
def trimming(arr):
    temp=[]
    for line in arr:
        if 10 < line[0]/au < 2100:
            temp.append(line)
    
    return np.vstack((temp))

def coorderz_trans_2(x_gal,y_gal):
    
    
    l = y_gal  # Galactic longitude in degrees
    b = x_gal  # Galactic latitude in degrees
    
    # Convert galactic coordinates to RA/DEC
    galactic_coords = SkyCoord(l=l*u.degree, b=b*u.degree, frame=Galactic)
    ra_dec_coords = galactic_coords.icrs
    
    # Convert RA/DEC to physical units (x, y in meters)
    # Assuming distance to the source is 8.1 kpc
    distance = 8.1 * u.kpc
    x = (ra_dec_coords.ra.deg * distance).to(u.meter, u.dimensionless_angles())
    y = (ra_dec_coords.dec.deg * distance).to(u.meter, u.dimensionless_angles())
    
    
    x = x-x[24]
    y = y-y[25]
    # Apply inclination and position angle transformations
    # Replace these with your actual angles
    inclination_angle = 38 # in degrees
    position_angle = 170   # in degrees
    
    # Convert angles to radians for calculations
    inclination_rad = np.radians(inclination_angle)
    position_rad = np.radians(position_angle)
    
    # Transformation matrix
    transformation_matrix = np.array([
    [np.cos(position_rad), np.sin(position_rad)],
    [-np.sin(position_rad) / np.cos(inclination_rad), np.cos(position_rad) * np.cos(inclination_rad)]
    ])
    
    # Applying the transformation
    xy_coords = np.vstack((x.value, y.value))
    transformed_coords = np.dot(transformation_matrix, xy_coords)
    
    # Extracting transformed x and y
    transformed_x, transformed_y = transformed_coords
    
    # Print or process the transformed coordinates
    print(transformed_x, transformed_y)
    

    return transformed_x, transformed_y


def coorderz_trans_3(x_gal, y_gal):
    # Galactic longitude (l) and latitude (b)
    l = y_gal  # in degrees
    b = x_gal  # in degrees
    
    # Convert galactic coordinates to physical units (x, y in meters)
    # Assuming distance to the source is 8.1 kpc
    distance = 8.1 * u.kpc
    l_rad = np.radians(l)
    b_rad = np.radians(b)

    # Calculating x, y positions
    # This assumes a simple projection, which may not be accurate for all purposes
    x = (distance * np.cos(b_rad) * np.sin(l_rad)).to(u.meter)
    y = (distance * np.sin(b_rad)).to(u.meter)

    # Centering the coordinates around a specific point (assuming [24] is the center)
    x = x - x[24]
    y = y - y[25]
  
    # Inclination and position angle in degrees
    inclination_angle = 38
    position_angle = 170
    
    # Convert angles to radians
    inclination_rad = np.radians(inclination_angle)
    position_rad = np.radians(position_angle)
    #position_rad = np.pi - position_rad
    
    # Transformation matrix
    transformation_matrix = np.array([
        [np.cos(position_rad), -np.sin(position_rad)],
        [np.sin(position_rad) / np.cos(inclination_rad), np.cos(position_rad) / np.cos(inclination_rad)]
    ])
    
    # Applying the transformation
    xy_coords = np.vstack((x.value, y.value))
    transformed_coords = np.dot(transformation_matrix, xy_coords)
    
    # Extracting transformed x and y
    transformed_x, transformed_y = transformed_coords
    
    # Print or process the transformed coordinates
    print(transformed_x, transformed_y)

    return transformed_x, transformed_y


get_data()
