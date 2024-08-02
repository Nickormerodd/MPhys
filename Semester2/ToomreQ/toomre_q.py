# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 14:08:05 2023
Edited June 2024

@author: Christopher Andrews
path: /mnt/c/users/christopher/onedrive/documents/4th_year_physics/MPhys_Project/s2/w4-6/ToomreQ/Toomre_Q_3.py

Input the continuum and temperature filepath and adjust the user-defined inputs
at the top. The temperature map is reprojected to match the shape of the 
continuum image, as the temp map has more pixels over the same area.
This code assumes Keplerian motion and needs a constant of proportionality
which can be taken from 3DBarolo.
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.wcs import WCS
from reproject import reproject_interp
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize
from astropy.coordinates import SkyCoord, Galactic
from matplotlib.patches import Ellipse
import copy

G = 6.6743e-11 # grav const
au = 1.495979e+11 # AU to m
M_o = 1.98847e+30 # solar mass
pc = 3.086e+16 # pc to m
R = 8.3144621 # molar const


##############################################################################
# User defined inputs

XO_PIX,YO_PIX = 25,25 # Define central pixel to normalise around, 
# or set as 'true' for max continuum flux to be the cental pixel point:
MAX_PIX_FLUX = False

# The constant in v = A * r **(-0.5) in SI units (take from 3DBarolo or velocity_curve_err.py)
A_kep_SI = 28544680664.453144

# x and y length (in AU) to be cut on either side from centre (it has to be defined)
X_CUT = 1000 # AU
Y_CUT = 1000 # AU
# This is also used to calcualte the disc mass within this distance scale

# Define a radial distance to calcualte the Q > 2 disc mass within
R_DISC = 800 # AU

# Compute coord trans into plane of disc? true/false
TRANSFORMATION = False

# Distance to source in kpc
DISTANCE_TO_SOURCE = 8.1

##############################################################################

def get_data():
    
    # Continuum data is in Jy/beam and is 4-dimensions 
    # Temperature data is in Kelvin and is 2 dimensions
    continuum_path = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/new_data/continuum_Xreg.fits'
    temperature_path = '/mnt/c/users/christopher/onedrive/documents/4th_year_physics/mphys_project/s2/w4-6/xclass/CH3-13CN_6it/T_rot_CH3-13CN_test_8.fits'
    continuum = fits.open(continuum_path)

    data_temp = fits.open(temperature_path)
    main(continuum, data_temp)
    data_temp.close()
    continuum.close()
    return

def main(data_contin, data_temp):
    
    # Assuming the relevant data is in the last two dimensions of the continuum image
    continuum_data_2d = data_contin[0].data[0, 0, :, :]
    
    # Extract WCS for the 2D slice of the continuum
    wcs_continuum_2d = WCS(data_contin[0].header, naxis=2)
    
    # Reproject temperature map to match continuum image's WCS and pixel scale
    data_temp, footprint = reproject_interp(data_temp[0], wcs_continuum_2d, shape_out=continuum_data_2d.shape)
    
    # Create a new header for the reprojected temperature map
    new_header = data_contin[0].header.copy()
    new_header['NAXIS'] = 2  # Update the number of axes
    for key in ['NAXIS3', 'NAXIS4', 'CRPIX3', 'CRPIX4', 'CDELT3', 'CDELT4', 'CRVAL3', 'CRVAL4', 'CTYPE3', 'CTYPE4']:
        if key in new_header:
            del new_header[key]

    # Save the reprojected temperature map
    data_temp = fits.PrimaryHDU(data=data_temp, header=new_header)
    
    # Matching shapes with NaN so can perfom coord trans
    data_contin[0].data = match_shapes_with_nan(data_contin[0].data)
    data_temp.data = match_shapes_with_nan(data_temp.data)
    
    m = data_contin[0].shape[2] # y shape
    p = data_contin[0].shape[3] # x shape

    header = data_contin[0].header
 
    # Pixel area from the continuum header
    pix_1 = np.deg2rad(header['CDELT1'])
    pix_2 = np.deg2rad(header['CDELT2'] )
    pixel_size = pix_1*pix_2
    pixel_area = 0.5 * np.abs(pixel_size**1) * (8.1 * 10**3 * pc)**2
    
    F_v = jy_per_beam_to_jy(data_contin)
    print(f'Sum flux in Jansky: {np.nansum(F_v):.3f}')

    # Either user defined central pixels or max continuum peak as centre
    if MAX_PIX_FLUX:
        yo_pix,xo_pix = find_max_value_coordinates(data_contin[0].data[0][0])
    else:
        xo_pix, yo_pix = XO_PIX, YO_PIX

    print('centre:')
    print(xo_pix,yo_pix)
    
    x_m, y_m, x_gal, y_gal = galactic_coords(data_contin[0], header)
    xo_m = x_m[xo_pix] #central xo point
    yo_m = y_m[yo_pix] #central yo point
    x_norm = x_m - xo_m #normalising around the xo point
    y_norm = y_m - yo_m #normalising around the yo point
    
    # Compute coord transformation
    if TRANSFORMATION:
        # F_v = F_v * np.sin(np.deg2rad(38)) # Not sure whether this is needed
        # We observe larger flux per pixel. When transformed, observe disc over 
        # larger distance, so flux per pixel changes as inclination angle ?
        
        x_trans, y_trans = coord_trans(x_norm, y_norm)
    else:
        x_trans, y_trans = x_norm, y_norm
    
    # Trimming
    y_temp = []
    for y_value in y_trans:
        if -Y_CUT < y_value/au < Y_CUT:
            y_temp.append(y_value)
        else:
            y_temp.append(np.nan)
    y_trans = np.vstack(y_temp)
    
    x_temp = []
    for x_value in x_trans:
        if -X_CUT < x_value/au < X_CUT:
            x_temp.append(x_value)
        else:
            x_temp.append(np.nan)
    x_trans = np.vstack(x_temp)
    
    xx, yy = np.meshgrid(x_trans, y_trans)
    
    ang_vel = np.zeros((m, p))
    M = np.zeros((m, p))
    
    # Calculating the mass
    wave = 1.36
    for i in range(m): # i is y
        for j in range(p): #  j is x
            ang_vel[i,j] = ang_vel_calc(x_trans[j], y_trans[i])
            M[i,j] = np.abs(mass_calc(wave, data_temp.data[i,j], F_v[0,i,j], xx[i,j],yy[i,j]))

    # Total mass
    print(f'\n{np.nansum(M):.3f} sum sol mass\n')

    surf_density = (M * M_o)/pixel_area
    print(f'{pixel_area:.3g} pixel area in m^2')
    print('Total surf density = '+ str(np.nansum(surf_density/M_o)))
    print('Max surf density = '+ str(np.nanmax(surf_density/M_o)))
    
    # Radial distance
    radial_dist = np.sqrt(xx**2 + yy**2)
    radial_dist[radial_dist == 0] = np.nan
    radial_valid = (radial_dist/au < R_DISC)
    
    # Trimming
    valid_indices = (-X_CUT < xx/au) & (xx/au < X_CUT)
    valid_indices_2 = (-Y_CUT < yy/au) & (yy/au < Y_CUT)
    xx[~valid_indices] = 0
    yy[~valid_indices_2] = 0
    
    ang_vel_flat = ang_vel.flatten()
    r_flat = (np.sqrt(xx**2+yy**2)).flatten()
    ang_dist = np.stack((r_flat,ang_vel_flat)).T
    ang_dist_trimmed = trimming(ang_dist)
    ang_vel_av = np.nanmean(ang_dist_trimmed[:,1])
    ang_vel_med = np.nanmedian(ang_dist_trimmed[:,1])
    
    print(str(ang_vel_av)+' = mean ang vel')
    print(str(ang_vel_med)+' = median ang vel')    
    cs = speed_sound(data_temp.data)
    
    # Calculating the mass within area of interest
    M_mask = copy.deepcopy(M)
    M_mask[~valid_indices] = np.nan
    M_mask[~valid_indices_2] = np.nan
    print(f'\n{np.nansum(M_mask):.3f} sum sol mass for x < {X_CUT}, y < {Y_CUT} AU\n')
    
    # Toomre Q
    Q = Toomre_Q(cs, ang_vel, surf_density)

    # Total mass for Q > 2 in the area of interest
    mask_Q_2 = (Q < 2)
    M_2 = copy.deepcopy(M)
    M_2[mask_Q_2] = np.nan
    M_2[~valid_indices] = np.nan
    M_2[~valid_indices_2] = np.nan
    total_mass = np.nansum(M_2) # total mass for Q > 2 in the area of interest
    print(f'{total_mass:.3f} sum sol mass for Q > 2 (disc) and x < {X_CUT}, y < {Y_CUT} AU')
    M_mask_radial = copy.deepcopy(M)
    M_mask_radial[mask_Q_2] = np.nan
    M_mask_radial[~radial_valid] = np.nan
    total_radial_mass = np.nansum(M_mask_radial) # total mass Q > 2 + r < R_DISC
    print(f'{total_radial_mass:.3f} sum sol mass for Q > 2 (disc) and r < {R_DISC} AU')
    
    print('\nQ max = ' +str(np.nanmax(Q)))
    print('Q min = ' +str(np.nanmin(Q)))

    combined_matrix = np.stack(((Q), xx, yy), axis=-1)
    plotting(combined_matrix, header)
 
    return

def jy_per_beam_to_jy(hdul):
    """
    Convert 2D FITS data from Jy/beam to Jy using beam information from the FITS header.
    
    Parameters:
    - fits_file_path (str): Path to the input FITS file.
    
    Returns:
    - 2D numpy array: Converted data in Jy.
    """
    # Open the FITS file and extract data and header
    header = hdul[0].header
    data = hdul[0].data[0]
    # Extract beam information from header
    beam_major = header['BMAJ'] * u.deg
    beam_minor = header['BMIN'] * u.deg
    pixel_scale = np.abs(header['CDELT2']) * u.deg
    
    # Convert beam major and minor axis to radians for beam area calculation
    beam_major_rad = beam_major.to(u.rad)
    beam_minor_rad = beam_minor.to(u.rad)
    
    # Calculate the beam area in steradians
    fwhm_to_sigma = 1. / (8 * np.log(2))**0.5
    beam_sigma_major = beam_major_rad * fwhm_to_sigma
    beam_sigma_minor = beam_minor_rad * fwhm_to_sigma
    beam_area_sr = 2 * np.pi * beam_sigma_major * beam_sigma_minor
    
    # Convert Jy/beam to Jy/sr
    data_jy_per_sr = (data * u.Jy / u.beam).to(u.Jy / u.sr, equivalencies=u.beam_angular_area(beam_area_sr))
    
    # Convert pixel scale to steradians
    pixel_scale_sr = (pixel_scale.to(u.rad))**2
    
    # Multiply Jy/sr by the solid angle per pixel to get Jy
    jy_data = data_jy_per_sr * pixel_scale_sr
    
    return jy_data.value

           
def plotting(matrix, header_temp):
        
    no2 = 30 # Set to your maximum value of interest for the colourbar
    over_value_color = 'red'  # Colour to use for values over no2
        
    # Generate a colourmap with no2/2 colours from the 'jet' colourmap
    colors = plt.cm.jet(np.linspace(0, 1, no2//2))
    colors_2 = plt.cm.jet(np.linspace(1, 1, 2))
    # Create a new colourmap from the list of colours
    cmap = ListedColormap(colors)
    cmap_2 = ListedColormap(colors_2)
        
    # Define the boundaries for each colour
    bounds = np.linspace(0, no2, cmap.N+1)
        
    # Create a 'Normalize' object which maps the data values to the interval [0, 1]
    norm = Normalize(vmin=0, vmax=no2)
    norm_3 = Normalize(vmin=30, vmax=1000)
    contour_levels = np.arange(0, 10, 2)
        
    # Create the plot
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    epic_2 = plt.contourf(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], norm=norm_3, cmap=cmap_2)
    epic = plt.contourf(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], levels=bounds, cmap=cmap, norm=norm)
        
    epic2 = plt.contour(matrix[:,:,1]/au, matrix[:,:,2]/au, matrix[:,:,0], levels=contour_levels, cmap='Purples', alpha=0.8)
    # Create the colourbar with ticks up to no2
    cbar = plt.colorbar(epic, ticks=np.arange(0, no2, 2))  # Ticks every 2 units
        
    # Set the tick labels manually, add '+' to the last label
    tick_labels = [str(tick) for tick in range(0, no2-2, 2)] + ['28+']
    cbar.set_ticklabels(tick_labels, size=13)
    ax.tick_params(axis='both', labelsize=13)
    ax.grid(alpha=0.2)
        
    # Optionally, change the color of the last tick label
    cbar.ax.get_yticklabels()[-1].set_color(over_value_color)
    cbar.ax.get_yticklabels()[-1].set_weight('bold')  # Optionally set the last label to bold
        
    # Beam in the lower left corner
    bmaj = header_temp['BMAJ']
    bmin = header_temp['BMIN']
    bmaj_au = bmaj * DISTANCE_TO_SOURCE * 10**3 * 206264.81 * np.pi/180
    bmin_au = bmin * DISTANCE_TO_SOURCE * 10**3 * 206264.81 * np.pi/180
    beam = Ellipse(xy=(-X_CUT + X_CUT/7, -Y_CUT + Y_CUT/7), width=bmin_au, 
                   height=bmaj_au, angle=header_temp['BPA'], color='grey',
                   alpha=0.5)
    plt.gca().add_patch(beam)
    # Optionally remove this:
    plt.gca().text(-X_CUT + X_CUT/7, -Y_CUT + Y_CUT/7,'beam', color='white', ha='center', va='center', alpha = 0.6)
    
    # Optional - for the valid region on the temperature map
    radius = 800
    center = (0, 0)
    theta = np.linspace(0, 2*np.pi, 100)
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)
    plt.plot(x, y, color='red', ls='--', alpha=0.7,label='Valid Region: r $\\lesssim$ '+str(radius)+'AU')
        
    
    plt.xlabel("AU", size=14)
    plt.ylabel("AU", size=14)
    plt.title('Toomre Q - Keplerian - CH$_{3}^{\:13}$CN', size=15)
    plt.legend(fontsize=14, loc='upper right')

    plt.show()


def mass_calc(wave,T,F_v,x,y):
    """
    wavelength : in mm
    T : in K
    k_v : 1.99 cm^2 g^-1
    F_v : From contintuum in Jy
    d : in pc, (8.1 kpc)
    M : solar masses
    """
    d = DISTANCE_TO_SOURCE * 10 ** 3
    g=100 # dust-gas ratio
    k_v = 1.99 
    
    r = np.sqrt((x**2+y**2)) / au
    
    part1 = (np.exp(1.439 * (1.36*(T/10))**(-1) ) -1 )
    
    part2 = ((k_v/0.01)**(-1))
    part3 = (F_v)*((d/100)**2)
    part4 = ((wave)**3)
    
    if not r < 10:
        M = 0.12*part1*part2*part3*part4 * g
    else:
        M = np.nan
        
    
    return M


def galactic_coords(data, header):

    ctype1 = header['CTYPE1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    crpix1 = header['CRPIX1']

    ctype2 = header['CTYPE2']
    crval2 = header['CRVAL2']
    cdelt2 = header['CDELT2']
    crpix2 = header['CRPIX2']

    # Create the coordinate arrays
    x_gal = (crval1 + (np.arange(data.shape[3]) - crpix1) * cdelt1)
    y_gal = (crval2 + (np.arange(data.shape[2]) - crpix2) * cdelt2)

    x_m = x_gal * DISTANCE_TO_SOURCE * 10**3 * 206264.81 * np.pi/180 * au
    y_m = y_gal * DISTANCE_TO_SOURCE * 10**3 * 206264.81 * np.pi/180 * au

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
    except RuntimeWarning:
        ang_vel = np.nan

    return ang_vel

def speed_sound(T):

    Molar = 2.8 * 10**(-3)
    gam = 1.4
    cs = np.sqrt(gam*R*T/Molar)

    return cs

def match_shapes_with_nan(data):
    
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
    """
    Have left this out as you'd probably want to do this coord trans 
    more carefully using astropy packages and more carefully consider
    the coord units used etc. 
    
    """
    phi = np.deg2rad(170.63)
    inc = np.deg2rad(38)
    alpha = phi #+ np.pi
    
    #x_prime = x*np.cos(alpha) + y*np.sin(alpha)
    #y_prime = x*np.sin(alpha)/np.cos(inc) + y*np.cos(alpha)/np.cos(inc)
    
    x_prime,y_prime = x,y

    return x_prime,y_prime
    
def trimming(arr):
    temp=[]
    for line in arr:
        if 10 < line[0]/au < 2100:
            temp.append(line)
    
    return np.vstack((temp))


get_data()
