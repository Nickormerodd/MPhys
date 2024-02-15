"""
Created on Tue Nov  7 14:14:52 2023

@author: Christopher

path:
    /mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/w1-3/casa_immoments.py
"""

# CASA script to automate the immoments process and export to FITS

# /mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/Week_4-6/Immoments_automation/Automating_casa_immoments_k=7-8_temp.py

# Importing necessary CASA packages
from casatasks import immoments, exportfits
from astropy.io import fits
from spectral_cube import DaskSpectralCube
from spectral_cube import SpectralCube as sc
import numpy as np

k_number = [0,1,2,3,4,5,6,7,8]
#k_number = [7]

for k in k_number:
    

    # Define the file paths as variables
    cube_file_path = f'/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/new_data/CH3CN_K/CH3CN_small_k={k}_cal.fits'
    # Update with your actual file path
    region_file_path = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/new_data/smallregion.crtf'  # Update with your actual file path
    output_file_path = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/new_data/CH3CN_moments/'  # Update with your actual file path
    
    # Open the FITS file and read the data cube
    data = fits.open(cube_file_path)
    cube = sc.read(data)
    
    # Compute statistics and RMS
    """
    stats = cube.statistics()
    rms = stats['rms']
    print(f'RMS: {rms}')
    
    cube_data = cube.filled_data[:]
    
    # Calculate the RMS using NumPy
    rms = np.std(cube_data)
    print(f'RMS: {rms}')
    """
    # Define the includepix parameter based on the RMS value
    includepix = [30, 100]
    
    # Moments to calculate: 0 (integrated intensity), 1 (velocity field), and 2 (velocity dispersion)
    moments = [0, 1, 2, 8]
    #moments = [1]
    
    for moment in moments:
        # Define the output image name for this moment
        moment_image_path = f"{output_file_path}k={k}_mom{moment}_30K_small_cal.image"
        moment_fits_path = f"{output_file_path}k={k}_mom{moment}_30K_small_cal.fits"
        
        # Run the immoments CASA task
        immoments(imagename=cube_file_path,
                  moments=[moment],  # Calculate this moment
                  includepix=includepix,
                  outfile=moment_image_path,
                  region=region_file_path)
    
        print(f"Moment {moment} image has been saved to {moment_image_path}")
    
        # Export the moment image to a FITS file
        exportfits(imagename=moment_image_path,
                   fitsimage=moment_fits_path,
                   overwrite=True,
                   velocity = False)
        
        print(f"Moment {moment} FITS file has been saved to {moment_fits_path}")

