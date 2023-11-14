# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 18:17:57 2023

@author: nickl
"""
from astropy.io import fits
import numpy as np
from scipy.constants import c

########## params ################
calibration_freq = 220.70908  # GHz
transition_freq = 220.538635 # GHz
obs_freq = 220.517542       # GHz
input_filename = 'C:/Users/nickl/linux/week7/trims/CH3CN_nopix_k=7.fits'
output_filename = 'C:/Users/nickl/linux/week7/trims/CH3CN_cal_k=7.fits'
####################################

VELREF = (c/1000) * -(obs_freq - transition_freq) / transition_freq
#scale_factor = (c/1000) * (transition_freq - calibration_freq) / calibration_freq
#print(scale_factor)

def change_header_info(input_fits_file, output_fits_file):
    """
    Update the header in the output FITS file based on the input FITS file.
    """
    # Read data and header from the input file
    with fits.open(input_fits_file) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        # Update header values
        header['RESTFRQ'] = transition_freq * 10**9  # Converted to Hz
        header['VELREF'] = VELREF

        # Write the original data with the updated header to the output file
        fits.writeto(output_fits_file, data, header, overwrite=True)

# Example usage
change_header_info(input_filename, output_filename)
