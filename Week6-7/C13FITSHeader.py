# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#import numpy as np
from astropy.io import fits
import os
#from spectral_cube import SpectralCube as sc
#from spectral_cube import DaskSpectralCube
from tkinter import Tk
from tkinter import filedialog
#import matplotlib.pyplot as plt
#import radio_beam as beam
from scipy.constants import c
#from astropy import units as u
#from scipy.optimize import curve_fit
#from matplotlib.ticker import FuncFormatter
#import uncertainties.unumpy as unp
#from uncertainties import ufloat

initial_dir = '/raid/scratch/normerod'
output_filename = "Calibrated_C13.fits"
transition_freq = 220.398684 # GHz
obs_freq = 220.3765      # GHz

def get_data():
    # file picking with a temp window
    root = Tk()
    root.attributes('-topmost', True)  # Set the window to be always on top
    file = filedialog.askopenfilename(initialdir=initial_dir, title="Select File")
    root.destroy() # removes temp window
    
    # Extract the directory of the input file
    file_directory = os.path.dirname(file)
    
    # Combine the directory with the new filename to get the full path of the output file
    outfile = os.path.join(file_directory, output_filename)
    
    # Now you can use 'output_file' as the path for your output file
    print("Input file = ", file)
    print("Output file = ", outfile)

    return file, outfile

def change_header_info(file, outfile):
    """
    Update the header in the output FITS file based on the input FITS file.
    """
    # Read data and header from the input file
    with fits.open(file) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        print("OLD VELREF = ", header['VELREF'])
        print("OLD RESTFRQ = ", header['RESTFRQ'])
        # Update header values
        header['RESTFRQ'] = transition_freq * 10**9  # Converted to Hz
        VELREF = round((c/1000) * -(obs_freq - transition_freq) / transition_freq, 2)
        print(type(header['RESTFRQ']))
        print(type(VELREF))
        header['VELREF'] = VELREF
        #print(header)
        print("NEW VELREF = ", header['VELREF'])
        print("NEW RESTFRQ = ", header['RESTFRQ'])
        # Write the original data with the updated header to the output file
        fits.writeto(outfile, data, header, overwrite=True)
        return outfile
    
def execute():
    file, outfile = get_data()
    outfile = change_header_info(file, outfile)
    with fits.open(outfile) as hdul:
        header = hdul[0].header
        print("NEW VELREF = ", header['VELREF'])
    
execute()
