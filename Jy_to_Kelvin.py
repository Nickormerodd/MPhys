"""
@author: Christopher
This code coverts all the Jy/beam to Kelvin. This is what I believe
is right.

This code creates a new fits file with the name change of 
.replace(.'fits','_kelvin.fits')

Be careful, this code overwrites any files you have that have the same name
--> to stop this put data.writeto(new_filename, overwrite=False) !!!
"""


import numpy as np
from astropy.io import fits
from spectral_cube import SpectralCube as sc
from tkinter import Tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import os
import radio_beam as beam
import scipy.constants as const
from astropy import units as u

def get_data():

    Tk().withdraw()  
    file_paths = filedialog.askopenfilenames()
    for file_path in file_paths:
        data = fits.open(file_path)
        data_2 = sc.read(data)
        filename = (os.path.basename(file_path).replace(".fits", "")).replace("_"," ").replace("."," ")
        main(data, filename, file_path)
        data.close()    
    return

def main(data, filename, path):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    header = data[0].header
    ref_freq = header['CRVAL3'] #reference frequency
    incr_freq = header['CDELT3'] #frequency increment
    header['BUNIT'] = 'K'

    beam_major = np.deg2rad(header['BMAJ'])  # Major axis of the beam 
    beam_minor = np.deg2rad(header['BMIN'])  # Minor axis of the beam
    beam_size = (np.pi/(4 *np.log(2))) * beam_major * beam_minor
    
    num_channels = data[0].shape[0]  
    frequencies = ((ref_freq + incr_freq * (np.arange(num_channels))) / 10**9 ) *u.GHz
    beam_size = beam_size* u.sr
    
    for i in range(num_channels):
        equiv = u.brightness_temperature(frequencies[i])
        data[0].data[i] = (data[0].data[i]*u.Jy/beam_size).to(u.K,equivalencies=equiv)
    
    m = data[0].shape[1]
    p = data[0].shape[2]

    x = np.linspace(0, p - 1, p)
    y = np.linspace(0, m - 1, m)
    X, Y = np.meshgrid(x, y)
    
    contourf = ax.contourf(X, Y, data[0].data[0] , cmap='magma',
                           levels = 30) # 'viridis')
    plt.colorbar(contourf, label='Brightness Temperature, K')

    new_filename = os.path.basename(path).replace('.fits','_kelvin.fits') 
    print(new_filename)
    
    data.writeto(new_filename, overwrite=True) # Be careful, overwrites any 
                                               # files you have!!!!
         
    return
    
get_data()
