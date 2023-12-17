# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 13:24:44 2023

@author: nickl
"""

from astropy.io import fits
from astrodendro import Dendrogram
from tkinter import Tk, filedialog
import numpy as np
import matplotlib.pyplot as plt
from astrodendro import Dendrogram, pp_catalog
from astrodendro.analysis import PPStatistic
from astropy import units as u

######################################
minValue = 0.0008 # minimum temperature of source
minDelta = 0.00003 # the minimum difference in intensity between any two structures
minPix = 100 #minimum number of pixels for a source
SaveFig = True
initial_dir = 'C:/Users/nickl/linux/data/Fits_Jansky' # full continuum image
SaveFigName = initial_dir + '/../../Week11/dendrogram.png'
SaveFigName2 = initial_dir + '/../../Week11/continuum_sources.png'
#####################################

def get_data(initial_dir):
    # file picking with a temp window
    root = Tk()
    root.attributes('-topmost', True)  # Set the window to be always on top

    # Prompt the user to select files, starting in the specified directory
    file = filedialog.askopenfilename(initialdir=initial_dir, title="Select Folder, then select rings_final2.txt")

    root.destroy() # removes temp window

    return file

def dendrogram(file):
    data = fits.getdata(file)
    data_2d = data[0, 0, :, :]
    dendrogram = Dendrogram.compute(data_2d, min_value=minValue, min_delta=minDelta, min_npix=minPix)
    v = dendrogram.viewer()
    v.show()

def execute():
    file = get_data(initial_dir)
    dendrogram(file)

execute()
