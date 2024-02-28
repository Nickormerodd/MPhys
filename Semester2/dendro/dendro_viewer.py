# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 13:24:44 2023

@author: nickl

This code makes a continuum image overlaid with leaves and also
produces an interactive viewer straight from astrodendro library
"""

from astropy.io import fits
from astrodendro import Dendrogram
from tkinter import Tk, filedialog
import numpy as np
import matplotlib.pyplot as plt
#from astrodendro import Dendrogram, pp_catalog
#from astrodendro.analysis import PPStatistic
#from astropy import units as u

######################################
minValue = 0.0001 # minimum temperature of source rms is 1.787510269496e-5
minDelta = 0.00002 # the minimum difference in intensity between any two structures
minPix = 13 #minimum number of pixels for a source
SaveFig = False
initial_dir = 'C:/Users/nickl/linux/s2/new_data' # full continuum image
SaveFigName = initial_dir + '/../W1-4/dendro/dendrogram.png'
SaveFigName2 = initial_dir + '/../W1-4/dendro/continuum_sources.png'
#####################################

def get_data(initial_dir):
    # file picking with a temp window
    root = Tk()
    root.attributes('-topmost', True)  # Set the window to be always on top

    # Prompt the user to select files, starting in the specified directory
    file = filedialog.askopenfilename(initialdir=initial_dir, title="Select Folder, then select continuum")

    root.destroy() # removes temp window

    return file

def dendrogram(file):

    data = fits.getdata(file)
    dendrogram = Dendrogram.compute(data, min_value=minValue, min_delta=minDelta, min_npix=minPix)

    p = dendrogram.plotter()

    fig, ax = plt.subplots(figsize=(10, 6))
    p.plot_tree(ax, color= 'black')
    #p.plot_tree(ax, structure=2, color='red', lw=2, alpha=0.5)
    ax.set_xlabel('Structure')
    ax.set_ylabel('Flux / Jy/beam')
    ax.set_title('Dendrogram of Continuum Image')

    if SaveFig:
        plt.savefig(fname=SaveFigName, bbox_inches='tight', dpi=400)

    plt.show()

    fig, ax2 = plt.subplots(figsize=(10, 10))
    data_2d = data[0, 0, :, :]
    im = ax2.imshow(data_2d, origin='lower', cmap=plt.cm.Blues, vmax= np.nanmax(data_2d)/10)

    for leaf in dendrogram.leaves:
        position, peak = leaf.get_peak() # this is in pixel coords so will need converting eventually
        #position is in form (0,0,y,x), peak is a number
        y, x = position[2], position[3]

        ax2.plot(x, y, 'x', color = 'red', markersize=3, markeredgewidth=1) #circle at each leaf

    plt.colorbar(im, ax=ax2, label='Flux / Jy/beam')
    ax2.set_xlabel('X Pixel')
    ax2.set_ylabel('Y Pixel')
    ax2.set_title('Continuum Image with Dendrogram Leaves')

    if SaveFig:
        plt.savefig(fname=SaveFigName2, bbox_inches='tight', dpi=400)

    plt.show()

    data = fits.getdata(file)
    data_2d = data[0, 0, :, :]
    dendrogram = Dendrogram.compute(data_2d, min_value=minValue, min_delta=minDelta, min_npix=minPix)
    v = dendrogram.viewer()
    v.show()

def execute():
    file = get_data(initial_dir)
    dendrogram(file)

execute()
