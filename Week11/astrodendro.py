# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 11:30:26 2023

@author: nickl
"""

from astropy.io import fits
from astrodendro import Dendrogram
from tkinter import Tk, filedialog
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import PowerNorm
from astrodendro.analysis import PPStatistic

######################################
minValue = 0.0008 # minimum temperature of source
minDelta = 0.00001 # the minimum difference in intensity between any two structures
minPix = 5 #minimum number of pixels for a source
SaveFig = False
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
    dendrogram = Dendrogram.compute(data, min_value=minValue, min_delta=minDelta, min_npix=minPix)

    p = dendrogram.plotter()

    fig, ax = plt.subplots(figsize=(10, 6))
    p.plot_tree(ax)
    ax.set_xlabel('Structure')
    ax.set_ylabel('Flux / Jy/beam')
    ax.set_title('Dendrogram of Continuum Image')

    if SaveFig:
        plt.savefig(fname=SaveFigName, bbox_inches='tight', dpi=400)

    plt.show()

    fig, ax2 = plt.subplots(figsize=(10, 10))
    data_2d = data[0, 0, :, :]

    im = ax2.imshow(data_2d, origin='lower', cmap=plt.cm.Blues, vmax= np.nanmax(data_2d)/4.3)

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

    #this stuff below doesnt work

    p2 = dendrogram.plotter()
    fig = plt.figure()
    ax3 = fig.add_subplot(1, 1, 1)
    ax3.imshow(data_2d, origin='lower', interpolation='nearest', cmap=plt.cm.Blues, vmax=np.nanmax(data_2d)/4.3)
    data_2d[np.isnan(data_2d)] = 0


    # Create a plotter for the 2D dendrogram
    p2= dendrogram.plotter()

    p2.plot_contour(ax3, color='black')
    p2.plot_contour(ax3, structure=8, lw=3, colors='red') #highlighting
    p2.plot_contour(ax3, structure=24, lw=3, colors='orange')
def execute():
    file = get_data(initial_dir)
    dendrogram(file)

execute()
