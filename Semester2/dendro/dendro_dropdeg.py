# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 11:30:26 2023

@author: nickl

This code uses a fits file and produces a continuum image with overlaid 'leaves'
using astrodendro library. It plots highlighted contours on the continuum image
too, with an assosciated dendrogram
"""

from astropy.io import fits
from astrodendro import Dendrogram
from tkinter import Tk, filedialog
import numpy as np
import matplotlib.pyplot as plt

######################################
minValue = 0.0001 # minimum temperature of source rms is 1.787510269496e-5
minDelta = 0.00002 # the minimum difference in intensity between any two structures
# minpix is defined later minPix = 13 #minimum number of pixels for a source
SaveFig = False
initial_dir = 'C:/Users/nickl/linux/s2/new_data' # full continuum image
SaveFigName = initial_dir + '/../W1-4/dendro/dendrogram_contours.png'
#SaveFigName2 = initial_dir + '/../W1-4/dendro/dendrogram_leaves.png'
SaveFigName3 = initial_dir + '/../W1-4/dendro/dendrogram_tree.png'
#####################################

def get_data(initial_dir):
    # file picking with a temp window
    root = Tk()
    root.attributes('-topmost', True)  # Set the window to be always on top

    # Prompt the user to select files, starting in the specified directory
    file = filedialog.askopenfilename(initialdir=initial_dir, title="Select Folder, then select continuum")

    root.destroy() # removes temp window

    return file

def calculate_npix(file):
    with fits.open(file) as hdul:  # Open the FITS file and get the header
        header = hdul[0].header

    print(hdul[0].shape)
    theta_maj = header['BMAJ']  # Beam major axis in degrees
    theta_min = header['BMIN']  # Beam minor axis in degrees
    #print("Maj: ", theta_maj, "Min: ", theta_min)

    cdelt1 = header['CDELT1'] # Pixel scale in radians (from deg) along axis 1
    cdelt2 = header['CDELT2'] # Pixel scale in radians along axis 2
    A_pix = np.abs(cdelt1 * cdelt2)  # Area per pixel in square
    #print("Apix: ", A_pix)

    # Now we can calculate N_pix(min)
    npix_min = (2 * np.pi * theta_maj * theta_min) / (8 * np.log(2) * A_pix)
    roundPix = round(npix_min)+1
    print("nPix: ", roundPix)

    return roundPix


def dendrogram(file, minPix):
    data = fits.getdata(file)
    data_2d = data[0, 0, :, :]
    """
    dendrogram = Dendrogram.compute(data, min_value=minValue, min_delta=minDelta, min_npix=minPix)

    fig, ax2 = plt.subplots(figsize=(10, 10))
    im = ax2.imshow(data_2d, origin='lower', cmap=plt.cm.Blues, vmax= np.nanmax(data_2d)/10)

    for leaf in dendrogram.leaves:
        position, peak = leaf.get_peak() # this is in pixel coords so will need converting eventually
        #position is in form (0,0,y,x), peak is a number
        y, x = position[2], position[3]

        ax2.plot(x, y, 'x', color = 'red', markersize=3, markeredgewidth=1) #circle at each leaf


    if SaveFig:
        plt.savefig(fname=SaveFigName2, bbox_inches='tight', dpi=400)

    plt.show()
    """
    ## Next part is for continuum image with contours
    fig = plt.figure()
    ax3 = fig.add_subplot(1, 1, 1)
    im = ax3.imshow(data_2d, origin='lower', interpolation='nearest', cmap=plt.cm.Blues, vmax=np.nanmax(data_2d)/4.3)
    data_2d[np.isnan(data_2d)] = 0

    dendrogram_2d = Dendrogram.compute(data_2d, min_value=minValue, min_delta=minDelta, min_npix=minPix)
    p2= dendrogram_2d.plotter()

    for leaf in dendrogram_2d.leaves:
        position, peak = leaf.get_peak() # this is in pixel coords so will need converting eventually
        #position is in form (0,0,y,x), peak is a number but data_2d is (y,x)
        y, x = position[0], position[1]
        ax3.plot(x, y, 'x', color = 'red', markersize=3, markeredgewidth=1) #circle at each leaf
    # Create a plotter for the 2D dendrogram

    #p2.plot_contour(ax3, color='black')
    p2.plot_contour(ax3, structure=11, lw=3, color='red') #highlighting
    p2.plot_contour(ax3, structure=8, lw=3, color='orange')
    p2.plot_contour(ax3, structure=14, lw=3, color='green')
    p2.plot_contour(ax3, structure=4, lw=3, color='blue')

    plt.colorbar(im, ax=ax3, label='Flux / Jy/beam')
    ax3.set_xlabel('X Pixel')
    ax3.set_ylabel('Y Pixel')
    ax3.set_title('Continuum Image with Dendrogram Leaves')

    if SaveFig:
        plt.savefig(fname=SaveFigName, bbox_inches='tight', dpi=400)

    plt.show()

    fig = plt.figure() # For plotting the tree
    ax4 = fig.add_subplot(1, 1, 1)
    p2.plot_tree(ax4, color='black')
    p2.plot_tree(ax4, structure=[4], color='blue', lw=2)
    p2.plot_tree(ax4, structure=[11], color='red', lw=2) #highlighting
    p2.plot_tree(ax4, structure=[8], color='orange', lw=2)
    p2.plot_tree(ax4, structure=[14], color='green', lw=2)

    ax4.set_xlabel("Structure")
    ax4.set_ylabel("Flux")

    if SaveFig:
        plt.savefig(fname=SaveFigName3, bbox_inches='tight', dpi=400)

    plt.show()

def execute():
    file = get_data(initial_dir)
    npix = calculate_npix(file)
    dendrogram(file, npix)

execute()
