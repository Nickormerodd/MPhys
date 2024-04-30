# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 11:43:36 2024

@author: nickl

Code which tells you about the minimum spanning tree and nearest neighbour info
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from astropy.io import fits
import networkx as nx
from scipy.constants import Boltzmann, m_p

const_pc = 3.086e+16
const_au = 1.496e+11
SaveFig = True
SaveFigName = 'Histogram_seperation_50K.png'

data1 = "3 sigma"
data2 = "5 sigma"

def get_data():


    data_3sig = np.loadtxt('core_positions_lots.dat')
    data_5sig = np.loadtxt('core_positions9.dat')

    fits_image_path = 'C://Users//nickl//Linux//s2//new_data//continuum_gal_1.fits'
    hdulist = fits.open(fits_image_path)
    image_data = hdulist[0].data[0, 0, :, :] # Select the spatial slice

    return data_3sig, data_5sig, image_data

def filtering(data, name):
    """
    Filters the core positions to include only those within the specified x and y ranges.
    """
    x_range = (1900, 2070)
    y_range = (1900, 2070)
    # Apply the filtering based on the provided ranges for x and y
    filtered_data = data[(data[:, 0] >= x_range[0]) & (data[:, 0] <= x_range[1]) &
                         (data[:, 1] >= y_range[0]) & (data[:, 1] <= y_range[1])]
    print(f"Number of cores in {name} : {len(filtered_data):.0f}")
    return filtered_data

def AU(pixels):
    """
    Convert pixel size to astronomical units (AU) given the pixel scale and distance to the object.
    """
    cdelt1_deg = 2.138888886286E-06  # degrees/pixel, CDELT1 from the FITS header
    distance_kpc = 8.1  # Distance to the object in kiloparsecs

    # Convert the pixel scale from degrees to radians
    cdelt_rad = cdelt1_deg * (np.pi / 180.0)

    distance_au = distance_kpc * 1e3 * const_pc/const_au #kpc to au (1 parsec = 206265 AU)
    size_au = pixels * cdelt_rad * distance_au # Calculate the size in AU

    return size_au

def Jeans_length(vel, n_dens):

    length = 0.4 * const_pc * (vel/0.2) * (n_dens/10**2)**(-1/2)
    length_au = length / const_au
    return length_au

def Jeans_mass(vel, n_dens):

    mass = 2 * (vel/0.2)**3 * (n_dens/10**3)**(-1/2)

    return mass


def plot_tree(data, image_data):
    # neigh = NearestNeighbors(n_neighbors=2, algorithm='auto')
    neigh = NearestNeighbors(n_neighbors=len(data))

    neigh.fit(data) # Fit the model with the data
    distances, indices = neigh.kneighbors(data) # Find the distances and indices of the nearest neighbors for each point

    # Since we're interested in the distance to the nearest neighbor (excluding itself), we take the second column of the distances array
    nearest_neighbor_distances = distances[:, 1]

    ################## plotting on a fits file ############
    A = neigh.kneighbors_graph(data, mode='distance')
    csr_matrix_A = csr_matrix(A) # Convert to a CSR (Compressed Sparse Row) matrix, then calculate the MST
    mst = minimum_spanning_tree(csr_matrix_A)


    mst_dense = mst.toarray() # Convert the MST to a dense matrix and then to a NetworkX graph for easy plotting
    graph = nx.from_numpy_matrix(mst_dense)

    fig, ax = plt.subplots(1, 1, figsize=(10, 10)) # Plotting

    ax.imshow(image_data, cmap='gray', origin='lower', vmax=np.nanmax(image_data)/4.3) # Display the FITS image

    # Overlay the MST and core positions
    pos = {node: (data[node][0], data[node][1]) for node in graph.nodes()}
    nx.draw_networkx_edges(graph, pos, ax=ax, edge_color='red', width=2)
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_size=50, node_color='blue')

    plt.axis('off')
    plt.show()
    return nearest_neighbor_distances

def print_data(nearest_neighbor_distances):

    print("\nMinimum Spanning Tree Params:")
    print(f"\nMean distance to nearest neighbor (pix): {np.mean(nearest_neighbor_distances)}")
    print(f"Standard deviation of distances (pix): {np.std(nearest_neighbor_distances)}")
    print(f"Minimum distance to nearest neighbor (pix): {np.min(nearest_neighbor_distances)}")
    print(f"Maximum distance to nearest neighbor (pix): {np.max(nearest_neighbor_distances)}\n")

    print(f"Mean distance to nearest neighbor (AU): {AU(np.mean(nearest_neighbor_distances))}")
    print(f"Standard deviation of distances (AU): {AU(np.std(nearest_neighbor_distances))}")
    print(f"Minimum distance to nearest neighbor (AU): {AU(np.min(nearest_neighbor_distances))}")
    print(f"Maximum distance to nearest neighbor (AU): {AU(np.max(nearest_neighbor_distances))}\n")

    ## correction of 2.3 from nnn to nearest seperation, then extra 0.5 from 2D projection
    print(f"Mean seperation length (AU): {2.8*AU(np.mean(nearest_neighbor_distances))}\n")

    return

def Jeans_calculations():

    n_dens = 10**6 #cm^(-3) +/- 10**5

    T = 50 #K +/-20
    mu = 2.8 * 10**(-3)
    c_s = np.sqrt(1.4*8.314*T/mu) / 1000 #is 0.012 kms^-1
    print(f"c_s: {c_s:.2f} kms^-1")
    print("\nJeans Mass and Length Calculations:\n")
    print("Jeans Length:")

    vel_disp = 2.1 #kms^-1 +/- 0.2

    var_thermT = ((1/2*(np.sqrt(1.4*8.314/mu * T)) * 0.4/(0.2*1000) * const_pc/const_au * (n_dens)**(-1/2))**2 * (10)**2)
    var_thermn = ((0.4 * const_pc/const_au * (c_s/0.2) * (-1/2)*(1/10**2)**(-1/2)*(n_dens)**(-3/2))**2 * (10**5)**2)
    err_thermal_jeans = np.sqrt(var_thermT + var_thermn)
    print(f"err_thermT: {np.sqrt(var_thermT):.2f} + err_thermn: {np.sqrt(var_thermn):.2f}\nerr_thermal_jeans: {err_thermal_jeans:.2f}")

    var_turbv = ((0.4/(0.2) * const_pc/const_au * (n_dens/10**1)**(-1/2))**2 * (2.1)**2)
    var_turbn = ((0.4 * const_pc/const_au * (vel_disp/0.2) * (-1/2)*(1/10**3)**(-1/2)*(n_dens)**(-3/2))**2 * (10**4)**2)
    err_turb_jeans = np.sqrt(var_turbv+var_turbn)
    print(f"err_turbv: {np.sqrt(var_turbv):.2f} + err_turbn: {np.sqrt(var_turbn):.2f}\nerr_turb_jeans: {err_turb_jeans:.2f}\n")

    arberror = 1
    nerror = 0.5

    thermal_jeans = Jeans_length(c_s, n_dens), err_thermal_jeans * arberror  ### value and error
    turbulent_jeans = Jeans_length(vel_disp, n_dens), err_turb_jeans * arberror### value and error

    print(f"Thermal Jeans length (AU): {thermal_jeans[0]:.0f} +/- {thermal_jeans[1]:.0f}")
    print(f"Turbulent Jeans length (AU): {turbulent_jeans[0]:.0f} +/- {turbulent_jeans[1]:.0f}\n")
    print("Jeans Mass:")
    var_thermT_mass = nerror**2 * ((2* (3/2) * (1.4*8.314*T/mu)**(1/2) * (1/(0.2*1000))**(3) * const_pc/const_au * (n_dens)**(-1/2))**2 * (20)**2)
    var_thermn_mass = nerror * ((2 * const_pc/const_au * (c_s/0.2)**(3) * (-1/2)*(1/10**1)**(-1/2)*(n_dens)**(-3/2))**2 * (10**2)**2)
    err_thermal_jmass = np.sqrt(var_thermT_mass + var_thermn_mass)
    print(f"err_thermT_mass: {np.sqrt(var_thermT_mass):.2f} + err_thermn_mass: {np.sqrt(var_thermn_mass):.2f}\nerr_thermal_jeans_mass: {err_thermal_jmass:.2f}")

    var_turbv_mass = nerror * ((2*3*vel_disp**2/((0.2*1000)**(3)) * const_pc/const_au * (n_dens/10**1)**(-1/2))**2 * (2.1*1000)**2)
    var_turbn_mass = nerror**2 * ((2 * const_pc/const_au * (vel_disp/(0.2))**3 * (-1/2)*(1/10**3)**(-1/2)*(n_dens)**(-3/2))**2 * (10**1)**2)
    err_turb_jmass = np.sqrt(var_turbv_mass + var_turbn_mass)
    print(f"err_turbv_mass: {np.sqrt(var_turbv_mass):.2f} + err_turbn_mass: {np.sqrt(var_turbn_mass):.2f}\nerr_turb_jeans_mass: {err_turb_jmass:.2f}\n")

    thermal_jmass = Jeans_mass(c_s, n_dens), err_thermal_jmass
    turbulent_jmass = Jeans_mass(vel_disp, n_dens), err_turb_jmass

    print(f"Thermal Jeans mass (Sol Mass): {thermal_jmass[0]:.2f} +/- {thermal_jmass[1]:.2f}")
    print(f"Turbulent Jeans mass (Sol Mass): {turbulent_jmass[0]:.0f} +/- {turbulent_jmass[1]:.0f}\n")

    return thermal_jeans, turbulent_jeans


def plot_histogram(nn_dist_3, nn_dist_5, thermal_jeans, turb_jeans):
    # Calculate the logarithm of separations
    nn_dist_3 = AU(nn_dist_3) * 2.8
    nn_dist_5 = AU(nn_dist_5) * 2.8
    #print(nearest_neighbor_distances)
    log_distances_3 = np.log10(nn_dist_3)
    log_distances_5 = np.log10(nn_dist_5)
    #print(log_distances)

    # Create a histogram of log separations
    bin_edges = np.linspace(3, 4, num=21)  # 11 bins between 3 and 4
    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    counts_3, _ = np.histogram(log_distances_3, bins=bin_edges)
    counts_5, _ = np.histogram(log_distances_5, bins=bin_edges)

    mean_separation =  AU(np.mean(nn_dist_3))
    mean_err =  AU(np.std(nn_dist_3))

    mean_logsep = np.sum(bin_centres*counts_3)/np.sum(counts_3)
    mean_err_logsep = 1/(np.sum(counts_3)) * np.sum(np.abs(counts_3*(bin_centres-mean_logsep))) + mean_err/mean_separation * np.sqrt(1/len(log_distances_3))
    mean_err_sep = mean_err_logsep * 10**(mean_logsep)

    # Plot the histogram

    plt.figure(figsize=(8, 6))

    #for weighted mean
    plt.axvspan((mean_logsep - mean_err_logsep/2), (mean_logsep + mean_err_logsep/2), color='black', alpha=0.5, linestyle='--',
                label=f'Mean Separation: {round(10**(mean_logsep),-1):.0f} +/- {round(mean_err_sep,-1):.0f} AU')
    """ for the actual mean rather than weighted
    plt.axvspan(np.log10(mean_separation - mean_err/2), np.log10(mean_separation + mean_err/2), color='black', alpha=0.5, linestyle='--',
                label=f'Mean Separation: {mean_separation:.0f} +/- {mean_err:.0f} AU')
    """
    plt.grid(True)
    plt.axvspan(np.log10(thermal_jeans[0] - thermal_jeans[1]/2), np.log10(thermal_jeans[0] + thermal_jeans[1]/2), color='green', alpha=0.5,
                label='Thermal Jeans length')
    plt.axvspan(np.log10(turb_jeans[0] - turb_jeans[1]/2), np.log10(turb_jeans[0] + turb_jeans[1]/2), color='red', alpha=0.5,
                label='Turbulent Jeans length')
    plt.bar(bin_centres, counts_3, width=0.05, align='center', color = 'khaki', label = '3$\sigma$')
    plt.bar(bin_centres, counts_5, width=0.05, align='center', hatch='//', color = 'khaki', label = '5$\sigma$')
    plt.xlabel('Log Separation [AU]')
    plt.ylabel('Frequency')
    plt.title('Frequency vs Log Separation')
    print("Final weighted mean calculation:\n")
    print(f"Weighted log(mean) seperation:  {mean_logsep:.2f} +/- {mean_err_logsep:.2f}",
          f"\nWeighted mean seperation (AU): {10**(mean_logsep):.0f} +/- {mean_err_sep:.0f}")

    plt.legend()
    if SaveFig:
        plt.savefig(fname=SaveFigName, bbox_inches='tight', dpi=400)

    plt.show()

def main():
    data_3sig, data_5sig, image_data = get_data()
    data_3sig = filtering(data_3sig, data1)
    data_5sig = filtering(data_5sig, data2)
    nn_dist_3 = plot_tree(data_3sig, image_data)
    nn_dist_5 = plot_tree(data_5sig, image_data)
    print_data(nn_dist_3)
    thermal_jeans, turb_jeans = Jeans_calculations()
    plot_histogram(nn_dist_3, nn_dist_5, thermal_jeans, turb_jeans)
    return

main()

