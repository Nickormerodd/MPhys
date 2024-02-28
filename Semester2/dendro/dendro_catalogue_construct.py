# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 14:24:06 2023

@author: nickl

This code produces an astrodendro index catalogue for leaves in a dendrogram
(or can be adapted to do all structures by using 'd' instead of 'leaves').

It then manipulates the catalogue to clean it and also add 2 more columns,
exact area in pixels (rather than arsec^2) and a mass (calculated by assuming
a temperature and using a flux to mass conversion formula) in solar masses.

It also outputs a log mass distribution graph for the leaves
"""
from astropy.io import fits
from astrodendro import Dendrogram
from tkinter import Tk, filedialog
#from astropy.io import fits
from astrodendro import pp_catalog
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt


###################### values #########
pc = 3.086e+16
G = 6.6743e-11
SunMass = 1.98847e30
au = 1.496e+11
distance = 8.1*10**3  # Example distance in parsecs
kappa_nu = 1.99  # Dust opacity in cm^2/g
lambda_mm = 1.3  # Wavelength in millimeters
temperature_dust = 22  # Temperature dust estimate
temperature_x = 234  # Temperature in Kelvin
minValue = 0.0001 # minimum temperature of source rms is 1.787510269496e-5
minDelta = 0.00002 # the minimum difference in intensity between any two structures
minPix = 13 #minimum number of pixels for a source
SaveFig = False
initial_dir = 'C:/Users/nickl/linux/s2/new_data' # full continuum image
output_folder = 'C:/Users/nickl/linux/python/dendro/'
SaveFigName = initial_dir + '/../W1-4/dendro/logMassDistributions.png'
#####################################

def get_data(initial_dir):
    # file picking with a temp window
    root = Tk()
    root.attributes('-topmost', True)  # Set the window to be always on top

    # Prompt the user to select files, starting in the specified directory
    file = filedialog.askopenfilename(initialdir=initial_dir, title="Select Folder, then select continuum")

    root.destroy() # removes temp window

    return file

def catalogue(file):
    input_file = file
    output_filename = output_folder + "catalog_output.txt"

    image = fits.getdata(input_file)
    print(image.shape)
    image = image[0, 0, :, :]
    print(image.shape)

    d = Dendrogram.compute(image, minValue, minDelta, minPix)
    leaves = [structure for structure in d.all_structures if structure.is_leaf]
    print(len(leaves))
    hdu_list = fits.open(input_file)
    header = hdu_list[0].header
    hdu_list.close()

    metadata = {}
    metadata['data_unit'] = u.Jy / u.beam
    metadata['spatial_scale'] =  header['CDELT1']*3600 * u.arcsec
    metadata['beam_major'] =  header['BMAJ'] *3600 * u.arcsec
    metadata['beam_minor'] =  header['BMIN'] *3600 * u.arcsec

    cat = pp_catalog(leaves, metadata)

    with open(output_filename, "w") as file:
        for line in cat:
            file.write(str(line) + '\n')

    print(f"Catalog saved to {output_filename}")

    input_file2 = output_folder + 'catalog_output.txt'
    int_file = output_folder + 'int.txt'
    output_file2 = output_folder + 'cleaned_cat.txt'

    with open(input_file2, 'r') as file:
        lines = file.readlines()

    # Separate the header and units
    header = lines[0]
    units = lines[1]

    # Clean and keep only lines with numerical data
    cleaned_lines = [header, units]
    for line in lines[2:]:
        if any(char.isdigit() for char in line):
            cleaned_lines.append(line)

    with open(int_file, 'w') as file:
        file.writelines(cleaned_lines)

    print(f'Cleaned file written to {int_file}')

    with open(int_file, 'r') as file:
        lines = file.readlines()

    # Extracting the header and the units
    header = lines[0]
    units = lines[1]

    # Taking every other line starting from the third line (lines 3, 5, 7, 9, etc.)
    selected_lines = lines[2::2]

    # Combining the header, units, and the selected lines
    final_lines = [header, units] + selected_lines

    # Writing to the output file
    with open(output_file2, 'w') as file:
        file.writelines(final_lines)

    print(f'Catalogue columns written to {output_file2}')

    return output_file2

def mass_calc(Jy, temperature):
    mass_sol = 2* 100 * 0.12 * (np.exp(1.439 * (lambda_mm)**(-1)*(temperature/10)**(-1)) - 1) * (Jy/(1*u.Jy)) *(kappa_nu/ 0.01)**(-1) * (distance/ 100)**2 * (lambda_mm)**3
    return mass_sol.value

def MassDist(file_path, file):

    hdu_list = fits.open(file)
    fits_header = hdu_list[0].header
    # Lists to store the 'Column4' values for mass calculation
    column4_values = []
    # Area in arcsec values for area in pixels calculation
    area_arcsec_values = []

    with open(file_path, 'r') as file:
        # Read the entire file into memory
        lines = file.readlines()

    # Modify the first two lines to include 'Mass' and 'Exact Area' column headers
    lines[0] = lines[0].strip() + '      Mass     Exact_Area\n  idx       '
    lines[1] = lines[1].strip() + '       Sol_Mass    pix\n'

    # Extract pixel scale from the FITS header
    cdelt1 = abs(fits_header['CDELT1']) * 3600  # arcsec per pixel
    pixel_area = cdelt1**2  # Area of each pixel in arcsec^2

    for line in lines[2:]:
        columns = line.strip().split()
        column4 = float(columns[3])  # Assuming the fourth column contains the flux in Jy
        area_arcsec = float(columns[2])  # Assuming the third column contains the area in arcseconds
        column4_values.append(column4)
        area_arcsec_values.append(area_arcsec)

    # Calculate mass and area in pixels for each entry
    calculated_masses = [mass_calc(Jy, temperature_x) for Jy in column4_values]
    areas_pixels = [area / pixel_area for area in area_arcsec_values]

    # Combine original data with new mass and area in pixels values
    output_area_and_mass = output_folder + 'mass_catalogue.txt'
    with open(output_area_and_mass, 'w') as outfile:
        outfile.writelines(lines[:2])  # Write the modified headers

        for i, line in enumerate(lines[2:]):
            # Append the mass and area in pixels to each line
            new_line = line.strip() + f" {calculated_masses[i]:.6f} {areas_pixels[i]:.1f}\n"
            outfile.write(new_line)

    print(f'Updated catalogue with mass and area in pixels written to {output_area_and_mass}')

# Note: The function `mass_calc` and the fetching of `fits_header` from the FITS file will remain unchanged.
# Ensure `fits_header` is obtained from the FITS file before calling `MassDist`.

    """ #Print statements
    max_flux_index = np.argmax(column4_values)
    max_flux = column4_values[max_flux_index]
    corresponding_mass = calculated_masses[max_flux_index]

    min_flux_index = np.argmin(column4_values)
    min_flux = column4_values[min_flux_index]
    corresponding_massmin = calculated_masses[min_flux_index]

    print(f'Max Flux Value: {max_flux}')
    print(f'Corresponding Mass: {corresponding_mass:.6f}')

    print(f'Min Flux Value: {min_flux}')
    print(f'Corresponding Mass min: {corresponding_massmin:.6f}')
    """
    return calculated_masses

def plot_bins(masses):
    log_bins = np.logspace(np.log10(min(masses)), np.log10(max(masses)), num=15)
    plt.hist(masses, bins=log_bins, edgecolor='k')

    plt.xscale('log')
    plt.xlabel('Mass (Solar Masses)')
    plt.ylabel('Frequency')
    plt.title('Mass Histogram')
    plt.grid(True, which="both", ls="-")
    plt.show()

    if SaveFig:
        plt.savefig(fname=SaveFigName, bbox_inches='tight', dpi=400)

def execute():
    file = get_data(initial_dir)
    output = catalogue(file)
    masses = MassDist(output, file)
    plot_bins(masses)

execute()
