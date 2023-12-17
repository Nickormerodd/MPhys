# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 14:24:06 2023

@author: nickl
"""

from astropy.io import fits
from astrodendro import Dendrogram, pp_catalog
from astropy import units as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SaveFigName = 'MassDistributions.png'
SaveFig = True

pc = 3.086e+16
G = 6.6743e-11
SunMass = 1.98847e30
au = 1.496e+11
distance = 8.1*10**3  # Example distance in parsecs
kappa_nu = 1.99  # Dust opacity in cm^2/g
lambda_mm = 1.3  # Wavelength in millimeters
temperature_dust = 22  # Temperature dust estimate
temperature_x = 234  # Temperature in Kelvin

def catalogue():
    input_file = 'E_continuum.image.pbcor.galactic.fits'
    output_filename = "catalog_output.txt"

    image = fits.getdata(input_file)
    print(image.shape)
    image = image[0, 0, :, :]
    print(image.shape)

    d = Dendrogram.compute(image, min_value=0.001, min_delta=0.00003, min_npix=80)

    hdu_list = fits.open(input_file)
    header = hdu_list[0].header
    hdu_list.close()

    metadata = {}
    metadata['data_unit'] = u.Jy / u.beam
    metadata['spatial_scale'] =  header['CDELT1']*3600 * u.arcsec
    metadata['beam_major'] =  header['BMAJ'] *3600 * u.arcsec
    metadata['beam_minor'] =  header['BMIN'] *3600 * u.arcsec

    cat = pp_catalog(d, metadata)

    with open(output_filename, "w") as file:
        for line in cat:
            file.write(str(line) + '\n')

    print(f"Catalog saved to {output_filename}")

    input_file = 'catalog_output.txt'  # Replace with the path to your file
    int_file = 'int.txt'
    output_file = 'cleaned_cat.txt'

    with open(input_file, 'r') as file:
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
    with open(output_file, 'w') as file:
        file.writelines(final_lines)

    print(f'Catalogue columns written to {output_file}')

    return output_file

def mass_calc(Jy, temperature):
    mass_sol = 100 * 0.12 * (np.exp(1.439 * (lambda_mm)**(-1)*(temperature/10)**(-1)) - 1) * (Jy/(1*u.Jy)) *(kappa_nu/ 0.01)**(-1) * (distance/ 100)**2 * (lambda_mm)**3
    return mass_sol.value

def MassDist(file_path):
    # Lists to store the 'Column4' values
    column4_values = []

    with open(file_path, 'r') as file:
        # Skip the first two rows
        next(file)
        next(file)

        for line in file:
            # Split the line into columns and convert the fourth column to a float
            columns = line.strip().split()
            column4 = float(columns[3])  # Assuming the fourth column contains the flux in Jy
            column4_values.append(column4)

    # Calculate the masses using the mass_calc function
    calculated_masses = [mass_calc(Jy, temperature_dust) for Jy in column4_values]

    # Define bins for the histogram (e.g., specify the number of bins)
    num_bins = 15  # Adjust the number of bins as needed
    plt.hist(calculated_masses, bins=num_bins, edgecolor='k')
    plt.xlabel('Mass (Solar Masses)')
    plt.ylabel('Frequency')
    plt.title('Mass Histogram')
    plt.grid(True)

    if SaveFig:
        plt.savefig(fname=SaveFigName, bbox_inches='tight', dpi=400)

    # Create a new text file with an extra column for mass
    output_file_with_mass = 'mass_catalogue.txt'

    with open(file_path, 'r') as infile, open(output_file_with_mass, 'w') as outfile:
        # Copy the header lines
        header_lines = [next(infile) for _ in range(2)]
        outfile.writelines(header_lines)

        # Write lines with mass values
        for line, mass in zip(infile, calculated_masses):
            outfile.write(line.strip() + f' {mass:.6f}\n')  # Add mass as a new column

    max_flux_index = np.argmax(column4_values)
    max_flux = column4_values[max_flux_index]
    corresponding_mass = calculated_masses[max_flux_index]

    min_flux_index = np.argmin(column4_values)
    min_flux = column4_values[min_flux_index]
    corresponding_massmin = calculated_masses[min_flux_index]

    print(f'Updated catalogue with mass written to {output_file_with_mass}')
    print(f'Max Flux Value: {max_flux}')
    print(f'Corresponding Mass: {corresponding_mass:.6f}')

    print(f'Min Flux Value: {min_flux}')
    print(f'Corresponding Mass min: {corresponding_massmin:.6f}')

def execute():
    output = catalogue()
    MassDist(output)

execute()
