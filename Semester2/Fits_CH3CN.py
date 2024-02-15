"""
Created on Tue Nov  7 15:45:01 2023

@author: Christopher
filepath = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/Week_4-6/Fits_files_for_CH3CN/fitsfile_for_CH3CN.py'
"""

# Importing necessary CASA packages
from casatasks import imsubimage, exportfits, imstat
from astropy.io import fits
from spectral_cube import SpectralCube as sc

# Define the path to the original FITS file
input_fits_path = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/new_data/CH3CN_gal_1_kelvin.fits'

# Channel ranges for each k value
channel_ranges = {
    'k=0': (1267, 1284),
    'k=1': (1252, 1267),
    'k=2': (1223, 1250),
    'k=3': (1180, 1207),
    'k=4': (1116, 1150),
    'k=5': (1044, 1075),
    'k=6': (948, 965),
    'k=7': (834, 856),
    'k=8': (701, 732),
   #  'k=20': (250,320)
    
}

# Path to save the new FITS files
output_directory = '/mnt/c/Users/Christopher/onedrive/documents/4th_year_physics/mphys_project/s2/W1-3/CH3CN_K/'
"""
# Open the FITS file and create a DaskSpectralCube to calculate RMS
data = fits.open(input_fits_path)
cube = sc.read(data)
stats = cube.moment0().std()
rms = stats.value
print(f'RMS: {rms}', rms**2)
"""
# Define the includepix parameter based on the RMS value
#includepix = [0.154, 100]

# Process each channel range and create a new FITS file
for k, (start_chan, end_chan) in channel_ranges.items():
    # The CASA image name for the subset
    subset_image_name = f"{output_directory}CH3CN_0.154_{k}.image"
    #subset_image_name = f"{output_directory}SPW25_noiserms.image"
    # Run the imsubimage task to create a subset image
    imsubimage(imagename=input_fits_path,
               outfile=subset_image_name,
               chans=f'{start_chan}~{end_chan}',
               dropdeg=False,  # Drop the degenerate axis if there is one
               overwrite=True)
    
    # Define the output FITS name
    output_fits_name = f"{output_directory}CH3CN_{k}.fits"
    #output_fits_name = f"{output_directory}SPW25_noiserms.fits"
    # Run exportfits to save the subset image as a FITS file
    exportfits(imagename=subset_image_name,
               fitsimage=output_fits_name,
               overwrite=True,
               minpix=0.154,
               dropdeg=False,
               velocity=False)
               #minpix=3,
               #maxpix=100)  # Drop the degenerate axis in the FITS file

    print(f'Created FITS file for CH3CN in {output_fits_name}')

print('All FITS files created successfully.')
