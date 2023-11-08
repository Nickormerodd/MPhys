"""
Created on Wed Nov  8 12:59:29 2023

@author: Christopher

This code reads in a csv file which you have to specify below.
This csv file has the following:
1st column - file path to input file
2nd column - file path to region file
3rd coumn - channel range (format n-m)
4th column - moments (in format 0,1,2,...,n)

The first row isnt being read, so start from
2nd row onwards.

For any rows which you want the code to ignore,
put a # infront of the file path in column 1.
"""

import csv
from casatasks import immoments, exportfits, imsubimage
from astropy.io import fits
from spectral_cube import SpectralCube as sc
import ast
import os

# Function to write the subimage and compute moments
def file_writing(name, chan_1, chan_2, rms, file_path, region_path, moments):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    includepix = [rms, 100]

    subimage_path = f"{script_dir}/{name}_{chan_1}-{chan_2}.image"
    imsubimage(imagename=file_path,
               outfile=subimage_path,
               chans=f"{chan_1}~{chan_2}",
               overwrite=True)

    for moment in moments:
        moment_image_path = f"{script_dir}/{name}_moment{moment}.image"
        moment_fits_path = f"{script_dir}/{name}_moment{moment}.fits"

        # Run the immoments CASA task
        immoments(imagename=subimage_path,
                  moments=[moment],
                  includepix=includepix,
                  outfile=moment_image_path,
                  region=region_path,
                  overwrite=True)
        
        print(f"Moment {moment} image for {name} has been saved to {moment_image_path}")

        # Export the moment image to a FITS file
        exportfits(imagename=moment_image_path,
                   fitsimage=moment_fits_path,
                   overwrite=True,
                   velocity=False)

        print(f"Moment {moment} FITS file for {name} has been saved to {moment_fits_path}")

# Main function to orchestrate the process
def process_csv(csv_file):
    with open(csv_file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip the header row

        for row in reader:
            if row[0].startswith('#'):  # Skip rows where the file path starts with '#'
                continue

            file_path, region_path, channels_str, moments_str = row

            # Convert channels string to a tuple of integers
            chan_1, chan_2 = map(int, channels_str.split('-'))

            # Convert moments string to a list of integers
            moments = ast.literal_eval(f"[{moments_str}]")

            with fits.open(file_path) as data:
                cube = sc.read(data)
            
            stats = cube.statistics()
            rms = stats['rms'].value
            filename = os.path.splitext(os.path.basename(file_path))[0]
            name = filename.replace(".image.pbcor_line.galactic", "").replace("kelvin", "")
            
            file_writing(name, chan_1, chan_2, rms, file_path, region_path, moments)

if __name__ == "__main__":
    csv_file = '/path/to/your/instructions.csv'  # Replace with your actual CSV file path
    process_csv(csv_file)
