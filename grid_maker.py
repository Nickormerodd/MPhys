"""
Created on Sat Nov 25 20:22:08 2023

@author: Christopher
"""

import matplotlib.pyplot as plt
from PIL import Image


# Function to create a 3x3 grid of images without squeezing and save it
def create_image_grid(image_paths, output_path,dpi=30):
    # Open the first image to get its size
    with Image.open(image_paths[0]) as img:
        img_width, img_height = img.size
        # Calculate the size of the figure in inches for a 3x3 grid
        figsize = (3 * img_width / dpi, 3 * img_height / dpi)

    # Create the figure with the calculated size
    fig, axs = plt.subplots(3, 4, figsize=figsize, dpi=dpi, 
                            gridspec_kw={'wspace':0, 'hspace':0})

    for ax in axs.flat:
        ax.axis('off')  # Turn off the axis

    for i, ax in enumerate(axs.flat):
        with Image.open(image_paths[i]) as img:
            ax.imshow(img)
            ax.set_aspect('equal')  # Maintain the aspect ratio

    # Adjust the layout to remove any extra whitespace
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)

    # Save the figure
    #fig.savefig(output_path, bbox_inches='tight', pad_inches=0)
    plt.show()
    plt.close(fig)  # Close the figure to free memory

# Dummy list of image paths for demonstration purposes
#image_paths = ['/path/to/image1.png', '/path/to/image2.png', ...]  # Replace with actual paths
output_path = "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/Grid_with_iso.png"  # Replace with your desired output file path


# Dummy list of image paths for demonstration purposes
image_paths = ["C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3CN_k=0_gal.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3CN_k=1_gal.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3CN_k=2_gal.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3-13-CN_k=1.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3CN_k=3_gal.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3CN_k=4_gal.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3CN_k=5_gal.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3-13-CN_k=2.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3CN_k=6_gal.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3CN_k=7_gal.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3CN_k=8_gal.png",
               "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/Images/CH3CN gal/CH3-13-CN_k=3.png",]  # Replace with your actual paths

# Check if exactly 9 images were selected
if len(image_paths) == 12:
    output_path = "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/image_grid.png"  # Output path for the grid image
    create_image_grid(image_paths, output_path)
else:
    print("Please select exactly 9 images.")

    #output_path = "C:/Users/Christopher/OneDrive/Documents/4th_year_physics/MPhys_Project/Week_1-3/momment_mopping/image_grid.png"# Set the correct output path

