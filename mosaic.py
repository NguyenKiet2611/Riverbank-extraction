# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 14:35:08 2023

@author: kiet2
"""
import os
import rasterio
from rasterio.merge import merge
import shutil

def mosaic_images(text, input_folder, output_path):
    """
    Automatically merges image files or saves the input image if only one file is present.

    Parameters:
        text (str): The ending part of the image file name to be merged, for example: "B1.TIF".
        input_folder (str): Path to the folder containing the image files to be merged or to save the input image.
        output_path (str): Path to store the merged image or the input image.

    Returns:
        None
    """
    # Check if the input folder exists
    if not os.path.exists(input_folder):
        print("Input folder does not exist.")
        return

    # Get a list of all files in the input folder
    all_files = os.listdir(input_folder)

    # Create a list of image files with filenames ending with the specified text
    image_files = [os.path.join(input_folder, filename) for filename in all_files if filename.endswith(text)]

    # Check if there's only one image file for merging
    if len(image_files) == 1:
        # If there's only one file, copy it to the output path
        shutil.copy(image_files[0], output_path)
        print(f"Only one image file. Copied the input image to {output_path}")
        return

    # If there are multiple files, perform image merging
    src_files_to_mosaic = []

    # Open the source image files and add them to the list
    for fp in image_files:
        src = rasterio.open(fp)
        src_files_to_mosaic.append(src)

    try:
        # Call the merge function to merge the image files
        mosaic, out_trans = merge(src_files_to_mosaic)

        # Get metadata from one of the source files
        out_meta = src.meta.copy()

        # Update metadata
        out_meta.update({
            "driver": "GTiff",
            "height": mosaic.shape[1],
            "width": mosaic.shape[2],
            "transform": out_trans,
        })

        # Save the merged image to the output path
        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(mosaic)

        print(f"Images have been merged and saved at {output_path}")
    except Exception as e:
        print(f"Error occurred during image merging: {str(e)}")


# Define the province
province='Lam_Dong/2024'

# Set the root folder and output root based on the province
root_folder = f"Output/Clip/{province}"
output_root = f"Output/Mosaic/{province}"

# Create a list of years from the subfolders in the root folder
years = [folder_name for folder_name in os.listdir(root_folder) if os.path.isdir(os.path.join(root_folder, folder_name))]

# Iterate over each year and process
for year in years:
    input_folder = root_folder+'/'+ year + '/'
    output_folder = output_root+'/' + year + '/'

    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Iterate over files in the input folder
    for filename in os.listdir(input_folder):
        # Extract the last 10 characters from the filename
        text = filename[-10:]
        input_file_path = input_folder
        output_file_path = output_folder + text

        # Call the mosaic_images function to perform image merging or saving
        mosaic_images(text, input_folder, output_file_path)
