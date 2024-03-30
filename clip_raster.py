# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 09:27:57 2023

@author: kiet
"""

import os
from osgeo import gdal
import geopandas as gpd
import rasterio
from rasterio.mask import mask
#%%
def check_crs_shp2img(shapefile_path,tif_path):
    # Load a shapefile
    gdf = gpd.read_file(shapefile_path)
    print(gdf.crs)
    with rasterio.open(tif_path) as src:
        crs_target=src.crs
        print(crs_target)
    if(gdf.crs != crs_target):
        gdf_trans=gdf.to_crs(crs_target)
        output_path=shapefile_path[:-4]+'_transform'+shapefile_path[-4:]
        gdf_trans.to_file(output_path)
        return output_path
    return shapefile_path
#%%
def clip_img(input_path, output_path, shp_path):
    # Load the shapefile
    gdf = gpd.read_file(shp_path)

    # Iterate over each JP2 file in the input path
    for band_file in os.listdir(input_path):
        if band_file.endswith('.jp2'):
            input_file = os.path.join(input_path, band_file)
            output_file = os.path.join(output_path, band_file)

            with rasterio.open(input_file) as src:
                # Clip the raster using the shapefile
                out_image, out_transform = mask(src, gdf.geometry, crop=True)

                # Update the metadata
                out_meta = src.meta.copy()
                out_meta.update({
                    "driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform
                })

                # Write the clipped raster to the output path
                with rasterio.open(output_file, "w", **out_meta) as dest:
                    dest.write(out_image)

#%%
root_folder = "Data/Lam_Dong/2024"


output_root = "Output/Clip/Lam_Dong/2024"

shapefile_path='ShapeFiles/Lam_Dong1.shp'


years = [folder_name for folder_name in os.listdir(root_folder) if os.path.isdir(os.path.join(root_folder, folder_name))]


for year in years:
    input_folder = root_folder +'/'+ year 
    output_folder =output_root+'/'+year+'/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    
    clip_img(input_folder, output_folder,shapefile_path)