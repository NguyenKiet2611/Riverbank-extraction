# -*- coding: utf-8 -*-
"""
Created on Sat Aug 5 14:23:01 2023
@author: kiet
"""

#% Import Lib
import os
import rasterio
import geopandas as gpd
import matplotlib.pyplot as plt
from rasterio.mask import mask
from shapely.geometry import mapping, box


def get_image_info(image_path):
    try:
        with rasterio.open(image_path) as src:
            width = src.width
            height = src.height
            num_bands = src.count
            crs = src.crs
            pixel_width, pixel_height = src.res
            transform=src.transform

        image_info = {
            'width': width,
            'height': height,
            'num_bands': num_bands,
            'crs': crs,
            'pixel_width': pixel_width,
            'pixel_height': pixel_height,
            'transform':transform
        }
        
        return image_info
    except Exception as e:
        print(f"Lỗi: {e}")
        return None


def read_shapefile(shapefile_path):
    try:
        gdf = gpd.read_file(shapefile_path)
        return gdf
    except Exception as e:
        print(f"Lỗi: {e}")
        return None


def plot_band_histogram(image_path, band_index=1):
    try:
        with rasterio.open(image_path) as src:
            band_data = src.read(band_index)
            plt.hist(band_data.flatten(), bins=50, color='blue', alpha=0.7)
            plt.xlabel('Giá trị Pixel')
            plt.ylabel('Số lượng Pixel')
            plt.title(f'Biểu đồ histogram - Băng tần {band_index}')
            plt.show()
    except Exception as e:
        print(f'Lỗi: {e}')

