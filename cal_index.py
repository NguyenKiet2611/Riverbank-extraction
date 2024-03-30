# -*- coding: utf-8 -*-
"""
Created on Sat Jan 6 00:12:56 2024

@author: kiet2
"""

def calculate_NDVI(red_band, nir_band):
    """
    Calculate the NDVI index from red band and near-infrared (NIR) band data.
    
    :param red_band: Numpy array containing red band data
    :param nir_band: Numpy array containing near-infrared (NIR) band data
    :return: Numpy array containing NDVI values
    """
    # Convert pixel values to floating point format and avoid division by zero
    red_band = red_band.astype(float)
    nir_band = nir_band.astype(float)

    ndvi = (nir_band - red_band) / (nir_band + red_band)
    
    return ndvi

def calculate_SAVI(red_band, nir_band, L=0.5):
    """
    Calculate the SAVI index from red band and near-infrared (NIR) band data.
    
    :param red_band: Numpy array containing red band data
    :param nir_band: Numpy array containing near-infrared (NIR) band data
    :param L: Parameter L (default is 0.5)
    :return: Numpy array containing SAVI values
    """
    # Convert pixel values to floating point format and avoid division by zero
    red_band = red_band.astype(float)
    nir_band = nir_band.astype(float)

    savi = ((nir_band - red_band) / (nir_band + red_band + L)) * (1.0 + L)
    
    return savi

def calculate_NDMI(nir_band, swir_band):
    """
    Calculate the NDMI index from near-infrared (NIR) band and shortwave infrared (SWIR) band data.
    
    :param nir_band: Numpy array containing near-infrared (NIR) band data
    :param swir_band: Numpy array containing shortwave infrared (SWIR) band data
    :return: Numpy array containing NDMI values
    """
    # Convert pixel values to floating point format and avoid division by zero
    nir_band = nir_band.astype(float)
    swir_band = swir_band.astype(float)

    ndmi = (nir_band - swir_band) / (nir_band + swir_band)
    
    return ndmi

def calculate_MNDWI(green_band, swir_band):
    """
    Calculate the MNDWI index from green band and shortwave infrared (SWIR) band data.
    
    :param green_band: Numpy array containing green band data
    :param swir_band: Numpy array containing shortwave infrared (SWIR) band data
    :return: Numpy array containing MNDWI values
    """
    # Convert pixel values to floating point format and avoid division by zero
    green_band = green_band.astype(float)
    swir_band = swir_band.astype(float)
    
    mndwi = (green_band - swir_band) / (green_band + swir_band)
    
    return mndwi

def calculate_NDSI(green_band, swir_band):
    """
    Calculate the NDSI index from green band and shortwave infrared (SWIR) band data.
    
    :param green_band: Numpy array containing green band data
    :param swir_band: Numpy array containing shortwave infrared (SWIR) band data
    :return: Numpy array containing NDSI values
    """
    # Convert pixel values to floating point format and avoid division by zero
    green_band = green_band.astype(float)
    swir_band = swir_band.astype(float)

    ndsi = (green_band - swir_band) / (green_band + swir_band)
    
    return ndsi

def calculate_NDWI_Feeter(nir_band, green_band):
    """
    Calculate the NDWI Feeter index from near-infrared (NIR) band and green band data.
    
    :param nir_band: Numpy array containing near-infrared (NIR) band data
    :param green_band: Numpy array containing green band data
    :return: Numpy array containing NDWI Feeter values
    """
    # Convert pixel values to floating point format and avoid division by zero
    nir_band = nir_band.astype(float)
    green_band = green_band.astype(float)

    ndwi_feeter = (green_band - nir_band) / (green_band + nir_band)
    
    return ndwi_feeter

def calculate_NDWI_Gao(nir_band, swir_band):
    """
    Calculate the NDWI Gao index from near-infrared (NIR) band and shortwave infrared (SWIR) band data.
    
    :param nir_band: Numpy array containing near-infrared (NIR) band data
    :param swir_band: Numpy array containing shortwave infrared (SWIR) band data
    :return: Numpy array containing NDWI Gao values
    """
    # Convert pixel values to floating point format and avoid division by zero
    nir_band = nir_band.astype(float)
    swir_band = swir_band.astype(float)
    
    ndwi_gao = (nir_band - swir_band) / (nir_band + swir_band)
    
    return ndwi_gao

def calculate_NDBI(nir_band, swir_band):
    """
    Calculate the NDBI index from near-infrared (NIR) band and shortwave infrared (SWIR) band data.
    
    :param nir_band: Numpy array containing near-infrared (NIR) band data
    :param swir_band: Numpy array containing shortwave infrared (SWIR) band data
    :return: Numpy array containing NDBI values
    """
    # Convert pixel values to floating point format and avoid division by zero
    nir_band = nir_band.astype(float)
    swir_band = swir_band.astype(float)

    ndwi_gao = (swir_band - nir_band) / (nir_band + swir_band)
    
    return ndwi_gao

def calculate_custom_index(bands, formula, epsilon=1e-20):
    """
    Calculate a custom index based on the formula and a list of input bands.

    :param bands: List of numpy arrays containing data from input bands, along with the name of each band.
                  Example: [("band1", array1), ("band2", array2), ...]
    :param formula: Formula to calculate the index (e.g., "(band1 - band2) / (band1 + band2)")
    :param epsilon: A small value to avoid division by zero (default is 1e-20)
    :return: A numpy array containing the calculated index values based on the formula
    """
    # Convert pixel values to floating point format and avoid division by zero
    bands = [(name, band.astype(float)) for name, band in bands]
    
    # Create a dictionary of band variables from the corresponding names
    band_variables = {name: band for name, band in bands}
    
    try:
        # Add the epsilon variable to the computation scope
        band_variables['epsilon'] = epsilon
        
        # Use the formula and band variables to calculate the index
        custom_index = eval(formula, band_variables)
        
        return custom_index
    except ZeroDivisionError:
        print("Error: Division by zero in the formula.")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None
