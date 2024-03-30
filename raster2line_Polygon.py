# Import necessary libraries
from osgeo import gdal, ogr, osr  
import geopandas as gpd  
import pandas as pd  
import matplotlib.pyplot as plt 
import os 
from unidecode import unidecode  
from shapelysmooth import taubin_smooth 
from shapely.geometry import Polygon, MultiPolygon  
from shapely import wkt  

#%% Function to convert raster to shapefile
def convert_raster_to_shapefile(input_raster_path, output_shapefile_path):
    try:
        # Open the GeoTIFF raster
        raster = gdal.Open(input_raster_path)

        # Get the raster band
        band = raster.GetRasterBand(1)

        # Get the projection information from the raster
        proj = raster.GetProjection()

        # Create a SpatialReference object and import the projection
        shp_proj = osr.SpatialReference()
        shp_proj.ImportFromWkt(proj)

        # Create the shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')
        
        # If the output shapefile already exists, delete it first (optional)
        if os.path.exists(output_shapefile_path):
            driver.DeleteDataSource(output_shapefile_path)

        create_shp = driver.CreateDataSource(output_shapefile_path)

        # Create the shapefile layer with the same CRS as the raster
        shp_layer = create_shp.CreateLayer('layername', srs=shp_proj)

        # Create a new field 'ID' of type Integer
        new_field = ogr.FieldDefn('ID', ogr.OFTInteger)
        shp_layer.CreateField(new_field)

        # Polygonize the raster band and add it to the shapefile layer
        gdal.Polygonize(band, None, shp_layer, 0, [], callback=None)

        # Destroy the shapefile and raster data sources
        create_shp.Destroy()
        raster = None

        # Open the shapefile for writing (mode 1)
        shp_ds = ogr.Open(output_shapefile_path, 1)

        # Get the layer from the shapefile
        shp_layer = shp_ds.GetLayer()

        # Create an empty list to store features to delete
        features_to_delete = []

        # Loop through all features in the layer and find features with 'ID' field value of 0
        for feature in shp_layer:
            if feature.GetField('ID') == 0 or feature.GetField('ID') == -9999:
                features_to_delete.append(feature.GetFID())

        # Loop through the list of features to delete and delete them
        for fid in features_to_delete:
            shp_layer.DeleteFeature(fid)

        # Commit the changes and close the shapefile
        shp_ds = None
        return True  # Indicates success
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return False  # Indicates failure 

#%%
# Function to smooth and filter shapefile
def smooth_and_filter_shapefile(input_shapefile_path, output_shapefile_path, factor=0.2, mu=0.5, steps=10, min_perimeter=10, min_area=1000):
    # Read shapefile 
    gdf = gpd.read_file(input_shapefile_path)

    filtered_smoothed_polygons = []

    for geom in gdf.geometry:
        if geom.length >= min_perimeter and geom.area >= min_area:
            # Taubin_smooth smoothing 
            smoothed_geom = taubin_smooth(geom, factor, mu, steps)
        
            filtered_smoothed_polygons.append(smoothed_geom)

    multi_smoothed = MultiPolygon(filtered_smoothed_polygons)

    smoothed_gdf = gpd.GeoDataFrame({'geometry': [multi_smoothed]}, crs=gdf.crs)

    # Save shapefile
    smoothed_gdf.to_file(output_shapefile_path)

#%%
# Function to remove holes from shapefile
def remove_holes(input_shapefile_path, output_shapefile_path, eps=1000):
    # Read shapefile
    gdf = gpd.read_file(input_shapefile_path)

    new_geoms = []

    for geom in gdf.geometry:
        if isinstance(geom, MultiPolygon):
            new_geoms.extend([
                MultiPolygon([Polygon(p.exterior, holes=[
                    interior for interior in p.interiors if Polygon(interior).area > eps
                ]) for p in geom.geoms])
            ])
        elif isinstance(geom, Polygon):
            new_geoms.append(
                Polygon(geom.exterior, holes=[
                    interior for interior in geom.interiors if Polygon(interior).area > eps
                ])
            )

    gdf = gpd.GeoDataFrame(geometry=new_geoms, crs=gdf.crs)

    # Save shapefile
    gdf.to_file(output_shapefile_path)

#%%
# Function to create boundary shapefile
def create_boundary_shapefile(input_file, output_file):
    # Read the initial polygon file
    gdf = gpd.read_file(input_file)

    # Create a new GeoDataFrame to contain the boundaries of the polygons
    line_gdf = gpd.GeoDataFrame(columns=['geometry'], crs=gdf.crs)

    # Iterate through each polygon in the initial GeoDataFrame and convert them to boundaries
    for idx, row in gdf.iterrows():
        polygon = row['geometry']
        boundary = polygon.boundary
        line_gdf = pd.concat([line_gdf, gpd.GeoDataFrame({'geometry': [boundary]}, crs=gdf.crs)])

    # Save the new GeoDataFrame containing boundaries to a shapefile
    line_gdf.to_file(output_file)

    print("Created shapefile containing boundaries of polygons.")

#%%
# Function to calculate intersection lengths
def calculate_intersection_length(lines_path, polygons_path):
    # Read shapefile
    lines = gpd.read_file(lines_path)
    polygons = gpd.read_file(polygons_path)
    
    result = {}

    # Iterate through lines
    for i, line in lines.iterrows():
        # Iterate through polygons
        for j, polygon in polygons.iterrows():
            if line['geometry'].intersects(polygon['geometry']):
                # Calculate overlapping length
                intersection = line['geometry'].intersection(polygon['geometry'])
                length = intersection.length
                # Save the result
                result[line['id']] = length

    # Convert dictionary to array
    result_array = list(result.items())
    
    return result_array

#%% Define input parameters
d_m = ''
year = '2020'
method = 'Ảnh tỷ số'
province = "Thanh_Hoa_2A"

input_raster = f'Output/Img_water/{province}/{year}{d_m}/Water {year}.TIF'

output_path = f'Output/Shapfiles/{province}/{year}{d_m}'

if not os.path.exists(output_path):
    os.makedirs(output_path)
    
output_shapefile = f'{output_path}/Polygonize_{year}_{method}.shp'

#%%
# Call the function to convert the raster to a shapefile
convert_raster_to_shapefile(input_raster, output_shapefile)

#%% Smoothing and filtering
input_shapefile_path = output_shapefile
output_shapefile_path = f'{output_path}/Smooth_Polygon_{year}_{method}.shp'
factor = 0.01  # Smoothing factor
mu = 0.01  # Mu parameter
steps = 0  # Iterations
min_perimeter = 100  # Minimum perimeter to retain polygon 
min_area = 1000  # Minimum area to retain polygon

smooth_and_filter_shapefile(input_shapefile_path, output_shapefile_path, factor, mu, steps, min_perimeter, min_area)

#%% Remove Holes
shapefile_path = f'{output_path}/Smooth_Polygon_{year}_{method}.shp'
output_shapefile_path = f'{output_path}/Smooth_Holes_{year}_{method}.shp'

remove_holes(shapefile_path, output_shapefile_path, eps=500)

#%% Final Polygon Image
Img_report = f'Output/Img_report/{year}/{province}'
if not os.path.exists(Img_report):
    os.makedirs(Img_report)
    
gdf = gpd.read_file(output_shapefile_path)

# Create a figure and axis with specific size
fig, ax = plt.subplots(figsize=(5, 5))  

# Plot the shapefile
gdf.plot(ax=ax)

plt.title(f' Ma River ({year})')
plt.savefig(f'{Img_report}/Polygon Sông Mã xã Vĩnh Hòa {province} ({year}_{method}).png')
plt.axis('off')
plt.show()

#%%
# Create boundary shapefile
input_file = output_shapefile_path
output = f'Output/Line/{year}/{province}/{year}/Line_{year}_{method}.shp'
if not os.path.exists(f'Output/Line/{year}/{province}/{year}'):
    os.makedirs(f'Output/Line/{year}/{province}/{year}')
create_boundary_shapefile(input_file, output)

