# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 21:37:37 2023

@author: kiet2
"""


"""
Created on Mon Aug 14 19:03:58 2023

@author: kiet
"""
from osgeo import gdal, ogr, osr
from rasterio.plot import show
import os
import numpy as np
import rasterio
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import matplotlib.colors as mcolors
import cv2
from PIL  import Image
from scipy.spatial.distance import cdist
#%%
import sys
sys.path.append('Code/')  
from cal_index import calculate_custom_index
from get_Info import get_image_info,plot_band_histogram

#%% Function
def read_raster_as_array(raster_path, band_number=1):
    try:
        raster = gdal.Open(raster_path)
        if raster is None:
            raise Exception("Could not open the raster file.")
        
        band = raster.GetRasterBand(band_number)
        if band is None:
            raise Exception("Could not retrieve the specified band.")
        
        array = band.ReadAsArray()
        #array=array.astype(np.float)
        return array
    except Exception as e:
        print("An error occurred:", str(e))
        return None


    except Exception as e:
        print("An error occurred:", str(e))
        return None 
def save_tiff(image_array,file_path,transform,crs='EPSG:32648' ,n_band=1):
    """
    Save a NumPy array as a TIF file using the rasterio library.

    Parameters:
        image_array (numpy.ndarray): NumPy array containing image data.
        file_path (str): Path to the destination TIF file.

    Returns:
        None
    """
    height, width = image_array.shape
    dtype = image_array.dtype

    with rasterio.open(
        file_path,
        'w',
        driver='GTiff',
        height=height,
        width=width,
        count=n_band,  # Number of bands
        dtype=dtype,
        crs=crs,  # Coordinate reference system, adjust as needed
        transform = transform
    ) as dst:
        dst.write(image_array,1)  # Write array data to the first band

#%% Information Image
d_m='/12_07'
year='2023'
province='Thanh_Hoa'
path=f'Output/Mosaic/{province}/{year}{d_m}'
#%% Path Image

path_B=path+'/02_10m.jp2'
path_G=path+'/03_10m.jp2'
path_R=path+'/04_10m.jp2'

path_VRE1=path+'/05_20m.jp2'
path_VRE2=path+'/06_20m.jp2'
path_VRE3=path+'/07_20m.jp2'
path_VRE4=path+'/8A_20m.jp2'

path_NIR=path+'/08_10m.jp2'

path_WV=path+'/09_60m.jp2'

path_SWIR1=path+'/11_20m.jp2'
path_SWIR2=path+'/12_20m.jp2'

path_TCI=path+'/CI_10m.jp2'

image_info=get_image_info(path_G)
if image_info:
    print(f'Kích thước ảnh: {image_info["width"]} x {image_info["height"]} pixels')
    print(f'Số lượng băng tần: {image_info["num_bands"]}')
    print(f'Hệ tọa độ của ảnh: {image_info["crs"]}')
    print(f'Kích thước pixel: {image_info["pixel_width"]} x {image_info["pixel_height"]} meters')
#%% Read Image
# Read raster data and handle data type conversion
B=read_raster_as_array(path_B)

G = read_raster_as_array(path_G)

R= read_raster_as_array(path_R)

NIR = read_raster_as_array(path_NIR)

VRE1=read_raster_as_array(path_VRE1)
VRE1=cv2.resize(VRE1, dsize=(G.shape[1],G.shape[0] ), interpolation=cv2.INTER_NEAREST)

VRE2=read_raster_as_array(path_VRE2)
VRE2=cv2.resize(VRE2, dsize=(G.shape[1],G.shape[0] ), interpolation=cv2.INTER_NEAREST)

VRE3=read_raster_as_array(path_VRE3)
VRE3=cv2.resize(VRE3, dsize=(G.shape[1],G.shape[0] ), interpolation=cv2.INTER_NEAREST)

VRE4=read_raster_as_array(path_VRE4)
VRE4=cv2.resize(VRE4, dsize=(G.shape[1],G.shape[0] ), interpolation=cv2.INTER_NEAREST)

SWIR1 = read_raster_as_array(path_SWIR1)
SWIR1 =cv2.resize(SWIR1 , dsize=(G.shape[1],G.shape[0] ), interpolation=cv2.INTER_NEAREST)

SWIR2 = read_raster_as_array(path_SWIR2)
SWIR2 =cv2.resize(SWIR2 , dsize=(G.shape[1],G.shape[0] ), interpolation=cv2.INTER_NEAREST)

WV=read_raster_as_array(path_WV)
WV =cv2.resize(WV , dsize=(G.shape[1],G.shape[0] ), interpolation=cv2.INTER_NEAREST)

TCI=np.stack((read_raster_as_array(path_TCI,band_number=1),
             read_raster_as_array(path_TCI,band_number=2),
             read_raster_as_array(path_TCI,band_number=3)),axis=-1
             )
#%% Image Report
Img_report=f'Output/Img_report/{province}/{year}{d_m}'
if(os.path.exists(Img_report)==False):
    os.makedirs(Img_report)
if(os.path.exists(f'Output/Img_report/{province}/{year}{d_m}')==False):
    os.makedirs(f'Output/Img_report/{province}/{year}{d_m}') 

Img_water=f'Output/Img_water/{province}/{year}{d_m}'
if(os.path.exists(Img_water)==False):
    os.makedirs(Img_water)
if(os.path.exists(f'Output/Img_water/{province}/{year}')==False):
    os.makedirs(f'Output/Img_water/{province}/{year}{d_m}') 

#%%      ẢNH TỶ SỐ
from PIL import Image, ImageDraw, ImageFont
bands = [("G", G),
         ("NIR", NIR),
         ("SWIR1", SWIR1)]

# Định nghĩa công thức
formula = "(NIR *SWIR1)/((G+epsilon)*(G+epsilon))"

# Gọi hàm calculate_custom_index với danh sách bands và công thức
ty_so = calculate_custom_index(bands, formula=formula,epsilon=1e-20)

show=np.copy(ty_so)
show[G==0]=None

plt.title(f' 18/05/2023')
plt.imshow(show, cmap='gray')
plt.axis('off')
plt.savefig(f'{Img_report}/Ảnh tỷ số NIR_G ({year}).png')
plt.show()
save_tiff(show,file_path=f'Output/Img_water/{province}/{year}/Ảnh tỷ số {year}.TIF',transform=image_info['transform'],crs=image_info['crs'])    
#%%
data=ty_so.reshape(-1, 1)
plt.hist(data, bins=100, density=True, alpha=0.5, color='b')
# Tùy chỉnh đồ thị
plt.title(f'Histogram {year}')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()
#%%                                  K-MEAN clustering
# Prepearing data

#%%
mu=1/10
image=1/(1+np.exp(-(data)**(mu)))**(1/2)
# plt.axis('off')
# plt.imshow(image, cmap='gray')
# plt.show()
# Assuming you have already defined 'image' and 'year' variables
X = image.reshape(-1, 1)
#%%
plt.hist(X, bins=5000, density=True, alpha=0.5, color='b')
plt.title(f'Histogram {year}')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()
plt.savefig(f'{Img_report}/Histogram ({year}).png')
#%%
plt.imshow(X.reshape(G.shape[0], G.shape[1]), cmap='gray')
plt.show()
#%%
save_tiff(X.reshape(G.shape[0], G.shape[1]),file_path=f'Output/Img_water/{province}/{year}{d_m}/Sigmod_10{year}.TIF',transform=image_info['transform'],crs=image_info['crs'])    
#%%
from scipy.spatial.distance import cdist
k_max=10
# tạo biểu đồ khuỷu tay
distortions = []
K = range(1,k_max)
for k in K:
    kmeanModel = KMeans(n_clusters=k,init='k-means++',n_init='auto', random_state=11).fit(X)
    distortions.append(sum(np.min(cdist(X, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / X.shape[0])

# Vẽ biểu đồ
plt.plot(K, distortions, 'bx-')
plt.xlabel('k')
plt.ylabel('Distortion')
plt.title('Number of cluster')
plt.show()
plt.savefig(f'{Img_report}/Elbow ({year}).png')
#%% Training 
n_cluster = 5
kmeans= KMeans(n_clusters=n_cluster,init='k-means++',n_init='auto', random_state=11)
labels = kmeans.fit_predict(X)
print(kmeans.inertia_)
#%%%
labels=np.copy(labels.reshape(G.shape[0], G.shape[1]))
#%%

plt.figure()
cmap = plt.cm.colors.ListedColormap(['red','blue', 'green','yellow','white','gray',
                                     'purple'])
plt.title(f'K-Means ({year})')
temp=np.copy(labels).astype(float)
temp[G==0]=np.nan
plt.axis('off')
plt.imshow(temp, cmap=cmap)
plt.show()
plt.savefig(f'{Img_report}/Ảnh phân loại K-Means({year}).png')
#%%
labels[G==0]=-9999
save_tiff(labels.reshape(G.shape[0], G.shape[1]),file_path=f'Output/Img_water/{province}/{year}{d_m}/K-Means_sigmod_10.TIF',transform=image_info['transform'],crs=image_info['crs'])

#%%
labels_water=2
clus=np.copy(labels).astype(float)
#clus[np.logical_or(clus == labels_water, clus == 0)] = 255

clus[(clus != labels_water)] = 0
clus[clus == labels_water] = 255
clus[G==0]=-9999

#%%%
plt.figure()

plt.title(f'Water ({year} )')
# Tạo colormap tùy chỉnh với hai màu: đen và deepskyblue
custom_cmap = mcolors.ListedColormap(['black', 'deepskyblue'])

# Tạo hình ảnh với colormap tùy chỉnh
plt.imshow(clus, cmap=custom_cmap)

# Tạo bảng chú thích
cbar = plt.colorbar(ticks=[1, 255],shrink=0.6) 
cbar.set_ticklabels(['No-Water', 'Water'])  
plt.axis('off')
plt.show()
plt.savefig(f'{Img_report}/Ảnh xác định nước ({year}).png')


#%%
save_tiff(clus,file_path=f'Output/Img_water/{province}/{year}{d_m}/Water {year}.TIF',transform=image_info['transform'],crs=image_info['crs'])

