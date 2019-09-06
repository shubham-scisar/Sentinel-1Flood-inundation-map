# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import math
import gdal,ogr,os,osr
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import uniform_filter
from scipy.ndimage.measurements import variance
from skimage.filters import threshold_otsu
import os


"""
Opening SAR Image

"""

image=gdal.Open('/home/shubham/S1A_IW_GRDH_1SDV_20170829T002620_20170829T002645_018131_01E74D_D734.SAFE/measurement'+
                  '/s1a-iw-grd-vv-20170829t002620-20170829t002645-018131-01e74d-001.tiff')


"""
 Selecting a subset as selected in SNAP
 
"""
m1=4325      # start pixel row    
n1=6228       # start pixel column  

      
m2=7612       # end pixel row
n2=11072      # end pixel column




#Opening image which is by default an amplitude image
 
amp_vv=((image.ReadAsArray())[m1:m2+1,n1:n2+1]).astype('float64')

row,col=amp_vv.shape


# Extracting projection information

proj=image.GetProjection()


srs=osr.SpatialReference(wkt=proj)
if srs.IsProjected:
    projection=srs.GetAttrValue('AUTHORITY',1)
    




 # Calibration to sigma naught
 
from xml.dom import minidom

pp=list()
file="/home/shubham/S1A_IW_GRDH_1SDV_20170829T002620_20170829T002645_018131_01E74D_D734.SAFE/annotation/calibration/calibration-s1a-iw-grd-vv-20170829t002620-20170829t002645-018131-01e74d-001.xml"
xmldoc = minidom.parse(file)
itemlist = xmldoc.getElementsByTagName('line') 
sigma_n=xmldoc.getElementsByTagName('sigmaNought') 

for i in range(len(itemlist)):
            pp.append(int(itemlist[i].firstChild.nodeValue))
 
lin_int=pp[1]-pp[0]

def calib(m1,n1):
    
    
    global n2,col
    
    for i in range(len(pp)):
        if pp[i]<=m1<=pp[i+1]:
            index=i  
        

    tt1="".join(sigma_n[index].firstChild.nodeValue.split())   
    
    ls_sigma=[]
    for i in range(0,len(tt1),12):
         ls_sigma.append(float(tt1[i:i+12]))
        
    index_sigma_1=n1//40
    ind1=n1%40
    ar=np.linspace(ls_sigma[index_sigma_1],ls_sigma[index_sigma_1+1],num=40,endpoint=True)
    sigma_val_1=ar[ind1]



    index_sigma_2=n2//40
    ind2=n2%40
    ar2=np.linspace(ls_sigma[index_sigma_2],ls_sigma[index_sigma_2+1],num=40,endpoint=True)
    sigma_val_2=ar2[ind2]

    arr=np.linspace(sigma_val_1,sigma_val_2,num=col,endpoint=True)

    
    arr=np.array([arr]).T
    arr=arr.T

    return arr




c=0
rem=col%lin_int

for i in range(0,row,lin_int+1):
    
       
    if c==0:
        arr=calib(i+m1,n1) 
        arr=np.ones((lin_int,1))*arr
        c=1
        
    else :
        arr1=calib(i+m1,n1) 
        arr1=np.ones((lin_int,1))*arr1
        arr=np.concatenate((arr,arr1),axis=0)
        


sigma_gain=arr[0:row:]       # Array containing gains
        
sigma_naught=(amp_vv/sigma_gain)**2    # Sigma-naught array
 
"""
Application of Speckle Filter -Lee   
"""


def lee(image, size):
    image_mean = uniform_filter(image, (size, size))
    image_sqr_mean = uniform_filter(image**2, (size, size))
    image_var = image_sqr_mean - image_mean**2

    overall_var = variance(image)

    image_weights = image_var / (image_var + overall_var)
    image_dspk = image_mean + image_weights * (image - image_mean)
    return image_dspk    
       
       
dspk_img=lee(sigma_naught,7)

        
                
dspk_img_log=10*np.log10(dspk_img)        # Creation of Log scaled Image

     
"""           
#Plot despeckled Image in log scale  
plt.imshow(dspk_img_log,interpolation='nearest',cmap='gray')      
"""


#Histogram Plot
plt.hist(dspk_img_log.ravel(),bins=512)    # Histogram Plot to find threshold




"""
# Generation of Mask"

"""



#After checking the log10 histogram we get the middle value between two maximas as -11.98

threshold=threshold_otsu(dspk_img_log)



dspk_img_mask=dspk_img_log

dspk_img_mask=np.where(dspk_img_mask<threshold,1,np.nan)


plt.imshow(dspk_img_mask,interpolation='nearest',cmap='Blues_r')









# Georeferencing using annotation file which contains gcp's


file="/home/shubham/S1A_IW_GRDH_1SDV_20170829T002620_20170829T002645_018131_01E74D_D734.SAFE/annotation/s1a-iw-grd-vv-20170829t002620-20170829t002645-018131-01e74d-001.xml"
xmldoc = minidom.parse(file)
linelist = xmldoc.getElementsByTagName('line') 
pixelist = xmldoc.getElementsByTagName('pixel') 
latlist =  xmldoc.getElementsByTagName('latitude') 
lonlist =  xmldoc.getElementsByTagName('longitude') 

linl=list()
pixl=list()
latl=list()
lonl=list()


for i in range(len(linelist)):
            linl.append(int(linelist[i].firstChild.nodeValue))
            
            

for i in range(len(pixelist)):
            pixl.append(int(pixelist[i].firstChild.nodeValue)) 
           
            
            
for i in range(len(latlist)):
            latl.append(float(latlist[i].firstChild.nodeValue))
            
            
for i in range(len(lonlist)):
            lonl.append(float(lonlist[i].firstChild.nodeValue))
            

linl=np.array(linl).reshape((len(linl),1))
pixl=np.array(pixl).reshape((len(pixl),1))

latl=np.array(latl).reshape((len(latl),1))
lonl=np.array(lonl).reshape((len(lonl),1))


  
geo=np.concatenate((linl,pixl,latl,lonl),axis=1)



l_size=len(np.unique(linl))
p_size=len(np.unique(pixl))

rowsize=np.unique(linl)[-1]
colsize=np.unique(pixl)[-1]     
 
# lat-lon values of the four corners of the scene
  
ullat=geo[0,2]                         #top left
ullon=geo[0,3]

urlat=geo[p_size-1,2]                  # top right
urlon=geo[p_size-1,3]



aa=np.where(geo==rowsize)[0][0]

lllat=geo[aa,2]                         #lower left
lllon=geo[aa,3]



lrlat=geo[geo.shape[0]-1,2]                 #lower right
lrlon=geo[geo.shape[0]-1,3]






lat_range=np.linspace(ullat,lllat,num=rowsize+1,endpoint=True)            #linear interpolation along lat and lon
lon_range=np.linspace(ullon,urlon,num=colsize+1,endpoint=True)




top_lat=lat_range[m1].tolist()         #upper left Corner latitude of scene
top_lon=lon_range[n1].tolist()         #upper left Corner Longitude of scene


top_lat=round(top_lat,2)
top_lon=round(top_lon,2)

# creating intermediate tiff file

def array2raster(file_name,origin,width,height,array,projection):
    

    cols = array.shape[1]
    rows = array.shape[0]
    originX = origin[0]
    originY = origin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(file_name, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, width, 0, originY, 0, height))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(projection)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
    


file_name='out.tiff'                           # intermediate tiff file
origin=(top_lon,top_lat)
width=n2-n1+1
height=m2-m1+1
projection=int(projection)

array2raster(file_name,origin,width,height,dspk_img_mask,projection)



bottom_rlat=lat_range[m2].tolist()             #bottom right latitude of scene
bottom_rlon=lon_range[n2].tolist()           #bottom right longitude of scene


bottom_rlat=round(bottom_rlat,2)
bottom_rlon=round(bottom_rlon,2)

outf="flood.tiff"                                              #final output file

 
# georeferencing using top left and bottom right corner using system call gdal_translate utility

os.system("gdal_translate -of GTiff -a_ullr %s %s %s %s -a_srs EPSG:%s %s %s"%(top_lon,top_lat,bottom_rlon,bottom_rlat,projection,file_name,outf))
os.remove("out.tiff")



