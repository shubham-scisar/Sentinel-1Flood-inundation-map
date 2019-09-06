# Sentinel-1Flood-inundation-map
Contains a script that generates flood inundation map for a sentinel-1 SAR image

# Description
This repository contains a python file that creates a flood inundation map of a particular portion of a Sentinel-1 SAR Image.
Further, this script implements Lee filter for Speckle filtering. The paths of Image,annotation and calibration files have to be given by the user.
The script generates a flood inundation map file named "flood.tiff" in the parent directory.The user can also create a histogram and save it via matplotlib plot.


# Prerequisites
Python 3(3.5),GDAL(2.2.4),Numpy(1.12.1),matplotlib(2.2.2),scikit-image(0.13.0),scipy(1.1.0).


# More about the Author

Video link of my talk at PyCon India 2018,Hyderabad : https://www.youtube.com/watch?v=Yon0aJ4GBFA
