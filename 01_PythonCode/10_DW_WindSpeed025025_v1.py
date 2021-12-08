# -*- coding: utf-8 -*-
"""
Convert nc4 to GeoTiff 

For DP_10 GLDAS Noah Land Surface Model L4 from GLDAS_NOAH025_M

Band: "Wind_f_inst" 

ReadMe file: https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary

Resolution: 0.25 arc degree

Created on Wed Dec  8 11:25:02 2021

@author: M.L.
"""

import numpy as np
import os
import netCDF4
from osgeo import gdal

src_dataset = gdal.Open("D:/10_Article/08_MonthlyRaster/IDW_REM/200702.tif")
geotransform = (-180.0, 0.25, 0.0, -60.0, 0.0, 0.25)
spatialreference = src_dataset.GetProjection()
ncol = 1440
nrow = 600
nband = 1

inputNc4FileFolder = 'D:/10_Article/09_TempOutput/07_MonthlyPrecipitationTif/Temp'

outputWindSpeedGeoTiffFileFolder = 'D:/10_Article/09_TempOutput/09_WindSpeed/'

## List input raster files
os.chdir(inputNc4FileFolder)
rasterFilesRaw = os.listdir(os.getcwd())
rasterFiles = []
for raster in rasterFilesRaw:
    if raster[-3:] == "nc4":
        rasterFiles.append(raster)
#print(rasterFiles)

for raster in rasterFiles:
    nc4File = inputNc4FileFolder + '/' + raster
    readNc4File = netCDF4.Dataset(nc4File)
    
    totalWindSpeedRate = readNc4File["Wind_f_inst"][:] 
    totalWindSpeedRate = np.nanmean(totalWindSpeedRate, axis = 0)
    totalWindSpeedRateOutputRaster = outputWindSpeedGeoTiffFileFolder + 'totalWindSpeed' + raster[17:23] + ".tif"
    
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(totalWindSpeedRateOutputRaster, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(totalWindSpeedRate)
    dst_dataset = None
    