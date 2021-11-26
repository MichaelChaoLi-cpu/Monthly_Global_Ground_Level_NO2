# -*- coding: utf-8 -*-
"""
Convert nc4 to GeoTiff 

For DP_10 GLDAS Noah Land Surface Model L4 from GLDAS_NOAH025_M

Band: "Rainf_f_tavg" "Qair_f_tavg"

ReadMe file: https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1/summary

Resolution: 0.25 arc degree

Created on Thu Nov 25 16:48:36 2021

@author: li.chao.987@s.kyushu-u.ac.jp
"""

import numpy as np
import os
from scipy.io import netcdf
import netCDF4
from osgeo import gdal

src_dataset = gdal.Open("D:/10_Article/08_MonthlyRaster/IDW_REM/200702.tif")
geotransform = (-180.0, 0.25, 0.0, -60.0, 0.0, 0.25)
spatialreference = src_dataset.GetProjection()
ncol = 1440
nrow = 600
nband = 1

inputNc4FileFolder = 'D:/10_Article/09_TempOutput/07_MonthlyPrecipitationTif/Temp'
outputPrecipitationGeoTiffFileFolder = 'D:/10_Article/09_TempOutput/07_MonthlyPrecipitationTif/'
outputSpecificHumidityGeoTiffFileFolder = "D:/10_Article/09_TempOutput/06_MonthlyVaporTif/"

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
    
    totalPrecipitationRate = readNc4File["Rainf_f_tavg"][:] 
    totalPrecipitationRate = np.nanmean(totalPrecipitationRate, axis = 0)
    totalPrecipitationRateOutputRaster = outputPrecipitationGeoTiffFileFolder + 'totalPrecipitationRate' + raster[17:23] + ".tif"
    
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(totalPrecipitationRateOutputRaster, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(totalPrecipitationRate)
    dst_dataset = None
    
    specificHumidity = readNc4File["Qair_f_inst"][:]
    specificHumidity = np.nanmean(specificHumidity, axis = 0)
    specificHumidityOutputRaster = outputSpecificHumidityGeoTiffFileFolder + 'specificHumidity' + raster[17:23] + ".tif"
    
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(specificHumidityOutputRaster, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(specificHumidity)
    dst_dataset = None