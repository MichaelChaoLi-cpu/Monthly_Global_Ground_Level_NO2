# -*- coding: utf-8 -*-
"""
Extension of windspeed raster, because some cities are close to the sea which have no data.

# Avoid using the buffer to calculate the raw data, we make this.

Created on Tue Dec 14 15:27:50 2021

@author: M.L.
"""

import numpy as np
from osgeo import gdal
import glob

src_dataset = gdal.Open("D:/10_Article/09_TempOutput/09_WindSpeed/totalWindSpeed201501.tif", gdal.GA_ReadOnly)
geotransform = src_dataset.GetGeoTransform()
spatialreference = src_dataset.GetProjection()
ncol = 1440
nrow = 600
nband = 1

### wind speed
originalFileLocation = "D:/10_Article/09_TempOutput/09_WindSpeed"
outputFileLocation = "D:/10_Article/09_TempOutput/09_WindSpeed/Add025Outline/"

fileList = glob.glob(originalFileLocation + "/*.tif")

for file in fileList:
    raster = gdal.Open(file, gdal.GA_ReadOnly)
    rasterArray = raster.ReadAsArray()
    
    addTimes = 0
    while addTimes < 2:
        addRasterLayer = np.full((600, 1440), np.nan)
        for i  in np.linspace(1,598,598, dtype = int):
            for j  in np.linspace(1,1438,1438, dtype = int):
                if np.isnan(rasterArray[i, j]): #addRasterLayer[316, 95] is number
                    grid = np.nanmean(rasterArray[i-1:i+2,j-1:j+2])
                    if np.isnan(grid):
                        pass
                    else:
                        addRasterLayer[i, j] = grid
        rasterArray = np.array([rasterArray, addRasterLayer])
        rasterArray = np.nanmean(rasterArray, axis = 0)
        addTimes = addTimes + 1         
    
    outputFileName = "add_" + file[41:-4] + ".tif"
    finalOutputName = outputFileLocation + outputFileName
    
    # write the tif files
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(finalOutputName, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(rasterArray)
    dst_dataset = None

### precipitation 
originalFileLocation = "D:/10_Article/09_TempOutput/07_MonthlyPrecipitationTif"
outputFileLocation = "D:/10_Article/09_TempOutput/07_MonthlyPrecipitationTif/Add025Outline/"

fileList = glob.glob(originalFileLocation + "/*.tif")

for file in fileList:
    raster = gdal.Open(file, gdal.GA_ReadOnly)
    rasterArray = raster.ReadAsArray()
    
    addTimes = 0
    while addTimes < 2:
        addRasterLayer = np.full((600, 1440), np.nan)
        for i  in np.linspace(1,598,598, dtype = int):
            for j  in np.linspace(1,1438,1438, dtype = int):
                if np.isnan(rasterArray[i, j]): #addRasterLayer[316, 95] is number
                    grid = np.nanmean(rasterArray[i-1:i+2,j-1:j+2])
                    if np.isnan(grid):
                        pass
                    else:
                        addRasterLayer[i, j] = grid
        rasterArray = np.array([rasterArray, addRasterLayer])
        rasterArray = np.nanmean(rasterArray, axis = 0)
        addTimes = addTimes + 1 
    
    outputFileName = "add_" + file[55:-4] + ".tif"
    finalOutputName = outputFileLocation + outputFileName
    
    # write the tif files
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(finalOutputName, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(rasterArray)
    dst_dataset = None
    
### humidity
originalFileLocation = "D:/10_Article/09_TempOutput/06_MonthlyVaporTif"
outputFileLocation = "D:/10_Article/09_TempOutput/06_MonthlyVaporTif/Add025Outline/"

fileList = glob.glob(originalFileLocation + "/*.tif")

for file in fileList:
    raster = gdal.Open(file, gdal.GA_ReadOnly)
    rasterArray = raster.ReadAsArray()
    
    addTimes = 0
    while addTimes < 2:
        addRasterLayer = np.full((600, 1440), np.nan)
        for i  in np.linspace(1,598,598, dtype = int):
            for j  in np.linspace(1,1438,1438, dtype = int):
                if np.isnan(rasterArray[i, j]): #addRasterLayer[316, 95] is number
                    grid = np.nanmean(rasterArray[i-1:i+2,j-1:j+2])
                    if np.isnan(grid):
                        pass
                    else:
                        addRasterLayer[i, j] = grid
        rasterArray = np.array([rasterArray, addRasterLayer])
        rasterArray = np.nanmean(rasterArray, axis = 0)
        addTimes = addTimes + 1 
    
    outputFileName = "add_" + file[47:-4] + ".tif"
    finalOutputName = outputFileLocation + outputFileName
    
    # write the tif files
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(finalOutputName, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(rasterArray)
    dst_dataset = None