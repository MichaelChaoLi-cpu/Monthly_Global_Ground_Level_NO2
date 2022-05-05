# -*- coding: utf-8 -*-
"""
Convert h4 to GeoTiff 

For DP_10 OMI/Aura NO2 Total and Tropospheric Column Daily L2 from OMNO2G

Band: 0(totalNo2) 2(troposphericNo2) 29(terrainPressure) 

Resolution: 0.25 arc degree

Created on Wed May  4 16:23:17 2022

@author: li.chao.987@s.kyushu-u.ac.jp
"""

from osgeo import gdal
import os
import numpy as np
import glob

src_dataset = gdal.Open("D:/10_Article/08_MonthlyRaster/IDW_REM/200702.tif")
geotransform = (-180.0, 0.25, 0.0, -90.0, 0.0, 0.25)
spatialreference = src_dataset.GetProjection()
ncol = 1440
nrow = 720
nband = 1

## List input raster files
aimFolder = 'D:\\10_Article\\Test\\'
rasterFiles = glob.glob(aimFolder + "*.he5")
rasterFilesSelect = []
year = ['2015', '2016', '2017', '2018', '2019', '2020', '2021']
for raster in rasterFiles:
    if raster[39:43] in year:
        rasterFilesSelect.append(raster)
rasterFilesSelect = rasterFilesSelect[0:2479]
#print(rasterFiles)
#rasterFiles = rasterFiles[-61:]

troposphericNo2FileExtension = ".tif"
troposphericNo2OutputFolder = "D:/10_Article/09_TempOutput/16_DailyTroposphericNo2Tif/"

try:
    os.mkdir(troposphericNo2OutputFolder)
except:
    pass

for raster in rasterFiles:
    hdflayer = gdal.Open(raster, gdal.GA_ReadOnly)
    subhdflayer = hdflayer.GetSubDatasets()[9][0]
    rlayer = gdal.Open(subhdflayer, gdal.GA_ReadOnly)
    singleDayRasterArray = rlayer.ReadAsArray()
    singleDayRasterArray[singleDayRasterArray < 0] = np.nan
    singleDayRasterMeanArray = np.nanmean(singleDayRasterArray, axis = 0)
    
    rasterName = raster[39:48]
    troposphericNo2OutputNameFinal = troposphericNo2OutputFolder + \
        rasterName + troposphericNo2FileExtension
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(troposphericNo2OutputNameFinal, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(singleDayRasterMeanArray)
    dst_dataset = None
    
    rlayer = None
    hdflayer = None
    
    


