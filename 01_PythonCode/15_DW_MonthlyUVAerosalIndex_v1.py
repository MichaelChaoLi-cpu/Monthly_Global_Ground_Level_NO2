# -*- coding: utf-8 -*-
"""
for OMI/Aura Ozone (O3) Total Column Daily L2 Global Gridded 0.25 degree x 0.25 degree V3 (OMTO3G)

Band: 28 (UVAerosolIndex)

Resolution: 0.25 arc degree

Created on Thu Dec  9 11:32:06 2021

@author: M.L.
"""

from osgeo import gdal
import os
import numpy as np

src_dataset = gdal.Open("D:/10_Article/08_MonthlyRaster/IDW_REM/200702.tif")
geotransform = (-180.0, 0.25, 0.0, -90.0, 0.0, 0.25)
spatialreference = src_dataset.GetProjection()
ncol = 1440
nrow = 720
nband = 1

## List input raster files
os.chdir('D:/10_Article/09_TempOutput/14_MonthlyOzone/RawData')
rasterFilesRaw = os.listdir(os.getcwd())
rasterFiles = []
for raster in rasterFilesRaw:
    if raster[-3:] == "he5":
        rasterFiles.append(raster)
#print(rasterFiles)

rasterFilePre = rasterFiles[1][:20]
UVAerosolIndexFileExtension = "_WGS84_.25.25_monthly_UVAerosolIndex.tif"
UVAerosolIndexOutputFolder = "D:/10_Article/09_TempOutput/15_MonthlyUVAerosolIndex/"

#get the month and year in the raster names
yearMonthList = []
for raster in rasterFiles:
    if raster[20:27] in yearMonthList:
        pass
    else:
        yearMonthList.append(raster[20:27])
#print(yearMonthList)

for aimedDate in yearMonthList:
    #aimedDate = yearMonthList[1] #TestCode
    singleMonthRasterNameList = []
    for raster in rasterFiles:
        if raster[20:27] == aimedDate:
            singleMonthRasterNameList.append(raster)
    UVAerosolIndexSingleMonthRasterList = []
    for raster in singleMonthRasterNameList:
        #raster = singleMonthRasterNameList[1]
        hdflayer = gdal.Open(raster, gdal.GA_ReadOnly)
        #print(hdflayer.GetSubDatasets()[2])
        subhdflayer = hdflayer.GetSubDatasets()[28][0]
        rlayer = gdal.Open(subhdflayer, gdal.GA_ReadOnly)
        singleDayRasterArray = rlayer.ReadAsArray()
        singleDayRasterArray[singleDayRasterArray < 0] = np.nan
        singleDayRasterMeanArray = np.nanmean(singleDayRasterArray, axis = 0)
        UVAerosolIndexSingleMonthRasterList.append(singleDayRasterMeanArray)
        
    UVAerosolIndexMonthlyMeanArray = np.nanmean(UVAerosolIndexSingleMonthRasterList, axis = 0)
    
    #build the name and address
    UVAerosolIndexOutputNameFinal = rasterFilePre + aimedDate + UVAerosolIndexFileExtension
    UVAerosolIndexOutputRaster = UVAerosolIndexOutputFolder + UVAerosolIndexOutputNameFinal
    
    # write the tif files
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(UVAerosolIndexOutputRaster, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(UVAerosolIndexMonthlyMeanArray)
    dst_dataset = None