# -*- coding: utf-8 -*-
"""
Convert h4 to GeoTiff 

For DP_10 OMI/Aura NO2 Total and Tropospheric Column Daily L2 from OMNO2G

Band: 5(total no2) 9(troposphere) 29(terrian pressure) 

Resolution: 0.25 arc degree

Created on Wed Dec  8 15:40:56 2021

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
os.chdir('D:/10_Article/Test')
rasterFilesRaw = os.listdir(os.getcwd())
rasterFiles = []
for raster in rasterFilesRaw:
    if raster[-3:] == "he5":
        rasterFiles.append(raster)
#print(rasterFiles)

rasterFilePre = rasterFiles[1][:20]
cloudFractionFileExtension = "_WGS84_.25.25_monthly_cloud_fraction.tif"
cloudPressureFileExtension = "_WGS84_.25.25_monthly_cloud_pressure.tif"
cloudFractionOutputFolder = "D:/10_Article/09_TempOutput/11_MonthlyCloudFraction/"
cloudPressureOutputFolder = "D:/10_Article/09_TempOutput/12_MonthlyCloudPressure/"

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
    cloudFractionSingleMonthRasterList = []
    cloudPressureSingleMonthRasterList = []
    for raster in singleMonthRasterNameList:
        #raster = singleMonthRasterNameList[1]
        hdflayer = gdal.Open(raster, gdal.GA_ReadOnly)
        subhdflayer = hdflayer.GetSubDatasets()[0][0]
        rlayer = gdal.Open(subhdflayer, gdal.GA_ReadOnly)
        singleDayRasterArray = rlayer.ReadAsArray()
        singleDayRasterArray[singleDayRasterArray < 0] = np.nan
        singleDayRasterMeanArray = np.nanmean(singleDayRasterArray, axis = 0)
        cloudFractionSingleMonthRasterList.append(singleDayRasterMeanArray)
        subhdflayer = hdflayer.GetSubDatasets()[2][0]
        rlayer = gdal.Open(subhdflayer, gdal.GA_ReadOnly)
        singleDayRasterArray = rlayer.ReadAsArray()
        singleDayRasterArray[singleDayRasterArray < 0] = np.nan
        singleDayRasterMeanArray = np.nanmean(singleDayRasterArray, axis = 0)
        cloudPressureSingleMonthRasterList.append(singleDayRasterMeanArray)
        
    cloudFractionMonthlyMeanArray = np.nanmean(cloudFractionSingleMonthRasterList, axis = 0)
    cloudPressureMonthlyMeanArray = np.nanmean(cloudPressureSingleMonthRasterList, axis = 0)
    
    #build the name and address
    cloudFractionOutputNameFinal = rasterFilePre + aimedDate + cloudFractionFileExtension   
    cloudPressureOutputNameFinal = rasterFilePre + aimedDate + cloudPressureFileExtension  
    
    #print(totalNo2OutputNameFinal)  
    cloudFractionOutputRaster = cloudFractionOutputFolder + cloudFractionOutputNameFinal
    cloudPressureOutputRaster = cloudPressureOutputFolder + cloudPressureOutputNameFinal
    
    # write the tif files
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(cloudFractionOutputRaster, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(cloudFractionMonthlyMeanArray)
    dst_dataset = None
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(cloudPressureOutputRaster, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(cloudPressureMonthlyMeanArray)
    dst_dataset = None
