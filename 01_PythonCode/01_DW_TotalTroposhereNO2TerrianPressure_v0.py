# -*- coding: utf-8 -*-
"""
Convert h4 to GeoTiff 

For DP_10 OMI/Aura NO2 Total and Tropospheric Column Daily L2 from OMNO2G

Band: 5(total no2) 9(troposphere) 29(terrian pressure) 

Resolution: 0.25 arc degree

Created on Mon Nov 22 15:45:11 2021

@author: C.L.
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
totalNo2FileExtension = "_WGS84_.25.25_monthly_total_no2.tif"
troposphericNo2FileExtension = "_WGS84_.25.25_monthly_tropospheric_no2.tif"
terrainPressureFileExtension = "_WGS84_.25.25_monthly_terrain_pressure.tif"
totalNo2OutputFolder = "D:/10_Article/09_TempOutput/01_MonthlyTotalNo2Tif/"
troposphericNo2OutputFolder = "D:/10_Article/09_TempOutput/02_MonthlyTroposphericNo2Tif/"
terrainPressureOutputFolder = "D:/10_Article/09_TempOutput/03_MonthlyTerrainPressureTif/"

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
    totalNo2SingleMonthRasterList = []
    troposphericNo2SingleMonthRasterList = []
    terrainPressureSingleMonthRasterList = []
    for raster in singleMonthRasterNameList:
        hdflayer = gdal.Open(raster, gdal.GA_ReadOnly)
        subhdflayer = hdflayer.GetSubDatasets()[5][0]
        rlayer = gdal.Open(subhdflayer, gdal.GA_ReadOnly)
        singleDayRasterArray = rlayer.ReadAsArray()
        singleDayRasterArray[singleDayRasterArray < 0] = np.nan
        singleDayRasterMeanArray = np.nanmean(singleDayRasterArray, axis = 0)
        totalNo2SingleMonthRasterList.append(singleDayRasterMeanArray)
        subhdflayer = hdflayer.GetSubDatasets()[9][0]
        rlayer = gdal.Open(subhdflayer, gdal.GA_ReadOnly)
        singleDayRasterArray = rlayer.ReadAsArray()
        singleDayRasterArray[singleDayRasterArray < 0] = np.nan
        singleDayRasterMeanArray = np.nanmean(singleDayRasterArray, axis = 0)
        troposphericNo2SingleMonthRasterList.append(singleDayRasterMeanArray)
        subhdflayer = hdflayer.GetSubDatasets()[29][0]
        rlayer = gdal.Open(subhdflayer, gdal.GA_ReadOnly)
        singleDayRasterArray = rlayer.ReadAsArray()
        singleDayRasterArray[singleDayRasterArray < 0] = np.nan
        singleDayRasterMeanArray = np.nanmean(singleDayRasterArray, axis = 0)
        terrainPressureSingleMonthRasterList.append(singleDayRasterMeanArray)
        
    totalNo2MonthlyMeanArray = np.nanmean(totalNo2SingleMonthRasterList, axis = 0)
    troposphericNo2MonthlyMeanArray = np.nanmean(troposphericNo2SingleMonthRasterList, axis = 0)
    terrainPressureMonthlyMeanArray = np.nanmean(terrainPressureSingleMonthRasterList, axis = 0)
    
    #build the name and address
    totalNo2OutputNameFinal = rasterFilePre + aimedDate + totalNo2FileExtension   
    troposphericNo2OutputNameFinal = rasterFilePre + aimedDate + troposphericNo2FileExtension  
    terrainPressureOutputNameFinal = rasterFilePre + aimedDate + terrainPressureFileExtension  
    #print(totalNo2OutputNameFinal)  
    totalNo2OutputRaster = totalNo2OutputFolder + totalNo2OutputNameFinal
    troposphericNo2OutputRaster = troposphericNo2OutputFolder + troposphericNo2OutputNameFinal
    terrainPressureOutputRaster = terrainPressureOutputFolder + terrainPressureOutputNameFinal
    
    # write the tif files
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(totalNo2OutputRaster, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(totalNo2MonthlyMeanArray)
    dst_dataset = None
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(troposphericNo2OutputRaster, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(troposphericNo2MonthlyMeanArray)
    dst_dataset = None
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(terrainPressureOutputRaster, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(terrainPressureMonthlyMeanArray)
    dst_dataset = None
     
        
    