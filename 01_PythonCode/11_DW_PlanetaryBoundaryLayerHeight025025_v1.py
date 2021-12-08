# -*- coding: utf-8 -*-
"""
Change the resolution from 5 arc seconds to 0.25 arc degree

Created on Wed Dec  8 13:16:42 2021

@author: li.chao.987@s.kyushu-u.ac.jp
"""

from osgeo import gdal
import glob
import os
import netCDF4
import numpy as np

src_dataset = gdal.Open("D:/10_Article/09_TempOutput/07_MonthlyPrecipitationTif/totalPrecipitationRate201501.tif", gdal.GA_ReadOnly)
geotransform = (-180.0, 0.625, 0.0, -90.0, 0.0, 0.5)
spatialreference = src_dataset.GetProjection()
ncol = 576
nrow = 361
nband = 1

FileList = glob.glob("D:/10_Article/09_TempOutput/10_PlanetaryBoundaryLayerHeight/Temp/*.nc4")
outputPBLHGeoTiffFileFolder = "D:/10_Article/09_TempOutput/10_PlanetaryBoundaryLayerHeight/RawData/"

for nc4File in FileList:
    #nc4File = FileList[0]
    readNc4File = netCDF4.Dataset(nc4File)
    
    PBLH = readNc4File["PBLH"][:]
    PBLH = np.nanmean(PBLH, axis = 0)
    PBLHRateOutputRaster = outputPBLHGeoTiffFileFolder + 'PBLH_' + nc4File[92:98] + ".tif"
    
    driver = gdal.GetDriverByName("GTiff")
    dst_dataset = driver.Create(PBLHRateOutputRaster, ncol, nrow, nband, gdal.GDT_Float32)
    dst_dataset.SetGeoTransform(geotransform)
    dst_dataset.SetProjection(spatialreference)
    dst_dataset.GetRasterBand(1).WriteArray(PBLH)
    dst_dataset = None
    
    
#Resampling to 0.25 bilinear 
FileList = glob.glob("D:/10_Article/09_TempOutput/10_PlanetaryBoundaryLayerHeight/RawData/*.tif")
outputPBLHGeoTiffFileFolder = "D:/10_Article/09_TempOutput/10_PlanetaryBoundaryLayerHeight/Resample/"

for GeoTiffFile in FileList:
    
    raster = gdal.Open(GeoTiffFile,  gdal.GA_ReadOnly)
    
    outputFinalAddress = outputPBLHGeoTiffFileFolder + GeoTiffFile[73:79] + ".2525.tif"
    rasterChangeResolution = gdal.Warp(outputFinalAddress, raster, xRes = 0.25,  yRes = 0.25,
                                       resampleAlg = gdal.gdalconst.GRA_Bilinear, 
                                       outputType = gdal.gdalconst.GDT_Float32)
    rasterChangeResolution = None
    raster = None
    
