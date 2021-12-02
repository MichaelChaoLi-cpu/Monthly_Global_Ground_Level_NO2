# -*- coding: utf-8 -*-
"""
Merge the 0.25 arc degree GeoTiff to one

Created on Wed Dec  1 11:19:22 2021

@author: M.L.
"""

from osgeo import gdal
import glob


rasterMergedLocation = "D:\\10_Article\\09_TempOutput\\08_MonthlyNighttimeLightTif\\MergeTif\\"

year = 2016
month = ['01', '02', '03', '04', '05', '06',
         '07', '08', '09', '10', '11', '12']

i = 0
while i < 12:
    rasterFileList = glob.glob("D:/10_Article/09_TempOutput/08_MonthlyNighttimeLightTif/OutputTif/" + str(year) + month[i] + '*.tif')
    fileNameDate = str(year) + month[i]
    fileNamePre = "NTL"
    fileExtention = ".tif"
    outputFullName = rasterMergedLocation + fileNamePre  + fileNameDate + fileExtention
    print(outputFullName)
    vrt = gdal.BuildVRT("merged.vrt", rasterFileList)
    gdal.Translate(outputFullName, vrt)
    vrt = None
    i = i + 1