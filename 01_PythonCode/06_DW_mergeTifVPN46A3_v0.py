# -*- coding: utf-8 -*-
"""
Merge the Tiff from 10 python 05_DW into monthly version 

Created on Mon Nov 29 14:37:52 2021

@author: M.L.
"""

from osgeo import gdal
import os
import subprocess
import glob

mergeFolder = 'D:/10_Article/09_TempOutput/08_MonthlyNighttimeLightTif/MergeTif/'

os.chdir('D:/10_Article/09_TempOutput/08_MonthlyNighttimeLightTif/OutputTif')
fileList = glob.glob("VNP46A3.*.tif")


## List input raster files
g = os.walk('D:/10_Article/09_TempOutput/08_MonthlyNighttimeLightTif/OutputTif')
rasterFileList = []
for path, dir_list, file_list in g:
    for file_name in file_list:
        fileAddress = os.path.join(path, file_name).replace("\\","/")
        if fileAddress[-3:] == 'tif':
            if os.path.getsize(fileAddress)!=0:
                rasterFileList.append(fileAddress)

files_string = " ".join(rasterFileList[1:5])
print(files_string)

fileNameDate = rasterFileList[1][75:82]
fileNamePre = "NTL"
fileExtention = "t.tif"
outputFullName = mergeFolder + fileNamePre  + fileNameDate + fileExtention
vrt = gdal.BuildVRT("merged.vrt", fileList)
gdal.Translate(outputFullName, vrt)
vrt = None
