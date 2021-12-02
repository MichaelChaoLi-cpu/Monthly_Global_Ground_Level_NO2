# -*- coding: utf-8 -*-
"""
Change the resolution from 5 arc seconds to 0.25 arc degree

Created on Tue Nov 30 13:58:30 2021

@author: li.chao.987@s.kyushu-u.ac.jp
"""

import os, sys, tarfile
from osgeo import gdal
import glob

def extract(tar_url, extract_path='.'):
    print(tar_url)
    tar = tarfile.open(tar_url, 'r')
    for item in tar:
        tar.extract(item, extract_path)
        tar.close()
        break

rasterChangeResolutionLocation = "D:\\10_Article\\09_TempOutput\\08_MonthlyNighttimeLightTif\\OutputTif\\"
         
tgzFilesLocation =  "D:\\10_Article\\09_TempOutput\\08_MonthlyNighttimeLightTif\\download\\*.tgz"
tgzFileList = glob.glob("D:\\10_Article\\09_TempOutput\\08_MonthlyNighttimeLightTif\\download\\*.tgz")

for tar_location in tgzFileList:
    extract_path = "D:\\10_Article\\09_TempOutput\\08_MonthlyNighttimeLightTif\\download\\Temp"
    extract(tar_location, extract_path)
    
    
fileList = glob.glob(extract_path + "\\*.tif")
for file in fileList:
    raster = gdal.Open(file)
    outputFinalAddress = rasterChangeResolutionLocation + file[80:105] + ".2525.tif"
    rasterChangeResolution = gdal.Warp(outputFinalAddress, raster, xRes = 0.25,  yRes = 0.25,
                                       resampleAlg = gdal.GRA_Average)
    rasterChangeResolution = None
    raster = None


