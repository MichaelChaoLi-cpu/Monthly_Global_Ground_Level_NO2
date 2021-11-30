# -*- coding: utf-8 -*-
"""
Convert h4 to GeoTiff 

For DP_05 Nighttime Light from VNP46A3

Band: AllAngle_Composite_Snow_Free_Num 

Resolution: 15 arc sec ()

Created on Mon Nov 29 14:26:02 2021

@author: NASA, M.L.
"""

from osgeo import gdal
import os

## List input raster files
g = os.walk('D:/10_Article/09_TempOutput/08_MonthlyNighttimeLightTif/Temp')
h5FileList = []
for path, dir_list, file_list in g:
    for file_name in file_list:
        fileAddress = os.path.join(path, file_name).replace("\\","/")
        if fileAddress[-2:] == 'h5':
            if os.path.getsize(fileAddress)!=0:
                h5FileList.append(fileAddress)
                
loop = 0 

while (loop < 324): # len(h5FileList)):
    rasterFilePre = h5FileList[loop][70:-3]
    
    fileExtension = "_BBOX.tif"
    
    ## Open HDF file
    hdflayer = gdal.Open(h5FileList[loop], gdal.GA_ReadOnly)
    
    #print (hdflayer.GetSubDatasets())
    
    # Open raster layer
    #hdflayer.GetSubDatasets()[0][0] - for first layer
    #hdflayer.GetSubDatasets()[1][0] - for second layer ...etc
    subhdflayer = hdflayer.GetSubDatasets()[5][0]
    rlayer = gdal.Open(subhdflayer, gdal.GA_ReadOnly)
    #outputName = rlayer.GetMetadata_Dict()['long_name']
    
    #Subset the Long Name
    outputName = subhdflayer[97:]
    
    outputNameNoSpace = outputName.strip().replace(" ","_").replace("/","_")
    outputNameFinal = rasterFilePre + fileExtension
    print(outputNameFinal)
    
    outputFolder = "D:/10_Article/09_TempOutput/08_MonthlyNighttimeLightTif/OutputTif/"
    
    outputRaster = outputFolder + outputNameFinal
    
    #collect bounding box coordinates
    HorizontalTileNumber = int(rlayer.GetMetadata_Dict()["HorizontalTileNumber"])
    VerticalTileNumber = int(rlayer.GetMetadata_Dict()["VerticalTileNumber"])
        
    WestBoundCoord = (10*HorizontalTileNumber) - 180
    NorthBoundCoord = 90-(10*VerticalTileNumber)
    EastBoundCoord = WestBoundCoord + 10
    SouthBoundCoord = NorthBoundCoord - 10
    
    EPSG = "-a_srs EPSG:4326" #WGS84
    
    translateOptionText = EPSG+" -a_ullr " + str(WestBoundCoord) + " " + str(NorthBoundCoord) + " " + str(EastBoundCoord) + " " + str(SouthBoundCoord)
    
    translateoptions = gdal.TranslateOptions(gdal.ParseCommandLine(translateOptionText))
    gdal.Translate(outputRaster,rlayer, options=translateoptions)
    
    #Display image in QGIS (run it within QGIS python Console) - remove comment to display
    #iface.addRasterLayer(outputRaster, outputNameFinal)
    
    loop = loop + 1
    