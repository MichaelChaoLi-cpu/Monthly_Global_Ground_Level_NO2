# -*- coding: utf-8 -*-
"""
Adding missing data

For DP_05 Nighttime Light from VNP46A3

Created on Sun Nov 28 14:18:56 2021

@author: M.L.
"""

import os


nasaAddress = "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000/VNP46A3/"
## List input raster files
g = os.walk('D:/10_Article/09_TempOutput/08_MonthlyNighttimeLightTif/Temp')
h5FileList = []
for path, dir_list, file_list in g:
    for file_name in file_list:
        fileAddress = os.path.join(path, file_name).replace("\\","/")
        if os.path.getsize(fileAddress)==0:
            fileAddress = fileAddress[61:]
            fileAddress = nasaAddress + fileAddress
            h5FileList.append(fileAddress)
        
writingFile = open('D:/10_Article/09_TempOutput/08_MonthlyNighttimeLightTif/Temp/addingMissing/missingLog.txt', "w")
for line in h5FileList:
    writingFile.writelines(line + "\n") 
writingFile.close()


"""
run cmd.exe as the manager
cd D:aimFolder
NUL > .crs_cookies

wget --load-cookies D:\10_Article\09_TempOutput\08_MonthlyNighttimeLightTif\Temp\addingMissing\.urs_cookies 
--save-cookies D:\10_Article\09_TempOutput\08_MonthlyNighttimeLightTif\Temp\addingMissing\.urs_cookies 
--auth-no-challenge=on --keep-session-cookies --user=chaoli0394 --password=097680Li 
--content-disposition -i D:\10_Article\09_TempOutput\08_MonthlyNighttimeLightTif\Temp\addingMissing\missingLog.txt
"""