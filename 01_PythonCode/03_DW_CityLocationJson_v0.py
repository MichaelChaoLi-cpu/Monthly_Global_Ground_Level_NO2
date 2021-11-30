# -*- coding: utf-8 -*-
"""
from location json to csv

json from: https://aqicn.org/data-platform/covid19/airquality-covid19-cities.json

Created on Fri Nov 26 13:47:58 2021

@author: M.L.
"""

import json
import pandas as pd

jsonLocation = "D:/10_Article/01_RawData/12_LocationJson/airquality-covid19-cities.json"

jsonFile = open(jsonLocation, 'r', encoding='UTF-8')
airPollutionMeasurementLocation = json.load(jsonFile)

locationTable = []

i = 0
while i < len(airPollutionMeasurementLocation['data']):
    placeDict = airPollutionMeasurementLocation['data'][i]['Place']
    singleLine = [placeDict['geo'][0], placeDict['geo'][1],
                  placeDict['name'], placeDict['country'],
                  placeDict['feature'], placeDict['pop']]
    locationTable.append(singleLine)
    i = i + 1    

locationTable = pd.DataFrame(locationTable)
locationTable.to_csv("D:/10_Article/01_RawData/12_LocationJson/CityLocationOfficial.csv",
                     encoding = 'utf-8')
