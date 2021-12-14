# Monthly_Global_Ground_Level_NO2 (DP10)
This article tries to estimate monthly global ground-level NO2 concentration from 2015 to 2021 based on NASA data.
Originally, we use GOME-2 data to detect the global ground-level NO2, from 2007 to 2017. The data used in the analysis are from 2015 to 2017.
  
## Author  
Chao Li, Shunsuke Managi
  
## Data
### Used in 01_DW_MonthlyAverageNo2Panel_v1.R
The measured ground-level NO2 concentration is provided by aqicn.org <https://aqicn.org/data-platform/register/>. Over 540 cities are covered, globally. The data are from 2015/01/01 to 2021/11/12 (the date when we downloaded the 2021 data). We use the monthly average values of the median concentration values of the listed cities. After we check the json file of the location <https://aqicn.org/data-platform/covid19/airquality-covid19-cities.json>, we find that in some cites there are several measurement points. In this way, the median value of the several measurements points in a specific city could represents the city's situation. Because most satellite products are monthly, the data have been convert into monthly data by averaging.  
  
### Used in 02_DW_MonthlyRasterDataPanel_v0.R
[PanelNo2Dataset.Rdata](03_Rawdata/PanelNo2Dataset.Rdata) "no2" monthly average no2 concentration ppm.  
**"mg_m2_total_no2"** monthly average total column amount of no2 (mg/m2). This data are extracted from the OMI/Aura NO2 Total and Tropospheric Column Daily L2 from OMNO2G, band 0(totalNo2). Resolution: 0.25 arc degree.   
**"mg_m2_troposphere_no2"** monthly average tropospheric column amount of no2 (mg/m2). This data are extracted from the OMI/Aura NO2 Total and Tropospheric Column Daily L2 from OMNO2G, band 2(troposphericNo2). Resolution: 0.25 arc degree.  
**"ter_pressure"** monthly average terrain surface pressure (hpa). This data are extracted from the OMI/Aura NO2 Total and Tropospheric Column Daily L2 from **OMNO2G**, band 29(terrainPressure). Resolution: 0.25 arc degree.
**"dayTimeTemperature"** monthly average day time temperature (C). This is raw data, which should be multiply by the factor (0.02 K). **MOD11C3** 0.05 arc degree monthly. We average value convert them into 0.25 arc degree.  
**"nightTimeTemperature"** monthly average night time temperature (C). This is raw data, which should be multiply by the factor (0.02 K). **MOD11C3** 0.05 arc degree monthly. We use average value to convert them into 0.25 arc degree.  
**"ndvi"** NDVI -1 to 1. This is raw data, which should be multiply by the factor (0.0001). **MOD13C2** 0.05 arc degree monthly NDVI. We use average value to convert them into 0.25 arc degree.  
**"humidity"** g/kg means 1 gram water in the 1 kg air. From GLDAS Noah Land Surface Model L4 from **GLDAS_NOAH025_M**, "Qair_f_tavg" band, their resolution is 0.25 arc degree.
**"precipitation"** the precipitation unit is kg/(m2 * h). From GLDAS Noah Land Surface Model L4 from **GLDAS_NOAH025_M**, "Rainf_f_tavg" band, their resolution is 0.25 arc degree.  
**"speedwind"** the wind speed unit is m/s. From GLDAS Noah Land Surface Model L4 from **GLDAS_NOAH025_M**, **"Wind_f_inst" band, their resolution is 0.25 arc degree.
**"NTL" nighttime light are from <https://eogdata.mines.edu/>. Its resolution is 15 arc-second, we use the average value convert it into 0.25 arc degree. 
**"PBLR"** planetary boundary layer height unit is m. The data are from M2TMNXFLX ("PBLH") with the resolution is 0.5 * 0.625. we use bilinear interpolation to get 0.25 arc-degree resolution.  
Additionally, there are several other data are gathered, but not used due to low correlation with the dependent variable, including **UVAerosolIndex**, **ozone**, **cloudfraction**, and **cloudpressure**.  

## Code
### Python Code
In this project, the python codes are mainly used to process the he5 or nc4 data, since its high speed with GDAL.  
[01_DW_TotalTroposhereNO2TerrianPressure_v0.py](01_PythonCode/01_DW_TotalTroposhereNO2TerrianPressure_v0.py): This script is to extract the data, including **"mg_m2_total_no2"**, **"mg_m2_troposphere_no2"**, and **"ter_pressure"**, from **OMNO2G**.  
[02_DW_PrecipitationRateSpecificHumidity_v0.py](01_PythonCode/02_DW_PrecipitationRateSpecificHumidity_v0.py): This script is to extract the data, including **"humidity"** and **"precipitation"**, from **GLDAS_NOAH025_M**.  
[03_DW_CityLocationJson_v0](01_PythonCode/03_DW_CityLocationJson_v0): This script is to location json file into csv.  
[04_DW_AddingMissingDataVNP46A3_v1](01_PythonCode/04_DW_AddingMissingDataVNP46A3_v1): This script is to get the VNP46A3. **(Aborted, Maybe, could be used in the nighttime light remote sensing project)**  
python ...laads-data-download.py -s <https://ladsweb.modaps.eosdis.nasa.gov/archive/Science%20Domain/Atmosphere/Atmospheric%20Profiles/MODIS%20Terra%20C6.1%20-%20Temperature%20and%20Water%20Vapor%20Profiles%205-Min%20L2%20Swath%205km/2021/> -t TOKEN -d Dictionary
This code in the cmd should miss some file. The 04_DW python is designed to solve this problem.  
[05_DW_ConvertH5intoTiffVPN46A3_v1.py](01_PythonCode/05_DW_ConvertH5intoTiffVPN46A3_v1.py): This script is to convert the he5 into GeoTiff by GDAL.  
[06_DW_mergeTifVPN46A3_v0.py](01_PythonCode/06_DW_mergeTifVPN46A3_v0.py): This script is to merge the GeoTiff generated by [05_DW_ConvertH5intoTiffVPN46A3_v1.py](01_PythonCode/05_DW_ConvertH5intoTiffVPN46A3_v1.py).  
[07_DW_DownloadNighttimeLightFromEogdata_v0.py](01_PythonCode/07_DW_DownloadNighttimeLightFromEogdata_v0.py): This a spider script to grasp the NTL data. We do not know other batch downloading method, so using the **selenium** to automatically download.  
[08_DW_ChangeResolutionTo025_v1.py](01_PythonCode/08_DW_ChangeResolutionTo025_v1.py): This script is to change the resolution to 0.25 arc degree.  
[09_DW_MergeTiffsTo025_v1.py](01_PythonCode/09_DW_MergeTiffsTo025_v1.py): This script is to merge the GeoTiff from [08_DW_ChangeResolutionTo025_v1.py](01_PythonCode/08_DW_ChangeResolutionTo025_v1.py).  
[10_DW_WindSpeed025025_v1.py](01_PythonCode/10_DW_WindSpeed025025_v1.py): This script is to extract the data, **speedwind**, from **GLDAS_NOAH025_M**.  
[11_DW_PlanetaryBoundaryLayerHeight025025_v1.py](01_PythonCode/11_DW_PlanetaryBoundaryLayerHeight025025_v1.py): This script is to extract the data, **PBLH**, from **M2TMNXFLX**.  
[12_DW_CloudFractionCloudPressure_v1.py](01_PythonCode/12_DW_CloudFractionCloudPressure_v1.py): This script is to extract the data, includind **CloudFraction** and **CloudPressure**, from **OMNO2G**.  
[13_DW_AerosolOpticalDepth_v1.py](01_PythonCode/13_DW_AerosolOpticalDepth_v1.py): This script is to extract the data, **AOD**, from **OMAEROe**.  
[14_DW_MonthlyOzone_v1.py](01_PythonCode/14_DW_MonthlyOzone_v1.py): This script is to extract the data, **Ozone**, from **OMTO3G**.  
[15_DW_MonthlyUVAerosalIndex_v1.py](01_PythonCode/15_DW_MonthlyUVAerosalIndex_v1.py): This script is to extract the data, **UVAerosolIndex**, from **OMTO3G**.  
  
## Data Collecting and Preprocessing
This section mainly records the satellite data sources.  

