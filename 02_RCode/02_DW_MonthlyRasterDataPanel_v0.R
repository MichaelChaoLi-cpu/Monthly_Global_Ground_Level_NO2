# Author: M.L.

# input: PanelNo2Dataset.Rdata
# PanelNo2Dataset.Rdata: "no2" monthly average no2 concentration ppm.

# input: CityLocationOfficial.csv
# CityLocationOfficial.csv: "Country", "City", "Latitude", "Longitude"

# output: mergedDataset.Rdata
# mergedDataset.Rdata "no2" monthly average no2 concentration ppm.
# mergedDataset.Rdata "mg_m2_total_no2" monthly average total column amount of no2 (mg/m2)
# mergedDataset.Rdata "mg_m2_troposphere_no2" monthly average tropospheric column amount of no2 (mg/m2)
# mergedDataset.Rdata "ter_pressure" monthly average terrain surface pressure (hpa)
# mergedDataset.Rdata "dayTimeTemperature" monthly average day time temperature (C)
# mergedDataset.Rdata "nightTimeTemperature" monthly average night time temperature (C)
# mergedDataset.Rdata "ndvi" NDVI -1 to 1
# mergedDataset.Rdata "humidity" g/kg means 1 gram water in the 1 kg air.
# mergedDataset.Rdata "precipitation" the precipitation unit is kg/(m2 * h)
# mergedDataset.Rdata "speedwind" the wind speed unit is m/s
# mergedDataset.Rdata "NTL" nighttime light
# mergedDataset.Rdata "PBLR" planetary boundary layer height unit is m.
# mergedDataset.Rdata "CityCode" identity index.
# mergedDataset.Rdata "period" year * 100 + month, time index. 

# end

library(tidyverse)
library(dplyr)
library(raster)
library(sp)
library(stringr)
library(rgdal)
library(lubridate)

extractPointDataFromRaster <- function(RasterFolder, filelist, cityLocationSpatialPoint,
                                       year_start_location, month_start_location, flip_reverse = T,
                                       aimed_column_name = "raw", year_end_location = year_start_location + 3,
                                       month_end_location = month_start_location + 1
                                       ){
  RasterDataset <- 
    data.frame(Doubles=double(),
               Ints=integer(),
               Factors=factor(),
               Logicals=logical(),
               Characters=character(),
               stringsAsFactors=FALSE)
  for (filename in filelist){
    test_tiff <- raster::raster(paste0(RasterFolder, filename))
    if(flip_reverse){
      test_tiff <- flip(test_tiff, direction = 'y')
    }
    crs(test_tiff) <- proj
    Year <- str_sub(filename, year_start_location, year_end_location) %>% as.numeric()
    Month <- str_sub(filename, month_start_location, month_end_location) %>% as.numeric()
    
    data_ext <- raster::extract(test_tiff, cityLocationSpatialPoint)
    cityLocationSpatialPoint@data$raw <- data_ext
    monthly_data <- cityLocationSpatialPoint@data %>%
      dplyr::select(Country, City, raw)
    monthly_data <- monthly_data %>%
      mutate(year = Year,
             month = Month)
    RasterDataset <- rbind(RasterDataset, monthly_data)
  }
  colnames(RasterDataset) <- c("Country", "City", aimed_column_name, "year", "month")
  return(RasterDataset)
}

extractBufferDataFromRaster <- function(RasterFolder, filelist, cityLocationSpatialPoint,
                                       year_start_location, month_start_location, flip_reverse = T,
                                       aimed_column_name = "raw", year_end_location = year_start_location + 3,
                                       month_end_location = month_start_location + 1
){
  RasterDataset <- 
    data.frame(Doubles=double(),
               Ints=integer(),
               Factors=factor(),
               Logicals=logical(),
               Characters=character(),
               stringsAsFactors=FALSE)
  for (filename in filelist){
    test_tiff <- raster::raster(paste0(RasterFolder, filename))
    if(flip_reverse){
      test_tiff <- flip(test_tiff, direction = 'y')
    }
    crs(test_tiff) <- proj
    Year <- str_sub(filename, year_start_location, year_end_location) %>% as.numeric()
    Month <- str_sub(filename, month_start_location, month_end_location) %>% as.numeric()
    
    data_ext <- raster::extract(test_tiff, cityLocationSpatialPoint, fun = mean, na.rm = TRUE)
    cityLocationSpatialPoint@data$raw <- data_ext
    monthly_data <- cityLocationSpatialPoint@data %>%
      dplyr::select(Country, City, raw)
    monthly_data <- monthly_data %>%
      mutate(year = Year,
             month = Month)
    RasterDataset <- rbind(RasterDataset, monthly_data)
  }
  colnames(RasterDataset) <- c("Country", "City", aimed_column_name, "year", "month")
  return(RasterDataset)
}

cityLocation <- read.csv("D:/10_Article/01_RawData/12_LocationJson/CityLocationOfficial.csv",
                     encoding="UTF-8") %>%
  dplyr::select(X0, X1, X2, X3)
colnames(cityLocation) <- c("Latitude", "Longitude", "City", "Country")
cityLocation <- cityLocation %>% 
  mutate(City = ifelse(City == "Washington, D.C.", "Washington D.C.", City))

load("03_Rawdata/PanelNo2Dataset.Rdata")
nametest <- totalNo2Dataset %>% dplyr::select(Country, City) %>% unique()
cityLocation <- left_join(nametest, cityLocation, by = c("Country", "City"))
rm(nametest)
cityLocation$CityCode <- 1:nrow(cityLocation)
totalNo2Dataset <- left_join(totalNo2Dataset, cityLocation, by = c("Country", "City"))

proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
xy <- cityLocation %>% dplyr::select(Longitude, Latitude)
cityLocationSpatialPoint <- SpatialPointsDataFrame(coords = xy, data = cityLocation[,c(1, 2, 3, 4, 5)],
                                                   proj4string = CRS(proj))
rm(xy)

#get monthly total no2 from the OMNO2G, band 5(total no2)
totalNo2RasterFolder <- "D:/10_Article/09_TempOutput/01_MonthlyTotalNo2Tif/"
filelist <- list.files(totalNo2RasterFolder)
totalNo2RasterDataset <- 
  extractPointDataFromRaster(totalNo2RasterFolder, filelist, cityLocationSpatialPoint,
                             21, 26, T, "raw_no2")
# convert molecular / cm2 to ug / m2
mol_g = 6.022140857 * 10^23  # mol
totalNo2RasterDataset$g_cm2 <- totalNo2RasterDataset$raw_no2 / mol_g * 46.0055 # convert mol to g
totalNo2RasterDataset$mg_m2_total_no2 <- totalNo2RasterDataset$g_cm2 * 10000 * 1000 # conver /cm2 to /m2 and g to mg
totalNo2RasterDataset <- totalNo2RasterDataset %>% dplyr::select("Country", "City", "year", "month", "mg_m2_total_no2")
totalNo2RasterDataset$Date <- as.Date(
  paste0(as.character(totalNo2RasterDataset$year),"-",as.character(totalNo2RasterDataset$month),"-01")
  )
totalNo2RasterDataset.lag <- totalNo2RasterDataset
totalNo2RasterDataset.lag$Date <- totalNo2RasterDataset.lag$Date %m+% months(1)
totalNo2RasterDataset.lag <- totalNo2RasterDataset.lag %>% 
  dplyr::select(Country, City, mg_m2_total_no2, Date) %>%
  rename(mg_m2_total_no2_lag = mg_m2_total_no2)
totalNo2RasterDataset <- left_join(totalNo2RasterDataset, totalNo2RasterDataset.lag,
                                   by = c("Country", "City", "Date"))
totalNo2RasterDataset <- totalNo2RasterDataset %>%
  dplyr::select(-raw_no2, -g_cm2, -Date)

#get monthly troposphere no2 from the OMNO2G, band 9(troposphere no2)
troposphereNo2RasterFolder <- "D:/10_Article/09_TempOutput/02_MonthlyTroposphericNo2Tif/"
filelist <- list.files(troposphereNo2RasterFolder)
troposphereNo2RasterDataset <- 
  extractPointDataFromRaster(troposphereNo2RasterFolder, filelist, cityLocationSpatialPoint,
                             21, 26, T, "raw_no2")
# convert molecular / cm2 to ug / m2
mol_g = 6.022140857 * 10^23  # mol
troposphereNo2RasterDataset$g_cm2 <- troposphereNo2RasterDataset$raw_no2 / mol_g * 46.0055 # convert mol to g
troposphereNo2RasterDataset$mg_m2_troposphere_no2 <- troposphereNo2RasterDataset$g_cm2 * 10000 * 1000 # conver /cm2 to /m2 and g to mg
troposphereNo2RasterDataset <- troposphereNo2RasterDataset %>% dplyr::select("Country", "City", "year", "month", "mg_m2_troposphere_no2")
troposphereNo2RasterDataset$Date <- as.Date(
  paste0(as.character(troposphereNo2RasterDataset$year),"-",as.character(troposphereNo2RasterDataset$month),"-01")
)
troposphereNo2RasterDataset.lag <- troposphereNo2RasterDataset
troposphereNo2RasterDataset.lag$Date <- troposphereNo2RasterDataset.lag$Date %m+% months(1)
troposphereNo2RasterDataset.lag <- troposphereNo2RasterDataset.lag %>% 
  dplyr::select(Country, City, mg_m2_troposphere_no2, Date) %>%
  rename(mg_m2_troposphere_no2_lag = mg_m2_troposphere_no2)
troposphereNo2RasterDataset <- left_join(troposphereNo2RasterDataset, troposphereNo2RasterDataset.lag,
                                   by = c("Country", "City", "Date"))
troposphereNo2RasterDataset <- troposphereNo2RasterDataset %>%
  dplyr::select(-Date)

#get monthly terrain pressure from the OMNO2G, band 29 (terrain pressure)
terrainPressureRasterFolder <- "D:/10_Article/09_TempOutput/03_MonthlyTerrainPressureTif/"
filelist <- list.files(terrainPressureRasterFolder)
terrainPressureRasterDataset <- 
  extractPointDataFromRaster(terrainPressureRasterFolder, filelist, cityLocationSpatialPoint,
                             21, 26, T, "ter_pressure")

#get monthly Daytime temperature from the MOD11C3 0.25 arc degree
dayTimeTemperatureRasterFolder <- "D:/10_Article/09_TempOutput/04_MonthlyTemperatureTif/Surf_Temp_Monthly_005dg_v6/LST_Day_CMG/"
filelist <- list.files(dayTimeTemperatureRasterFolder)
dayTimeTemperatureRasterDataset <- 
  extractPointDataFromRaster(dayTimeTemperatureRasterFolder, filelist, cityLocationSpatialPoint,
                             21, month_start_location = 26, F,
                             "dayTimeTemperature", month_end_location = 28)
dayTimeTemperatureRasterDataset <- aggregate(dayTimeTemperatureRasterDataset$dayTimeTemperature,
                                             by = list(dayTimeTemperatureRasterDataset$Country, 
                                                       dayTimeTemperatureRasterDataset$City,
                                                       dayTimeTemperatureRasterDataset$year, 
                                                       dayTimeTemperatureRasterDataset$month), 
                                             FUN = "mean", na.rm = T
                                             )
colnames(dayTimeTemperatureRasterDataset) <- c("Country", "City", "year", "month", "dayTimeTemperature")
dayTimeTemperatureRasterDataset$date <- 
  as.Date((dayTimeTemperatureRasterDataset$month - 1),
          origin = paste0(dayTimeTemperatureRasterDataset$year,"-01-01")) %>% as.character()
dayTimeTemperatureRasterDataset$month <- str_sub(dayTimeTemperatureRasterDataset$date, 6, 7) %>% as.numeric()
dayTimeTemperatureRasterDataset <- dayTimeTemperatureRasterDataset %>% dplyr::select(-date)
dayTimeTemperatureRasterDataset$dayTimeTemperature <- dayTimeTemperatureRasterDataset$dayTimeTemperature *
  0.02 - 273.16  #convert into c degree temperature

#get monthly Nighttime temperature from the MOD11C3 0.25 arc degree
nightTimeTemperatureRasterFolder <- "D:/10_Article/09_TempOutput/04_MonthlyTemperatureTif/Surf_Temp_Monthly_005dg_v6/LST_Night_CMG/"
filelist <- list.files(nightTimeTemperatureRasterFolder)
nightTimeTemperatureRasterDataset <- 
  extractPointDataFromRaster(nightTimeTemperatureRasterFolder, filelist, cityLocationSpatialPoint,
                             23, month_start_location = 28, F,
                             "nightTimeTemperature", month_end_location = 30)
nightTimeTemperatureRasterDataset <- aggregate(nightTimeTemperatureRasterDataset$nightTimeTemperature,
                                             by = list(nightTimeTemperatureRasterDataset$Country, 
                                                       nightTimeTemperatureRasterDataset$City,
                                                       nightTimeTemperatureRasterDataset$year, 
                                                       nightTimeTemperatureRasterDataset$month), 
                                             FUN = "mean", na.rm = T
)
colnames(nightTimeTemperatureRasterDataset) <- c("Country", "City", "year", "month", "nightTimeTemperature")
nightTimeTemperatureRasterDataset$date <- 
  as.Date((nightTimeTemperatureRasterDataset$month - 1),
          origin = paste0(nightTimeTemperatureRasterDataset$year,"-01-01")) %>% as.character()
nightTimeTemperatureRasterDataset$month <- str_sub(nightTimeTemperatureRasterDataset$date, 6, 7) %>% as.numeric()
nightTimeTemperatureRasterDataset <- nightTimeTemperatureRasterDataset %>% dplyr::select(-date)
nightTimeTemperatureRasterDataset$nightTimeTemperature <- nightTimeTemperatureRasterDataset$nightTimeTemperature *
  0.02 - 273.16  #convert into c degree temperature

#get monthly NDVI temperature from the MOD13C3 0.25 arc degree
ndviRasterFolder <- "D:/10_Article/09_TempOutput/05_MonthlyNDVITif/VI_Monthly_005dg_v6/NDVI/"
filelist <- list.files(ndviRasterFolder)
ndviRasterDataset <- 
  extractPointDataFromRaster(ndviRasterFolder, filelist, cityLocationSpatialPoint,
                             14, month_start_location = 19, F,
                             "ndvi", month_end_location = 21)
ndviRasterDataset <- aggregate(ndviRasterDataset$ndvi,
                                               by = list(ndviRasterDataset$Country, 
                                                         ndviRasterDataset$City,
                                                         ndviRasterDataset$year, 
                                                         ndviRasterDataset$month), 
                                               FUN = "mean", na.rm = T
)
colnames(ndviRasterDataset) <- c("Country", "City", "year", "month", "ndvi")
ndviRasterDataset$date <- 
  as.Date((ndviRasterDataset$month - 1),
          origin = paste0(ndviRasterDataset$year,"-01-01")) %>% as.character()
ndviRasterDataset$month <- str_sub(ndviRasterDataset$date, 6, 7) %>% as.numeric()
ndviRasterDataset <- ndviRasterDataset %>% dplyr::select(-date)
ndviRasterDataset$ndvi <- ndviRasterDataset$ndvi / 10000  #convert into from 1 to -1 

#get monthly water vapor from the GLDAS_NOAH025_M 0.25 arc degree
cityLocationSpatialBuffer <- rgeos::gBuffer(cityLocationSpatialPoint, byid = T, width = 0.5)
humidityRasterFolder <- "D:/10_Article/09_TempOutput/06_MonthlyVaporTif/"
filelist <- list.files(humidityRasterFolder)
humidityRasterDataset <- 
  extractBufferDataFromRaster(humidityRasterFolder, filelist, cityLocationSpatialBuffer,
                             17, 21, T, "humidity")
humidityRasterDataset$humidity <- humidityRasterDataset$humidity * 1000 #convert the unit into g/kg
humidityRasterDataset <- humidityRasterDataset %>% as.data.frame()
# 1 g/kg means 1 gram water in the 1 kg air.

#get monthly precipitation from the GLDAS_NOAH025_M 0.25 arc degree
precipitationRasterFolder <- "D:/10_Article/09_TempOutput/07_MonthlyPrecipitationTif/"
filelist <- list.files(precipitationRasterFolder)
filelist <- filelist[2:length(filelist)]
precipitationRasterDataset <- 
  extractBufferDataFromRaster(precipitationRasterFolder, filelist, cityLocationSpatialBuffer,
                             23, 27, T, "precipitation")
precipitationRasterDataset$precipitation <- precipitationRasterDataset$precipitation * 3600 
precipitationRasterDataset <- precipitationRasterDataset %>% as.data.frame()
# now, the precipitation unit is kg/(m2 * h)  

#get monthly speed wind from the GLDAS_NOAH025_M 0.25 arc degree
speedWindRasterFolder <- "D:/10_Article/09_TempOutput/09_WindSpeed/"
filelist <- list.files(speedWindRasterFolder)
speedWindRasterDataset <- 
  extractBufferDataFromRaster(speedWindRasterFolder, filelist, cityLocationSpatialBuffer,
                              15, 19, T, "speedwind")
speedWindRasterDataset <- speedWindRasterDataset %>% as.data.frame()
# now, the speed wind unit is m/s  

#get monthly NTL 0.25 arc degree
NTLRasterFolder <- "D:/10_Article/09_TempOutput/08_MonthlyNighttimeLightTif/MergeTif/"
filelist <- list.files(NTLRasterFolder)
NTLRasterDataset <- 
  extractPointDataFromRaster(NTLRasterFolder, filelist, cityLocationSpatialPoint,
                             4, 8, F, "NTL")
#get monthly NTL 0.25 arc degree

#get monthly planetary boundary layer height, 0.25 * 0.25 
PBLHRasterLayer <- "D:/10_Article/09_TempOutput/10_PlanetaryBoundaryLayerHeight/Resample/"
filelist <- list.files(PBLHRasterLayer)
PBLHRasterDataset <- 
  extractPointDataFromRaster(PBLHRasterLayer, filelist, cityLocationSpatialPoint,
                             1, 5, F, "PBLH")
#get monthly planetary boundary layer height, 0.25 * 0.25 

# get monthly cloud fraction, 0.25 * 0.25
cloudFractionRasterLayer <- "D:/10_Article/09_TempOutput/11_MonthlyCloudFraction/"
filelist <- list.files(cloudFractionRasterLayer)
cloudFractionRasterDataset <- 
  extractPointDataFromRaster(cloudFractionRasterLayer, filelist, cityLocationSpatialPoint,
                             21, 26, T, "cloudfraction")
# get monthly cloud fraction, 0.25 * 0.25

# get monthly cloud pressure, 0.25 * 0.25
cloudPressureRasterLayer <- "D:/10_Article/09_TempOutput/12_MonthlyCloudPressure/"
filelist <- list.files(cloudPressureRasterLayer)
cloudPressureRasterDataset <- 
  extractPointDataFromRaster(cloudPressureRasterLayer, filelist, cityLocationSpatialPoint,
                             21, 26, T, "cloudpressure")
# get monthly cloud pressure, 0.25 * 0.25

# get monthly Ozone, 0.25 * 0.25
ozoneRasterLayer <- "D:/10_Article/09_TempOutput/14_MonthlyOzone/"
filelist <- list.files(ozoneRasterLayer)
filelist <- filelist[1:77]
ozoneRasterDataset <- 
  extractPointDataFromRaster(ozoneRasterLayer, filelist, cityLocationSpatialPoint,
                             21, 26, T, "ozone")
# get monthly cloud pressure, 0.25 * 0.25

# get monthly UV Aerosol Index, 0.25 * 0.25
UVAerosolIndexRasterLayer <- "D:/10_Article/09_TempOutput/15_MonthlyUVAerosolIndex/"
filelist <- list.files(UVAerosolIndexRasterLayer)
UVAerosolIndexRasterDataset <- 
  extractPointDataFromRaster(UVAerosolIndexRasterLayer, filelist, cityLocationSpatialPoint,
                             21, 26, T, "UVAerosolIndex")
# get monthly UV Aerosol Index, 0.25 * 0.25

#break point
mergedDataset <- left_join(totalNo2Dataset, totalNo2RasterDataset, by = c("Country", "City", "year", "month"))
mergedDataset <- left_join(mergedDataset, troposphereNo2RasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$mg_m2_total_no2, mergedDataset$no2)
cor.test(mergedDataset$mg_m2_troposphere_no2, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, terrainPressureRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$ter_pressure, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, dayTimeTemperatureRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$dayTimeTemperature, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, nightTimeTemperatureRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$nightTimeTemperature, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, ndviRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$ndvi, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, humidityRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$humidity, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, precipitationRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$precipitation, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, NTLRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$NTL, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, speedWindRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$speedwind, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, PBLHRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$PBLH, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, cloudFractionRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$cloudfraction, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, cloudPressureRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$cloudpressure, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, ozoneRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$ozone, mergedDataset$no2)
mergedDataset <- left_join(mergedDataset, UVAerosolIndexRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(mergedDataset$UVAerosolIndex, mergedDataset$no2)

# convert into ug/m3
Pcoef = 0.00750061683
MW = 46.0055
mergedDataset$no2_measured_mg.m3 <-
  Pcoef * mergedDataset$ter_pressure * MW * mergedDataset$no2 /
  (62.4 * (273.2 + mergedDataset$dayTimeTemperature/2 + mergedDataset$nightTimeTemperature/2))

mergedDataset$month <- mergedDataset$month %>% as.factor()
mergedDataset$year <- mergedDataset$year %>% as.factor()
mergedDataset <- mergedDataset %>%
  mutate(
    Y2016 = ifelse(year == "2016", 1, 0),
    Y2017 = ifelse(year == "2017", 1, 0),
    Y2018 = ifelse(year == "2018", 1, 0),
    Y2019 = ifelse(year == "2019", 1, 0),
    Y2020 = ifelse(year == "2020", 1, 0),
    Y2021 = ifelse(year == "2021", 1, 0)
  )

mergedDataset %>% summary()

na.test <- mergedDataset %>% na.omit()
na.test$count <- 1
na.test <- aggregate(na.test$count, by = list(na.test$City, na.test$Country), FUN=sum)

save(mergedDataset, file = "C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/03_Rawdata/mergedDataset.Rdata")
