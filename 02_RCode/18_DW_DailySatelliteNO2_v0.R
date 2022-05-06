# Author: M.L.

# end

library(tidyverse)
library(dplyr)
library(stringr)
library(raster)
library(doSNOW)
library(foreach)
library(sp)

extractPointDataFromRasterPara <- function(RasterFolder, filelist, cityLocationSpatialPoint,
                                           year_start_location, month_start_location, flip_reverse = T,
                                           aimed_column_name = "raw", year_end_location = year_start_location + 3,
                                           month_end_location = month_start_location + 1, core = 1
){
  RasterDataset <- 
    data.frame(Doubles=double(),
               Ints=integer(),
               Factors=factor(),
               Logicals=logical(),
               Characters=character(),
               stringsAsFactors=FALSE)
  cl <- makeSOCKcluster(core)
  registerDoSNOW(cl)
  getDoParWorkers()
  RasterDataset <- foreach (filename = filelist, .combine = base::rbind,
                            .packages='tidyverse') %dopar% {
                              test_tiff <- raster::raster(paste0(RasterFolder, filename))
                              if(flip_reverse){
                                test_tiff <- raster::flip(test_tiff, direction = 'y')
                              }
                              raster::crs(test_tiff) <- proj
                              Year <- stringr::str_sub(filename, year_start_location, year_end_location) %>% as.numeric()
                              Month <- stringr::str_sub(filename, month_start_location, month_end_location) %>% as.numeric()
                              
                              data_ext <- raster::extract(test_tiff, cityLocationSpatialPoint)
                              cityLocationSpatialPoint@data$raw <- data_ext
                              monthly_data <- cityLocationSpatialPoint@data %>%
                                dplyr::select(CityCode, raw)
                              monthly_data <- monthly_data %>%
                                dplyr::mutate(year = Year,
                                              month = Month)
                              monthly_data <- monthly_data %>% na.omit()
                            }
  stopCluster(cl)
  colnames(RasterDataset) <- c("CityCode", aimed_column_name, "year", "month")
  return(RasterDataset)
}

load("03_Rawdata/cityLocationSpatialPoint.Rdata")
cityLocationSpatialPoint.todaily <- cityLocationSpatialPoint
cityLocationSpatialPoint.todaily@data <- cityLocationSpatialPoint.todaily@data %>%
  dplyr::select(CityCode)
proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

#### extract daily tropospheric NO2 
troposphereNo2RasterFolder <- "D:/10_Article/09_TempOutput/16_DailyTroposphericNo2Tif/"
filelist <- list.files(troposphereNo2RasterFolder)
troposphereNo2RasterDataset <- 
  extractPointDataFromRasterPara(troposphereNo2RasterFolder, filelist, cityLocationSpatialPoint.todaily,
                                 year_start_location = 1, month_start_location = 6, T, "raw_no2",
                                 month_end_location = 9, core = 8)
troposphereNo2RasterDataset$day <- troposphereNo2RasterDataset$month %% 100
troposphereNo2RasterDataset$month <-floor(troposphereNo2RasterDataset$month / 100)

setwd("D:/10_Article/01_RawData/11_AirPollutionRawTable/")

#2015.01 - 2021.1123 global AQI: https://aqicn.org/data-platform/covid19/verify/57c94d81-c396-481a-a12c-93a0065c705f
filelist <- list.files()
filelist <- filelist[2:length(filelist)] #the first file always the reminder.

totalAirPollutionDataset <- data.frame(Doubles=double(),
                                       Ints=integer(),
                                       Factors=factor(),
                                       Logicals=logical(),
                                       Characters=character(),
                                       stringsAsFactors=FALSE)
for (singlefile in filelist){
  Global_AQI <- read.csv(file = singlefile, skip = 4, encoding = "UTF-8")
  Global_AQI$year <- str_sub(Global_AQI$Date, 1, 4) %>% as.numeric()
  Global_AQI$month <- str_sub(Global_AQI$Date, 6, 7) %>% as.numeric()
  Global_AQI$day <- str_sub(Global_AQI$Date, 9, 10) %>% as.numeric()
  totalAirPollutionDataset <- rbind(totalAirPollutionDataset, Global_AQI)
  rm(Global_AQI)
}
totalNo2Dataset <- totalAirPollutionDataset %>%
  filter(Specie == "no2") %>% 
  filter(median < 500) %>% 
  dplyr::select(Country, City, median, year, month, day) 

cityCode.list <- cityLocationSpatialPoint@data %>% dplyr::select(CityCode, City, Country)
totalNo2Dataset <- left_join(totalNo2Dataset, cityCode.list, by = c("Country", "City"))
totalNo2Dataset <- totalNo2Dataset %>% filter(!is.na(CityCode)) 
totalNo2Dataset <- totalNo2Dataset %>% dplyr::select(CityCode, median, year, month, day) 
### dataset from csv

mergedDailyDataset <- full_join(troposphereNo2RasterDataset, totalNo2Dataset,
                                 by = c("CityCode", "year", "month", "day"))

overlapDailyDataset <- mergedDailyDataset %>%
  filter(!is.na(raw_no2), !is.na(median))
overlapDailyDataset <- overlapDailyDataset %>%
  aggregate(by = list(overlapDailyDataset$CityCode, overlapDailyDataset$year, 
                      overlapDailyDataset$month), FUN = mean)
overlapDailyDataset <- overlapDailyDataset %>%
  dplyr::select("CityCode", "year", "month", "median", "raw_no2")

setwd("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub")

# convert molecular / cm2 to ug / m2
mol_g = 6.022140857 * 10^23  # mol
overlapDailyDataset$g_cm2 <- overlapDailyDataset$raw_no2 / mol_g * 46.0055 # convert mol to g
overlapDailyDataset$mg_m2_troposphere_no2_overlap <- overlapDailyDataset$g_cm2 * 10000 * 1000 # conver /cm2 to /m2 and g to mg
overlapDailyDataset$ug_m2_troposphere_no2_overlap <- overlapDailyDataset$mg_m2_troposphere_no2_overlap * 1000
overlapDailyDataset <- overlapDailyDataset %>% 
  dplyr::select("CityCode", "year", "month", "median", "ug_m2_troposphere_no2_overlap")
overlapDailyDataset <- overlapDailyDataset %>% 
  rename(no2_overlap = median)

differenceTestDataset <- inner_join(overlapDailyDataset, 
                                    usedDataset %>% dplyr::select("CityCode", "year", "month", "no2", "ug_m2_troposphere_no2"), 
                                    by = c("CityCode", "year", "month"))
RMSE_no2 <- sqrt(mean((differenceTestDataset$no2 - differenceTestDataset$no2_overlap)^2))
MAE <- mean(abs(differenceTestDataset$no2 - differenceTestDataset$no2_overlap))
mean(differenceTestDataset$no2)
R2_no2 <- sum((differenceTestDataset$no2 - differenceTestDataset$no2_overlap)^2)/
  sum((differenceTestDataset$no2 - mean(differenceTestDataset$no2))^2)
cor.test(differenceTestDataset$no2, differenceTestDataset$no2_overlap)
## r 0.9802

R2_ug_m2_troposphere_no2 <- 1 - sum((differenceTestDataset$ug_m2_troposphere_no2 - differenceTestDataset$ug_m2_troposphere_no2_overlap)^2)/
  sum((differenceTestDataset$ug_m2_troposphere_no2 - mean(differenceTestDataset$ug_m2_troposphere_no2))^2)
cor.test(differenceTestDataset$ug_m2_troposphere_no2, differenceTestDataset$ug_m2_troposphere_no2_overlap)
## r 0.9618

### check Puebla 2021 sept 
Pueble.daily <- totalAirPollutionDataset %>%
  filter(Country == 'MX', 
         City == "Puebla",
         year == 2021) %>%
  filter(Specie == "no2")
 