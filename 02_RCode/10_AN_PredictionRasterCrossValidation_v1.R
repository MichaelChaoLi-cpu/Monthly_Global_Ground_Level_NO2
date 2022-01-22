# Author: M.L.

# input: COEF_raster.RData
# COEF_raster.RData: "CO_ug_m2_troposphere_no2.kriged.raster" interpolation of troposheric no2 coefficient
#                                                             based on ordinary kriging
# COEF_raster.RData: "CO_ndvi.kriged.raster" interpolation of ndvi coefficient based on ordinary kriging (OK)
# COEF_raster.RData: "CO_temp.kriged.raster" interpolation of night temperature coefficient (OK)
# COEF_raster.RData: "CO_PBLH.kriged.raster" interpolation of PBLH coefficient (OK)
# COEF_raster.RData: "CO_precipitation.kriged.raster" interpolation of precipitation coefficient (OK)
# COEF_raster.RData: "CO_ter_pressure.kriged.raster" interpolation of air pressure coefficient (OK)
# COEF_raster.RData: "CO_Y2016.kriged.raster" interpolation of 2016 dummy variable coefficient (OK)
# COEF_raster.RData: "CO_Y2017.kriged.raster" interpolation of 2017 dummy variable coefficient (OK)
# COEF_raster.RData: "CO_Y2018.kriged.raster" interpolation of 2018 dummy variable coefficient (OK)
# COEF_raster.RData: "CO_Y2019.kriged.raster" interpolation of 2019 dummy variable coefficient (OK)
# COEF_raster.RData: "CO_Y2020.kriged.raster" interpolation of 2020 dummy variable coefficient (OK)
# COEF_raster.RData: "CO_Y2021.kriged.raster" interpolation of 2021 dummy variable coefficient (OK)

# output: MEAN_raster.RData
# MEAN_raster.RData: "MEAN_ug_m2_troposphere_no2.kriged.raster" interpolation of troposheric no2 mean value
#                                                               in the usedDataset based on ordinary kriging (OK)
# MEAN_raster.RData: "MEAN_ndvi.kriged.raster" interpolation of ndvi mean value (OK)
# MEAN_raster.RData: "MEAN_temp.kriged.raster" interpolation of temperature mean value (OK)
# MEAN_raster.RData: "MEAN_NTL.kriged.raster" interpolation of night time light mean value (OK)
# MEAN_raster.RData: "MEAN_PBLH.kriged.raster" interpolation of PBLH mean value (OK)
# MEAN_raster.RData: "MEAN_precipitation.kriged.raster" interpolation of precipitation mean value (OK)
# MEAN_raster.RData: "MEAN_ter_pressure.kriged.raster" interpolation of air pressure mean value (OK)

# end

library(tidyverse)
library(dplyr)
library(raster)
library(lubridate)

makePredictRaster <- function(aim_year, aim_month){
  load("05_CoefficientRaster/COEF_raster.RData")
  load("05_CoefficientRaster/MEAN_raster.RData")
  
  raster_troposphere_no2_folder <- 
    "D:/10_Article/09_TempOutput/02_MonthlyTroposphericNo2Tif/OMI-Aura_L2G-OMNO2G_"
  raster_troposphere_no2_address <- 
    paste0(raster_troposphere_no2_folder, aim_year, "m", aim_month, 
           "_WGS84_.25.25_monthly_tropospheric_no2.tif")
  raster_troposphere_no2 <- raster::raster(raster_troposphere_no2_address)
  raster_troposphere_no2 <- flip(raster_troposphere_no2, direction = 'y')
  mol_g = 6.022140857 * 10^23  # mol
  raster_troposphere_no2 <- raster_troposphere_no2 / mol_g * 46.0055 # convert mol to g
  raster_troposphere_no2 <- raster_troposphere_no2 * 10000 * 1000000 # conver /cm2 to /m2 and g to ug
  predict_part_troposphere_no2 <- (raster_troposphere_no2 - MEAN_ug_m2_troposphere_no2.kriged.raster) *
    CO_ug_m2_troposphere_no2.kriged.raster
  
  raster_ter_pressure_folder <- 
    "D:/10_Article/09_TempOutput/03_MonthlyTerrainPressureTif/OMI-Aura_L2G-OMNO2G_"
  raster_ter_pressure_address <- 
    paste0(raster_ter_pressure_folder, aim_year, "m", aim_month, 
           "_WGS84_.25.25_monthly_terrain_pressure.tif")
  raster_ter_pressure <- raster::raster(raster_ter_pressure_address)
  raster_ter_pressure <- flip(raster_ter_pressure, direction = 'y')
  predict_part_ter_pressure <- (raster_ter_pressure - MEAN_ter_pressure.kriged.raster) *
    CO_ter_pressure.kriged.raster
  
  dateUse <- as.Date(paste0(aim_year,"-",aim_month,"-01"))
  day_number <- strftime(dateUse, format = "%j") 
  raster_day_temp_folder <- 
    "D:/10_Article/09_TempOutput/04_MonthlyTemperatureTif/Surf_Temp_Monthly_005dg_v6/LST_Day_CMG/MOD11C3_LST_Day_CMG_"
  raster_day_temp_address <- 
    paste0(raster_day_temp_folder, aim_year, "_", day_number, ".tif")
  raster_day_temp <- raster::raster(raster_day_temp_address)
  raster_night_temp_folder <- 
    "D:/10_Article/09_TempOutput/04_MonthlyTemperatureTif/Surf_Temp_Monthly_005dg_v6/LST_Night_CMG/MOD11C3_LST_Night_CMG_"
  raster_night_temp_address <- 
    paste0(raster_night_temp_folder, aim_year, "_", day_number, ".tif")
  raster_night_temp <- raster::raster(raster_night_temp_address)
  raster_temp <- ((raster_day_temp * 0.02 - 273.16) + (raster_night_temp * 0.02 - 273.16))
  predict_part_raster_temp <- (raster_temp - MEAN_temp.kriged.raster) *
    CO_temp.kriged.raster
  
  raster_ndvi_folder <- 
    "D:/10_Article/09_TempOutput/05_MonthlyNDVITif/VI_Monthly_005dg_v6/NDVI/MOD13C2_NDVI_"
  raster_ndvi_address <- 
    paste0(raster_ndvi_folder, aim_year, "_", day_number, ".tif")
  raster_ndvi <- raster::raster(raster_ndvi_address)
  raster_ndvi <- raster_ndvi / 10000
  predict_part_ndvi <- (raster_ndvi - MEAN_ndvi.kriged.raster) *
    CO_ndvi.kriged.raster
  
  raster_precipitation_folder <- 
    "D:/10_Article/09_TempOutput/07_MonthlyPrecipitationTif/Add025Outline/add_totalPrecipitationRate"
  raster_precipitation_address <- 
    paste0(raster_precipitation_folder, aim_year, aim_month, ".tif")
  raster_precipitation <- raster::raster(raster_precipitation_address)
  raster_precipitation <- flip(raster_precipitation, direction = 'y')
  raster_precipitation <- raster_precipitation / 3600 # now, the precipitation unit is kg/(m2 * h)  
  predict_part_precipitation <- (raster_precipitation - MEAN_precipitation.kriged.raster) *
    CO_precipitation.kriged.raster
  
  raster_PBLH_folder <-
    "D:/10_Article/09_TempOutput/10_PlanetaryBoundaryLayerHeight/Resample/"
  raster_PBLH_address <- 
    paste0(raster_PBLH_folder, aim_year, aim_month, ".2525.tif")
  raster_PBLH <- raster::raster(raster_PBLH_address)
  predict_part_PBLH <- (raster_PBLH - MEAN_PBLH.kriged.raster) *
    CO_PBLH.kriged.raster
  
  predict.raster <- predict_part_troposphere_no2 + predict_part_ter_pressure + predict_part_precipitation +
    predict_part_raster_temp + predict_part_ndvi + predict_part_PBLH
  
  if (aim_year == "2016")
  {
    predict.raster <- predict.raster + CO_Y2016.kriged.raster * 
      (1 - MEAN_Y2016.kriged.raster)
  }
  if (aim_year == "2017")
  {
    predict.raster <- predict.raster + CO_Y2017.kriged.raster * 
      (1 - MEAN_Y2017.kriged.raster)
  }
  if (aim_year == "2018")
  {
    predict.raster <- predict.raster + CO_Y2018.kriged.raster * 
      (1 - MEAN_Y2018.kriged.raster)
  }
  if (aim_year == "2019")
  {
    predict.raster <- predict.raster + CO_Y2019.kriged.raster * 
      (1 - MEAN_Y2019.kriged.raster)
  }
  if (aim_year == "2020")
  {
    predict.raster <- predict.raster + CO_Y2020.kriged.raster * 
      (1 - MEAN_Y2020.kriged.raster)
  }
  if (aim_year == "2021")
  {
    predict.raster <- predict.raster + CO_Y2021.kriged.raster * 
      (1 - MEAN_Y2021.kriged.raster)
  }
  predict.raster <- predict.raster + MEAN_no2_measured_ug.m3.kriged.raster
  values(predict.raster)[values(predict.raster) < 0] = 0
  return(predict.raster)
}

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

month.list <- c("01", "02", "03", "04", "05", "06",
                "07", "08", "09", "10", "11", "12")
year.list <- c("2015", "2016", "2017", "2018", "2019", "2020", "2021")
predict_raster_folder <- "D:/10_Article/11_PredictRaster/01_Test0104/"
predict_jpg_folder <- "D:/10_Article/11_PredictRaster/02_Test0104JPG/"
brks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 260)
pal <- colorRampPalette(c("blue","green","yellow","red"))

for (aim_year in year.list){
  for (aim_month in month.list){
    predict.raster <- makePredictRaster(aim_year, aim_month = aim_month)
    writeRaster(predict.raster, filename = paste0(predict_raster_folder,aim_year,aim_month,".tif"),
                format="GTiff", overwrite=TRUE)
    jpeg(paste0(predict_jpg_folder,aim_year,aim_month,".jpg"), 
         quality = 300, width = 1800, height = 1000)
    plot(predict.raster, breaks = brks, col = pal(12))
    title(paste0(aim_year, "-", aim_month))
    dev.off()
  }
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

# get the point prediction
predict_raster_folder <- "D:/10_Article/11_PredictRaster/01_Test0104/"
filelist <- list.files(predict_raster_folder)
predictGroundNo2 <- 
  extractPointDataFromRaster(predict_raster_folder, filelist, cityLocationSpatialPoint,
                             1, 5, F, "predict_no2")
load("03_Rawdata/usedDataset.Rdata")
usedDataset$month <- usedDataset$month %>% as.character() %>% as.numeric()
usedDataset$year <- usedDataset$year %>% as.character() %>% as.numeric()
testDataset <- left_join(usedDataset, predictGroundNo2, by = c("City", "Country", "year", "month"))
testDataset <- testDataset %>%
  filter(!is.na(no2_measured_ug.m3)) %>%
  filter(!is.na(predict_no2))
testDataset$residuals <- testDataset$predict_no2 - testDataset$no2_measured_ug.m3
testDataset <- testDataset %>%
  dplyr::select(Country, City, year, month, CityCode, no2_measured_ug.m3,
                predict_no2, residuals)
load("03_Rawdata/usedDataset.Rdata")
inFlag <- usedDataset %>% dplyr::select(CityCode) %>% distinct()
inFlag$flag <- 1
testDataset <- left_join(testDataset, inFlag, by = "CityCode")
testDataset <- testDataset %>%
  filter(flag == 1)

cor.test(testDataset$no2_measured_ug.m3, abs(testDataset$residuals))
cor.test(testDataset$no2_measured_ug.m3, testDataset$predict_no2)
lm(no2_measured_ug.m3 ~ predict_no2, testDataset) %>% summary()
r2 <- 1 - sum( (testDataset$no2_measured_ug.m3 - testDataset$predict_no2)^2 ) /
  sum((testDataset$no2_measured_ug.m3 - mean(testDataset$no2_measured_ug.m3))^2)
r2
(sum( (testDataset$no2_measured_ug.m3 - testDataset$predict_no2)^2 ) / nrow(testDataset)) %>% sqrt()
mean( abs( (testDataset$no2_measured_ug.m3 - testDataset$predict_no2) ) )
plot(testDataset$no2_measured_ug.m3, testDataset$predict_no2)
save(testDataset, file = "04_Results/FinalRasterCrossValidation.Rdata")

#jpg.list <- list.files(predict_jpg_folder)
#frames <- paste0(predict_jpg_folder, jpg.list)
#m <- image_read(frames)
#m <- image_animate(m, fps = 2)
#image_write(m, 
#            paste0("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/06_Animate/",
#                   "ani.gif"))
