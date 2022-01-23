# Author: M.L.

# input: idw_COEF_raster.RData
# idw_COEF_raster.RData: "CO_ug_m2_troposphere_no2.idw.raster" interpolation of troposheric no2 coefficient (IDW)
# idw_COEF_raster.RData: "CO_ndvi.idw.raster" interpolation of ndvi coefficient based on ordinary kriging (IDW)
# idw_COEF_raster.RData: "CO_temp.idw.raster" interpolation of night temperature coefficient (IDW)
# idw_COEF_raster.RData: "CO_PBLH.idw.raster" interpolation of PBLH coefficient (IDW)
# idw_COEF_raster.RData: "CO_precipitation.idw.raster" interpolation of precipitation coefficient (IDW)
# idw_COEF_raster.RData: "CO_ter_pressure.idw.raster" interpolation of air pressure coefficient (IDW)
# idw_COEF_raster.RData: "CO_Y2016.idw.raster" interpolation of 2016 dummy variable coefficient (IDW)
# idw_COEF_raster.RData: "CO_Y2017.idw.raster" interpolation of 2017 dummy variable coefficient (IDW)
# idw_COEF_raster.RData: "CO_Y2018.idw.raster" interpolation of 2018 dummy variable coefficient (IDW)
# idw_COEF_raster.RData: "CO_Y2019.idw.raster" interpolation of 2019 dummy variable coefficient (IDW)
# idw_COEF_raster.RData: "CO_Y2020.idw.raster" interpolation of 2020 dummy variable coefficient (IDW)
# idw_COEF_raster.RData: "CO_Y2021.idw.raster" interpolation of 2021 dummy variable coefficient (IDW)

# output: idw.MEAN_raster.RData
# idw.MEAN_raster.RData: "MEAN_ug_m2_troposphere_no2.idw.raster" interpolation of troposheric no2 mean value
#                                                               in the usedDataset based on ordinary kriging (IDW)
# idw.MEAN_raster.RData: "MEAN_ndvi.idw.raster" interpolation of ndvi mean value (IDW)
# idw.MEAN_raster.RData: "MEAN_temp.idw.raster" interpolation of temperature mean value (IDW)
# idw.MEAN_raster.RData: "MEAN_NTL.idw.raster" interpolation of night time light mean value (IDW)
# idw.MEAN_raster.RData: "MEAN_PBLH.idw.raster" interpolation of PBLH mean value (IDW)
# idw.MEAN_raster.RData: "MEAN_precipitation.idw.raster" interpolation of precipitation mean value (IDW)
# idw.MEAN_raster.RData: "MEAN_ter_pressure.idw.raster" interpolation of air pressure mean value (IDW)

# end

library(tidyverse)
library(dplyr)
library(raster)
library(lubridate)

makePredictRaster <- function(aim_year, aim_month){
  load("05_CoefficientRaster/idw_COEF_raster.RData")
  load("05_CoefficientRaster/idw.MEAN_raster.RData")
  
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
  predict_part_troposphere_no2 <- (raster_troposphere_no2 - MEAN_ug_m2_troposphere_no2.idw.raster) *
    CO_ug_m2_troposphere_no2.idw.raster
  
  raster_ter_pressure_folder <- 
    "D:/10_Article/09_TempOutput/03_MonthlyTerrainPressureTif/OMI-Aura_L2G-OMNO2G_"
  raster_ter_pressure_address <- 
    paste0(raster_ter_pressure_folder, aim_year, "m", aim_month, 
           "_WGS84_.25.25_monthly_terrain_pressure.tif")
  raster_ter_pressure <- raster::raster(raster_ter_pressure_address)
  raster_ter_pressure <- flip(raster_ter_pressure, direction = 'y')
  predict_part_ter_pressure <- (raster_ter_pressure - MEAN_ter_pressure.idw.raster) *
    CO_ter_pressure.idw.raster
  
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
  raster_temp <- ((raster_day_temp * 0.02 - 273.16) + (raster_night_temp * 0.02 - 273.16))/2
  predict_part_raster_temp <- (raster_temp - MEAN_temp.idw.raster) *
    CO_temp.idw.raster
  
  raster_ndvi_folder <- 
    "D:/10_Article/09_TempOutput/05_MonthlyNDVITif/VI_Monthly_005dg_v6/NDVI/MOD13C2_NDVI_"
  raster_ndvi_address <- 
    paste0(raster_ndvi_folder, aim_year, "_", day_number, ".tif")
  raster_ndvi <- raster::raster(raster_ndvi_address)
  raster_ndvi <- raster_ndvi / 10000
  predict_part_ndvi <- (raster_ndvi - MEAN_ndvi.idw.raster) *
    CO_ndvi.idw.raster
  
  raster_precipitation_folder <- 
    "D:/10_Article/09_TempOutput/07_MonthlyPrecipitationTif/Add025Outline/add_totalPrecipitationRate"
  raster_precipitation_address <- 
    paste0(raster_precipitation_folder, aim_year, aim_month, ".tif")
  raster_precipitation <- raster::raster(raster_precipitation_address)
  raster_precipitation <- flip(raster_precipitation, direction = 'y')
  raster_precipitation <- raster_precipitation / 3600 # now, the precipitation unit is kg/(m2 * h)  
  predict_part_precipitation <- (raster_precipitation - MEAN_precipitation.idw.raster) *
    CO_precipitation.idw.raster
  
  raster_PBLH_folder <-
    "D:/10_Article/09_TempOutput/10_PlanetaryBoundaryLayerHeight/Resample/"
  raster_PBLH_address <- 
    paste0(raster_PBLH_folder, aim_year, aim_month, ".2525.tif")
  raster_PBLH <- raster::raster(raster_PBLH_address)
  predict_part_PBLH <- (raster_PBLH - MEAN_PBLH.idw.raster) *
    CO_PBLH.idw.raster
  
  predict.raster <- predict_part_troposphere_no2 + predict_part_ter_pressure + predict_part_precipitation +
    predict_part_raster_temp + predict_part_ndvi + predict_part_PBLH
  
  Y2016_dummy = 0
  Y2017_dummy = 0
  Y2018_dummy = 0
  Y2019_dummy = 0
  Y2020_dummy = 0
  Y2021_dummy = 0
  
  if(aim_year == "2016") {Y2016_dummy = 1}
  if(aim_year == "2017") {Y2017_dummy = 1}
  if(aim_year == "2018") {Y2018_dummy = 1}
  if(aim_year == "2019") {Y2019_dummy = 1}
  if(aim_year == "2020") {Y2020_dummy = 1}
  if(aim_year == "2021") {Y2021_dummy = 1}
  
  predict_part_Y2016 <-  CO_Y2016.idw.raster * (Y2016_dummy - MEAN_Y2016.idw.raster)
  predict_part_Y2017 <-  CO_Y2017.idw.raster * (Y2017_dummy - MEAN_Y2017.idw.raster)
  predict_part_Y2018 <-  CO_Y2018.idw.raster * (Y2018_dummy - MEAN_Y2018.idw.raster)
  predict_part_Y2019 <-  CO_Y2019.idw.raster * (Y2019_dummy - MEAN_Y2019.idw.raster)
  predict_part_Y2020 <-  CO_Y2020.idw.raster * (Y2020_dummy - MEAN_Y2020.idw.raster)
  predict_part_Y2021 <-  CO_Y2021.idw.raster * (Y2021_dummy - MEAN_Y2021.idw.raster)
  
  predict.raster <- predict.raster + 
    predict_part_Y2016 + predict_part_Y2017 + predict_part_Y2018 +
    predict_part_Y2019 + predict_part_Y2020 + predict_part_Y2021
  
  predict.raster <- predict.raster + MEAN_no2_measured_ug.m3.idw.raster
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
predict_raster_folder <- "D:/10_Article/11_PredictRaster/03_Test0122/"
predict_jpg_folder <- "D:/10_Article/11_PredictRaster/04_Test0122JPG/"
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
predict_raster_folder <- "D:/10_Article/11_PredictRaster/03_Test0122/"
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
save(testDataset, file = "04_Results/idw.FinalRasterCrossValidation.Rdata")

#jpg.list <- list.files(predict_jpg_folder)
#frames <- paste0(predict_jpg_folder, jpg.list)
#m <- image_read(frames)
#m <- image_animate(m, fps = 2)
#image_write(m, 
#            paste0("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/06_Animate/",
#                   "ani.gif"))
