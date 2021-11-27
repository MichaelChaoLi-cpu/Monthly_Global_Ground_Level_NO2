

library(tidyverse)
library(dplyr)
library(raster)
library(sp)
library(stringr)
library(rgdal)

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
totalNo2RasterDataset$g_m2_total_no2 <- totalNo2RasterDataset$g_cm2 * 10000 # conver /cm2 to /m2
totalNo2RasterDataset <- totalNo2RasterDataset %>% dplyr::select("Country", "City", "year", "month", "g_m2_total_no2")



#get monthly troposphere no2 from the OMNO2G, band 9(troposphere no2)
troposphereNo2RasterFolder <- "D:/10_Article/09_TempOutput/02_MonthlyTroposphericNo2Tif/"
filelist <- list.files(troposphereNo2RasterFolder)
troposphereNo2RasterDataset <- 
  extractPointDataFromRaster(troposphereNo2RasterFolder, filelist, cityLocationSpatialPoint,
                             21, 26, T, "raw_no2")
# convert molecular / cm2 to ug / m2
mol_g = 6.022140857 * 10^23  # mol
troposphereNo2RasterDataset$g_cm2 <- troposphereNo2RasterDataset$raw_no2 / mol_g * 46.0055 # convert mol to g
troposphereNo2RasterDataset$g_m2_troposphere_no2 <- troposphereNo2RasterDataset$g_cm2 * 10000 # conver /cm2 to /m2
troposphereNo2RasterDataset <- troposphereNo2RasterDataset %>% dplyr::select("Country", "City", "year", "month", "g_m2_troposphere_no2")



#get monthly terrain pressure from the OMNO2G, band 29 (terrain pressure)
terrainPressureRasterFolder <- "D:/10_Article/09_TempOutput/03_MonthlyTerrainPressureTif/"
filelist <- list.files(terrainPressureRasterFolder)
terrainPressureRasterDataset <- 
  extractPointDataFromRaster(terrainPressureRasterFolder, filelist, cityLocationSpatialPoint,
                             21, 26, T, "ter_pressure")

#get monthly Daytime temperature from the MOD11C3
dayTimeTemperatureRasterFolder <- "D:/10_Article/09_TempOutput/04_MonthlyTemperatureTif/Surf_Temp_Monthly_005dg_v6/LST_Day_CMG/"
filelist <- list.files(dayTimeTemperatureRasterFolder)
dayTimeTemperatureRasterDataset <- 
  extractPointDataFromRaster(dayTimeTemperatureRasterFolder, filelist, cityLocationSpatialPoint,
                             21, month_start_location = 26, F,
                             "dayTimeTemperature", month_end_location = 28)


#break point
test <- left_join(totalNo2Dataset, totalNo2RasterDataset, by = c("Country", "City", "year", "month"))
test <- left_join(test, troposphereNo2RasterDataset, by = c("Country", "City", "year", "month"))
cor.test(test$g_m2_total_no2, test$no2)
cor.test(test$g_m2_troposphere_no2, test$no2)
test <- left_join(test, terrainPressureRasterDataset, by = c("Country", "City", "year", "month"))
cor.test(test$ter_pressure, test$no2)

test$period <- test$year * 100 + test$month

library(plm)

pdata <- pdata.frame(test, index = c("CityCode", "period"))
formula <- no2 ~ g_m2_total_no2 + ter_pressure
ols <- plm(formula, pdata, model = "pooling")
summary(ols)
fem <- plm(formula, pdata, model = "within")
summary(fem)
rem <- plm(formula, pdata, model = "random")
summary(rem)
