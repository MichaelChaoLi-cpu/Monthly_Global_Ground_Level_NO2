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
                                dplyr::select(GridID, raw)
                              monthly_data <- monthly_data %>%
                                dplyr::mutate(year = Year,
                                              month = Month)
                              monthly_data <- monthly_data %>% na.omit()
                            }
  stopCluster(cl)
  colnames(RasterDataset) <- c("GridID", aimed_column_name, "year", "month")
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


