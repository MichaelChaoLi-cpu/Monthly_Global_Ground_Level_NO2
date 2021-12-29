# end

library(tidyverse)
library(gstat)
library(sp) 
library(raster)
library(dplyr)
library(tmap)
library(automap)

# make the raster -180 -60 180 90
nx = 600                                       # number of cells in the x direction
ny = 1440                                     # number of cells in the y direction
xmin = -179.875                                     # x coordinate of lower, left cell center 
ymin = -59.875                                     # y coordinate of lower, left cell center 
xsize = 0.25                                   # extent of cells in x direction
ysize = 0.25                                   # extent of cells in y direction
proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

addcoord <- function(nx,xmin,xsize,ny,ymin,ysize) { # Michael Pyrcz, March, 2018                      
  # makes a 2D dataframe with coordinates based on GSLIB specification
  coords = matrix(nrow = nx*ny,ncol=2)
  ixy = 1
  for(iy in 1:nx) {
    for(ix in 1:ny) {
      coords[ixy,1] = xmin + (ix-1)*xsize  
      coords[ixy,2] = ymin + (iy-1)*ysize 
      ixy = ixy + 1
    }
  }
  coords.df = data.frame(coords)
  colnames(coords.df) <- c("X","Y")
  coordinates(coords.df) =~X+Y
  return (coords.df)
  
}
coords <- addcoord(nx,xmin,xsize,ny,ymin,ysize)
coords@proj4string <- CRS(proj)

### make the mean raster
load("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/03_Rawdata/mergedDataset.Rdata")

na.test <- mergedDataset %>% na.omit()
na.test$count <- 1
na.test <- aggregate(na.test$count, by = list(na.test$City, na.test$Country), FUN=sum)
colnames(na.test) <- c("City", "Country", "RecordCount")
na.test <- na.test %>% filter(RecordCount > 17) #freedom is 8

usedDataset <- left_join(mergedDataset, na.test, by = c("City", "Country"))
usedDataset <- usedDataset %>% filter(!is.na(RecordCount))
usedDataset <- usedDataset %>% dplyr::select(
  no2_measured_mg.m3,
  no2, mg_m2_total_no2, mg_m2_troposphere_no2,
  #mg_m2_total_no2_lag, mg_m2_troposphere_no2_lag,
  ter_pressure, dayTimeTemperature, nightTimeTemperature, ndvi,
  humidity, precipitation, NTL, speedwind, PBLH, 
  #UVAerosolIndex, ozone, cloudfraction, cloudpressure,
  CityCode, City, Country, 
  month, year, Date, Y2016, Y2017, Y2018, Y2019, Y2020, Y2021
) %>% na.omit()
usedDataset$humidity <- usedDataset$humidity %>% as.numeric()
usedDataset$precipitation <- usedDataset$precipitation %>% as.numeric()
usedDataset$year <- usedDataset$year %>% as.character() %>% as.numeric()
usedDataset$month <- usedDataset$month %>% as.character() %>% as.numeric()
usedDataset$period <- usedDataset$year * 100 + usedDataset$month
##### this is a strange value
usedDataset <- usedDataset %>% filter(CityCode != 499)
# preprocessing of the panel data set not we take the total column as the dependent variable

formula.CV.FEM <-
  no2_measured_mg.m3 ~ mg_m2_troposphere_no2 + ter_pressure + dayTimeTemperature + nightTimeTemperature +
  ndvi + humidity + precipitation + NTL + speedwind + PBLH + 
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021
rawCrossValidationDataset <- usedDataset %>% 
  dplyr::select("CityCode", "period", all.vars(formula.CV.FEM))
meanValueOfVariables <- stats::aggregate(rawCrossValidationDataset[,all.vars(formula.CV.FEM)],
                                         by = list(rawCrossValidationDataset$CityCode), mean)
colnames(meanValueOfVariables)[1] <- "CityCode"

cityLocation <- read.csv("D:/10_Article/01_RawData/12_LocationJson/CityLocationOfficial.csv",
                         encoding="UTF-8") %>%
  dplyr::select(X0, X1, X2, X3)
colnames(cityLocation) <- c("Latitude", "Longitude", "City", "Country")
cityLocation <- cityLocation %>% 
  mutate(City = ifelse(City == "Washington, D.C.", "Washington D.C.", City))
cityNameCode <- usedDataset %>% dplyr::select(CityCode, City, Country) %>% distinct()
cityLocation <- left_join(cityNameCode, cityLocation, by = c("City", "Country"))
rm(cityNameCode)

proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
xy <- cityLocation %>% dplyr::select(Longitude, Latitude)
cityLocationSpatialPoint <- SpatialPointsDataFrame(coords = xy, data = cityLocation[,c(1, 2, 3, 4, 5)],
                                                   proj4string = CRS(proj))
rm(xy)
# get the city points 

meanValueSpatialDataFrame <- cityLocationSpatialPoint
meanValueSpatialDataFrame@data <- left_join(meanValueSpatialDataFrame@data, meanValueOfVariables,
                                            by = "CityCode")

### no2_measured_mg.m3
no2_measured_mg.m3_emp_OK <- gstat::variogram(
  no2_measured_mg.m3 ~ 1, meanValueSpatialDataFrame
)
plot(no2_measured_mg.m3_emp_OK)
dat.fit  <- fit.variogram(no2_measured_mg.m3_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
no2_measured_mg.m3.kriged <- 
  krige(no2_measured_mg.m3~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
no2_measured_mg.m3.kriged@data$Ttest <- no2_measured_mg.m3.kriged@data$var1.pred / 
  no2_measured_mg.m3.kriged@data$var1.var
no2_measured_mg.m3.kriged.raster <- as(no2_measured_mg.m3.kriged, 'SpatialPixelsDataFrame')
no2_measured_mg.m3.kriged.raster <- as(no2_measured_mg.m3.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, no2_measured_mg.m3.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$no2_measured_mg.m3)^2) / 
  sum( ( meanValueSpatialDataFrame$no2_measured_mg.m3 - mean(meanValueSpatialDataFrame$no2_measured_mg.m3) )^2 )
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$no2_measured_mg.m3)
no2_measured_mg.m3.kriged.raster@data <- no2_measured_mg.m3.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_no2_measured_mg.m3.kriged.raster <- raster(no2_measured_mg.m3.kriged.raster)
rm(no2_measured_mg.m3_emp_OK, dat.fit, no2_measured_mg.m3.kriged, 
   no2_measured_mg.m3.kriged.raster)


### mg_m2_troposphere_no2
mg_m2_troposphere_no2_emp_OK <- gstat::variogram(
  mg_m2_troposphere_no2 ~ 1, meanValueSpatialDataFrame
)
plot(mg_m2_troposphere_no2_emp_OK)
dat.fit  <- fit.variogram(mg_m2_troposphere_no2_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
mg_m2_troposphere_no2.kriged <- 
  krige(mg_m2_troposphere_no2~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
mg_m2_troposphere_no2.kriged@data$Ttest <- mg_m2_troposphere_no2.kriged@data$var1.pred / 
  mg_m2_troposphere_no2.kriged@data$var1.var
mg_m2_troposphere_no2.kriged.raster <- as(mg_m2_troposphere_no2.kriged, 'SpatialPixelsDataFrame')
mg_m2_troposphere_no2.kriged.raster <- as(mg_m2_troposphere_no2.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, mg_m2_troposphere_no2.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$mg_m2_troposphere_no2)^2) / 
  sum( ( meanValueSpatialDataFrame$mg_m2_troposphere_no2 - mean(meanValueSpatialDataFrame$mg_m2_troposphere_no2) )^2 )
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$mg_m2_troposphere_no2)
mg_m2_troposphere_no2.kriged.raster@data <- mg_m2_troposphere_no2.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_mg_m2_troposphere_no2.kriged.raster <- raster(mg_m2_troposphere_no2.kriged.raster)
rm(mg_m2_troposphere_no2_emp_OK, dat.fit, mg_m2_troposphere_no2.kriged, 
   mg_m2_troposphere_no2.kriged.raster)

### ter_pressure
ter_pressure_emp_OK <- gstat::variogram(
  ter_pressure ~ 1, meanValueSpatialDataFrame
)
plot(ter_pressure_emp_OK)
dat.fit  <- fit.variogram(ter_pressure_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
ter_pressure.kriged <- 
  krige(ter_pressure~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
ter_pressure.kriged@data$Ttest <- ter_pressure.kriged@data$var1.pred / 
  ter_pressure.kriged@data$var1.var
ter_pressure.kriged.raster <- as(ter_pressure.kriged, 'SpatialPixelsDataFrame')
ter_pressure.kriged.raster <- as(ter_pressure.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, ter_pressure.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$ter_pressure)^2) / 
  sum( ( meanValueSpatialDataFrame$ter_pressure - mean(meanValueSpatialDataFrame$ter_pressure) )^2 )
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$ter_pressure)
ter_pressure.kriged.raster@data <- ter_pressure.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_ter_pressure.kriged.raster <- raster(ter_pressure.kriged.raster)
rm(ter_pressure_emp_OK, dat.fit, ter_pressure.kriged, 
   ter_pressure.kriged.raster)

### dayTimeTemperature
dayTimeTemperature_emp_OK <- gstat::variogram(
  dayTimeTemperature ~ 1, meanValueSpatialDataFrame
)
plot(dayTimeTemperature_emp_OK)
dat.fit  <- fit.variogram(dayTimeTemperature_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
dayTimeTemperature.kriged <- 
  krige(dayTimeTemperature~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
dayTimeTemperature.kriged@data$Ttest <- dayTimeTemperature.kriged@data$var1.pred / 
  dayTimeTemperature.kriged@data$var1.var
dayTimeTemperature.kriged.raster <- as(dayTimeTemperature.kriged, 'SpatialPixelsDataFrame')
dayTimeTemperature.kriged.raster <- as(dayTimeTemperature.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, dayTimeTemperature.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$dayTimeTemperature)^2) / 
  sum( ( meanValueSpatialDataFrame$dayTimeTemperature - mean(meanValueSpatialDataFrame$dayTimeTemperature) )^2 )
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$dayTimeTemperature)
dayTimeTemperature.kriged.raster@data <- dayTimeTemperature.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_dayTimeTemperature.kriged.raster <- raster(dayTimeTemperature.kriged.raster)
rm(dayTimeTemperature_emp_OK, dat.fit, dayTimeTemperature.kriged, 
   dayTimeTemperature.kriged.raster)

### nightTimeTemperature
nightTimeTemperature_emp_OK <- gstat::variogram(
  nightTimeTemperature ~ 1, meanValueSpatialDataFrame
)
plot(nightTimeTemperature_emp_OK)
dat.fit  <- fit.variogram(nightTimeTemperature_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
nightTimeTemperature.kriged <- 
  krige(nightTimeTemperature~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
nightTimeTemperature.kriged@data$Ttest <- nightTimeTemperature.kriged@data$var1.pred / 
  nightTimeTemperature.kriged@data$var1.var
nightTimeTemperature.kriged.raster <- as(nightTimeTemperature.kriged, 'SpatialPixelsDataFrame')
nightTimeTemperature.kriged.raster <- as(nightTimeTemperature.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, nightTimeTemperature.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$nightTimeTemperature)^2) / 
  sum( ( meanValueSpatialDataFrame$nightTimeTemperature - mean(meanValueSpatialDataFrame$nightTimeTemperature) )^2 )
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$nightTimeTemperature)
nightTimeTemperature.kriged.raster@data <- nightTimeTemperature.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_nightTimeTemperature.kriged.raster <- raster(nightTimeTemperature.kriged.raster)
rm(nightTimeTemperature_emp_OK, dat.fit, nightTimeTemperature.kriged, 
   nightTimeTemperature.kriged.raster)

### ndvi
ndvi_emp_OK <- gstat::variogram(
  ndvi ~ 1, meanValueSpatialDataFrame
)
plot(ndvi_emp_OK)
dat.fit  <- fit.variogram(ndvi_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
ndvi.kriged <- 
  krige(ndvi~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
ndvi.kriged@data$Ttest <- ndvi.kriged@data$var1.pred / 
  ndvi.kriged@data$var1.var
ndvi.kriged.raster <- as(ndvi.kriged, 'SpatialPixelsDataFrame')
ndvi.kriged.raster <- as(ndvi.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, ndvi.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$ndvi)^2) / 
  sum( ( meanValueSpatialDataFrame$ndvi - mean(meanValueSpatialDataFrame$ndvi) )^2 )
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$ndvi)
ndvi.kriged.raster@data <- ndvi.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_ndvi.kriged.raster <- raster(ndvi.kriged.raster)
rm(ndvi_emp_OK, dat.fit, ndvi.kriged, 
   ndvi.kriged.raster)
