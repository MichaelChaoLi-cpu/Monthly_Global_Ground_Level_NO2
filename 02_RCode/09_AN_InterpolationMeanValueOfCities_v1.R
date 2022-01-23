# Author: M.L.

# input: usedDataset.Rdata
# usedDataset.Rdata "no2" monthly average no2 concentration ppm.
# usedDataset.Rdata "mg_m2_total_no2" monthly average total column amount of no2 (mg/m2)
# usedDataset.Rdata "mg_m2_troposphere_no2" monthly average tropospheric column amount of no2 (mg/m2)
# usedDataset.Rdata "ter_pressure" monthly average terrain surface pressure (hpa)
# usedDataset.Rdata "temp" monthly average daily temperature (C)
# usedDataset.Rdata "ndvi" NDVI -1 to 1
# usedDataset.Rdata "precipitation" the precipitation unit is kg/(m2 * h)
# usedDataset.Rdata "NTL" nighttime light
# usedDataset.Rdata "PBLH" planetary boundary layer height unit is m.
# usedDataset.Rdata "ug_m2_total_no2" monthly average total column amount of no2 (ug/m2)
# usedDataset.Rdata "ug_m2_troposphere_no2" monthly average tropospheric column amount of no2 (ug/m2)
# usedDataset.Rdata "no2_measured_ug.m3" (ug/m3)

# output: MEAN_raster.RData
# MEAN_raster.RData: "MEAN_ug_m2_troposphere_no2.kriged.raster" interpolation of troposheric no2 mean value
#                                                               in the usedDataset based on ordinary kriging (OK)
# MEAN_raster.RData: "MEAN_ndvi.kriged.raster" interpolation of ndvi mean value (OK)
# MEAN_raster.RData: "MEAN_temp.kriged.raster" interpolation of temperature mean value (OK)
# MEAN_raster.RData: "MEAN_NTL.kriged.raster" interpolation of night time light mean value (OK)
# MEAN_raster.RData: "MEAN_PBLH.kriged.raster" interpolation of PBLH mean value (OK)
# MEAN_raster.RData: "MEAN_precipitation.kriged.raster" interpolation of precipitation mean value (OK)
# MEAN_raster.RData: "MEAN_ter_pressure.kriged.raster" interpolation of air pressure mean value (OK)

# output: krigingMeanResult.RData
# krigingMeanResult.RData: "N" number of observation
# krigingMeanResult.RData: "R2" R2
# krigingMeanResult.RData: "RMSE" root mean square error
# krigingMeanResult.RData: "MAE" mean absolute error
# krigingMeanResult.RData: "r" correlation coefficient
# krigingMeanResult.RData: "Intercept" intercept of regression
# krigingMeanResult.RData: "Slope" slope of regression
# krigingMeanResult.RData: "onPoint.R2" R2 after interpolation

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

load("03_Rawdata/usedDataset.Rdata")

formula.CV.FEM <-
  no2_measured_ug.m3 ~ ug_m2_troposphere_no2 + ter_pressure + temp +
  ndvi +  precipitation + PBLH + 
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

kriging.cv.mean.dataset <- data.frame(Doubles = double(), Ints = integer(),
                                        Factors = factor(), Logicals = logical(),
                                        Characters = character(), stringsAsFactors = FALSE)

### no2_measured_ug.m3
no2_measured_ug.m3_emp_OK <- gstat::variogram(
  no2_measured_ug.m3 ~ 1, meanValueSpatialDataFrame
)
plot(no2_measured_ug.m3_emp_OK)
dat.fit  <- fit.variogram(no2_measured_ug.m3_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(no2_measured_ug.m3~1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

no2_measured_ug.m3.kriged <- 
  krige(no2_measured_ug.m3~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
no2_measured_ug.m3.kriged@data$Ttest <- no2_measured_ug.m3.kriged@data$var1.pred / 
  no2_measured_ug.m3.kriged@data$var1.var
no2_measured_ug.m3.kriged.raster <- as(no2_measured_ug.m3.kriged, 'SpatialPixelsDataFrame')
no2_measured_ug.m3.kriged.raster <- as(no2_measured_ug.m3.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, no2_measured_ug.m3.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$no2_measured_ug.m3)^2) / 
  sum( ( meanValueSpatialDataFrame$no2_measured_ug.m3 - mean(meanValueSpatialDataFrame$no2_measured_ug.m3) )^2 )
cv.line <- c("no2_measured_ug.m3.mean", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$no2_measured_ug.m3)
no2_measured_ug.m3.kriged.raster@data <- no2_measured_ug.m3.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_no2_measured_ug.m3.kriged.raster <- raster(no2_measured_ug.m3.kriged.raster)
rm(no2_measured_ug.m3_emp_OK, dat.fit, no2_measured_ug.m3.kriged, 
   no2_measured_ug.m3.kriged.raster)


### ug_m2_troposphere_no2
ug_m2_troposphere_no2_emp_OK <- gstat::variogram(
  ug_m2_troposphere_no2 ~ 1, meanValueSpatialDataFrame
)
plot(ug_m2_troposphere_no2_emp_OK)
dat.fit  <- fit.variogram(ug_m2_troposphere_no2_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(ug_m2_troposphere_no2~1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

ug_m2_troposphere_no2.kriged <- 
  krige(ug_m2_troposphere_no2~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
ug_m2_troposphere_no2.kriged@data$Ttest <- ug_m2_troposphere_no2.kriged@data$var1.pred / 
  ug_m2_troposphere_no2.kriged@data$var1.var
ug_m2_troposphere_no2.kriged.raster <- as(ug_m2_troposphere_no2.kriged, 'SpatialPixelsDataFrame')
ug_m2_troposphere_no2.kriged.raster <- as(ug_m2_troposphere_no2.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, ug_m2_troposphere_no2.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$ug_m2_troposphere_no2)^2) / 
  sum( ( meanValueSpatialDataFrame$ug_m2_troposphere_no2 - mean(meanValueSpatialDataFrame$ug_m2_troposphere_no2) )^2 )
cv.line <- c("ug_m2_troposphere_no2.mean",
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$ug_m2_troposphere_no2)
ug_m2_troposphere_no2.kriged.raster@data <- ug_m2_troposphere_no2.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_ug_m2_troposphere_no2.kriged.raster <- raster(ug_m2_troposphere_no2.kriged.raster)
rm(ug_m2_troposphere_no2_emp_OK, dat.fit, ug_m2_troposphere_no2.kriged, 
   ug_m2_troposphere_no2.kriged.raster)

### ter_pressure
ter_pressure_emp_OK <- gstat::variogram(
  ter_pressure ~ 1, meanValueSpatialDataFrame
)
plot(ter_pressure_emp_OK)
dat.fit  <- fit.variogram(ter_pressure_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(ter_pressure~1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

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
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$ter_pressure)^2) / 
  sum( ( meanValueSpatialDataFrame$ter_pressure - mean(meanValueSpatialDataFrame$ter_pressure) )^2 )
cv.line <- c("ter_pressure.mean",
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$ter_pressure)
ter_pressure.kriged.raster@data <- ter_pressure.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_ter_pressure.kriged.raster <- raster(ter_pressure.kriged.raster)
rm(ter_pressure_emp_OK, dat.fit, ter_pressure.kriged, 
   ter_pressure.kriged.raster)

### temp
temp_emp_OK <- gstat::variogram(
  temp ~ 1, meanValueSpatialDataFrame
)
plot(temp_emp_OK)
dat.fit  <- fit.variogram(temp_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(temp~1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

temp.kriged <- 
  krige(temp~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
temp.kriged@data$Ttest <- temp.kriged@data$var1.pred / 
  temp.kriged@data$var1.var
temp.kriged.raster <- as(temp.kriged, 'SpatialPixelsDataFrame')
temp.kriged.raster <- as(temp.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, temp.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$temp)^2) / 
  sum( ( meanValueSpatialDataFrame$temp - mean(meanValueSpatialDataFrame$temp) )^2 )
cv.line <- c("temp.mean", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$temp)
temp.kriged.raster@data <- temp.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_temp.kriged.raster <- raster(temp.kriged.raster)
rm(temp_emp_OK, dat.fit, temp.kriged, 
   temp.kriged.raster)


### ndvi
ndvi_emp_OK <- gstat::variogram(
  ndvi ~ 1, meanValueSpatialDataFrame
)
plot(ndvi_emp_OK)
dat.fit  <- fit.variogram(ndvi_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(ndvi ~ 1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

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
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$ndvi)^2) / 
  sum( ( meanValueSpatialDataFrame$ndvi - mean(meanValueSpatialDataFrame$ndvi) )^2 )
cv.line <- c("ndvi.mean", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$ndvi)
ndvi.kriged.raster@data <- ndvi.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_ndvi.kriged.raster <- raster(ndvi.kriged.raster)
rm(ndvi_emp_OK, dat.fit, ndvi.kriged, 
   ndvi.kriged.raster)

### precipitation
precipitation_emp_OK <- gstat::variogram(
  precipitation ~ 1, meanValueSpatialDataFrame
)
plot(precipitation_emp_OK)
dat.fit  <- fit.variogram(precipitation_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(precipitation ~ 1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

precipitation.kriged <- 
  krige(precipitation~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
precipitation.kriged@data$Ttest <- precipitation.kriged@data$var1.pred / 
  precipitation.kriged@data$var1.var
precipitation.kriged.raster <- as(precipitation.kriged, 'SpatialPixelsDataFrame')
precipitation.kriged.raster <- as(precipitation.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, precipitation.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$precipitation)^2) / 
  sum( ( meanValueSpatialDataFrame$precipitation - mean(meanValueSpatialDataFrame$precipitation) )^2 )
cv.line <- c("precipitation.mean", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$precipitation)
precipitation.kriged.raster@data <- precipitation.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_precipitation.kriged.raster <- raster(precipitation.kriged.raster)
rm(precipitation_emp_OK, dat.fit, precipitation.kriged, 
   precipitation.kriged.raster)

fun <- F
if(run)
{
  ### NTL
  NTL_emp_OK <- gstat::variogram(
    NTL ~ 1, meanValueSpatialDataFrame
  )
  plot(NTL_emp_OK)
  dat.fit  <- fit.variogram(NTL_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                            vgm(model = "Sph"))
  out <- 
    krige.cv(NTL ~ 1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
             model = dat.fit)
  N <- nrow(meanValueSpatialDataFrame@data)
  R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
  MRSE <- sqrt(mean(out$residual^2))
  MAE <- mean(abs(out$residual))
  correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
  out$pred <- out$observed - out$residual
  coef.reg <- coef(lm(observed ~ pred, data = out))
  
  cv.line <- c("NTL.mean",              
               N, R2_interpolation, MRSE, MAE, 
               correlation_observed_predicted, coef.reg, onpoint.R2)
  print(cv.line)
  kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
  NTL.kriged <- 
    krige(NTL~1, meanValueSpatialDataFrame, coords, 
          model = dat.fit)
  NTL.kriged@data$Ttest <- NTL.kriged@data$var1.pred / 
    NTL.kriged@data$var1.var
  NTL.kriged.raster <- as(NTL.kriged, 'SpatialPixelsDataFrame')
  NTL.kriged.raster <- as(NTL.kriged.raster, "SpatialGridDataFrame")
  predict.value <- over(meanValueSpatialDataFrame, NTL.kriged.raster)
  meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
  # R2
  1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$NTL)^2) / 
    sum( ( meanValueSpatialDataFrame$NTL - mean(meanValueSpatialDataFrame$NTL) )^2 )
  cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$NTL)
  NTL.kriged.raster@data <- NTL.kriged.raster@data %>% dplyr::select(var1.pred)
  MEAN_NTL.kriged.raster <- raster(NTL.kriged.raster)
  rm(NTL_emp_OK, dat.fit, NTL.kriged, 
     NTL.kriged.raster)
}

### PBLH
PBLH_emp_OK <- gstat::variogram(
  PBLH ~ 1, meanValueSpatialDataFrame
)
plot(PBLH_emp_OK)
dat.fit  <- fit.variogram(PBLH_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(PBLH ~ 1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

PBLH.kriged <- 
  krige(PBLH~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
PBLH.kriged@data$Ttest <- PBLH.kriged@data$var1.pred / 
  PBLH.kriged@data$var1.var
PBLH.kriged.raster <- as(PBLH.kriged, 'SpatialPixelsDataFrame')
PBLH.kriged.raster <- as(PBLH.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, PBLH.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$PBLH)^2) / 
  sum( ( meanValueSpatialDataFrame$PBLH - mean(meanValueSpatialDataFrame$PBLH) )^2 )
cv.line <- c("PBLH.mean", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$PBLH)
PBLH.kriged.raster@data <- PBLH.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_PBLH.kriged.raster <- raster(PBLH.kriged.raster)
rm(PBLH_emp_OK, dat.fit, PBLH.kriged, 
   PBLH.kriged.raster)

### Y2016
Y2016_emp_OK <- gstat::variogram(
  Y2016 ~ 1, meanValueSpatialDataFrame
)
plot(Y2016_emp_OK)
dat.fit  <- fit.variogram(Y2016_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(Y2016 ~ 1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

Y2016.kriged <- 
  krige(Y2016~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
Y2016.kriged@data$Ttest <- Y2016.kriged@data$var1.pred / 
  Y2016.kriged@data$var1.var
Y2016.kriged.raster <- as(Y2016.kriged, 'SpatialPixelsDataFrame')
Y2016.kriged.raster <- as(Y2016.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, Y2016.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$Y2016)^2) / 
  sum( ( meanValueSpatialDataFrame$Y2016 - mean(meanValueSpatialDataFrame$Y2016) )^2 )
cv.line <- c("Y2016.mean", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$Y2016)
Y2016.kriged.raster@data <- Y2016.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_Y2016.kriged.raster <- raster(Y2016.kriged.raster)
rm(Y2016_emp_OK, dat.fit, Y2016.kriged, 
   Y2016.kriged.raster)

### Y2017
Y2017_emp_OK <- gstat::variogram(
  Y2017 ~ 1, meanValueSpatialDataFrame
)
plot(Y2017_emp_OK)
dat.fit  <- fit.variogram(Y2017_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(Y2017 ~ 1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

Y2017.kriged <- 
  krige(Y2017~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
Y2017.kriged@data$Ttest <- Y2017.kriged@data$var1.pred / 
  Y2017.kriged@data$var1.var
Y2017.kriged.raster <- as(Y2017.kriged, 'SpatialPixelsDataFrame')
Y2017.kriged.raster <- as(Y2017.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, Y2017.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$Y2017)^2) / 
  sum( ( meanValueSpatialDataFrame$Y2017 - mean(meanValueSpatialDataFrame$Y2017) )^2 )
cv.line <- c("Y2017.mean",
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$Y2017)
Y2017.kriged.raster@data <- Y2017.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_Y2017.kriged.raster <- raster(Y2017.kriged.raster)
rm(Y2017_emp_OK, dat.fit, Y2017.kriged, 
   Y2017.kriged.raster)

### Y2018
Y2018_emp_OK <- gstat::variogram(
  Y2018 ~ 1, meanValueSpatialDataFrame
)
plot(Y2018_emp_OK)
dat.fit  <- fit.variogram(Y2018_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(Y2018 ~ 1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

Y2018.kriged <- 
  krige(Y2018~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
Y2018.kriged@data$Ttest <- Y2018.kriged@data$var1.pred / 
  Y2018.kriged@data$var1.var
Y2018.kriged.raster <- as(Y2018.kriged, 'SpatialPixelsDataFrame')
Y2018.kriged.raster <- as(Y2018.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, Y2018.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$Y2018)^2) / 
  sum( ( meanValueSpatialDataFrame$Y2018 - mean(meanValueSpatialDataFrame$Y2018) )^2 )
cv.line <- c("Y2018.mean",              
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$Y2018)
Y2018.kriged.raster@data <- Y2018.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_Y2018.kriged.raster <- raster(Y2018.kriged.raster)
rm(Y2018_emp_OK, dat.fit, Y2018.kriged, 
   Y2018.kriged.raster)

### Y2019
Y2019_emp_OK <- gstat::variogram(
  Y2019 ~ 1, meanValueSpatialDataFrame
)
plot(Y2019_emp_OK)
dat.fit  <- fit.variogram(Y2019_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(Y2019 ~ 1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

Y2019.kriged <- 
  krige(Y2019~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
Y2019.kriged@data$Ttest <- Y2019.kriged@data$var1.pred / 
  Y2019.kriged@data$var1.var
Y2019.kriged.raster <- as(Y2019.kriged, 'SpatialPixelsDataFrame')
Y2019.kriged.raster <- as(Y2019.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, Y2019.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$Y2019)^2) / 
  sum( ( meanValueSpatialDataFrame$Y2019 - mean(meanValueSpatialDataFrame$Y2019) )^2 )
cv.line <- c("Y2019.mean", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$Y2019)
Y2019.kriged.raster@data <- Y2019.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_Y2019.kriged.raster <- raster(Y2019.kriged.raster)
rm(Y2019_emp_OK, dat.fit, Y2019.kriged, 
   Y2019.kriged.raster)

### Y2020
Y2020_emp_OK <- gstat::variogram(
  Y2020 ~ 1, meanValueSpatialDataFrame
)
plot(Y2020_emp_OK)
dat.fit  <- fit.variogram(Y2020_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(Y2020 ~ 1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

Y2020.kriged <- 
  krige(Y2020~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
Y2020.kriged@data$Ttest <- Y2020.kriged@data$var1.pred / 
  Y2020.kriged@data$var1.var
Y2020.kriged.raster <- as(Y2020.kriged, 'SpatialPixelsDataFrame')
Y2020.kriged.raster <- as(Y2020.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, Y2020.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$Y2020)^2) / 
  sum( ( meanValueSpatialDataFrame$Y2020 - mean(meanValueSpatialDataFrame$Y2020) )^2 )
cv.line <- c("Y2020.mean",
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$Y2020)
Y2020.kriged.raster@data <- Y2020.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_Y2020.kriged.raster <- raster(Y2020.kriged.raster)
rm(Y2020_emp_OK, dat.fit, Y2020.kriged, 
   Y2020.kriged.raster)

### Y2021
Y2021_emp_OK <- gstat::variogram(
  Y2021 ~ 1, meanValueSpatialDataFrame
)
plot(Y2021_emp_OK)
dat.fit  <- fit.variogram(Y2021_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(Y2021 ~ 1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
N <- nrow(meanValueSpatialDataFrame@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

Y2021.kriged <- 
  krige(Y2021~1, meanValueSpatialDataFrame, coords, 
        model = dat.fit)
Y2021.kriged@data$Ttest <- Y2021.kriged@data$var1.pred / 
  Y2021.kriged@data$var1.var
Y2021.kriged.raster <- as(Y2021.kriged, 'SpatialPixelsDataFrame')
Y2021.kriged.raster <- as(Y2021.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(meanValueSpatialDataFrame, Y2021.kriged.raster)
meanValueSpatialDataFrame$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$Y2021)^2) / 
  sum( ( meanValueSpatialDataFrame$Y2021 - mean(meanValueSpatialDataFrame$Y2021) )^2 )
cv.line <- c("Y2021.mean", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$Y2021)
Y2021.kriged.raster@data <- Y2021.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_Y2021.kriged.raster <- raster(Y2021.kriged.raster)
rm(Y2021_emp_OK, dat.fit, Y2021.kriged, 
   Y2021.kriged.raster)

colnames(kriging.cv.mean.dataset) <- c("Variable", "N", "R2", "MRSE", "MAE", 
                                       "r", "Intercept", "Slope", "onpoint.R2")

save(MEAN_no2_measured_ug.m3.kriged.raster,
     MEAN_ug_m2_troposphere_no2.kriged.raster,
     MEAN_ndvi.kriged.raster, MEAN_temp.kriged.raster, 
     MEAN_PBLH.kriged.raster, MEAN_precipitation.kriged.raster, 
     MEAN_ter_pressure.kriged.raster, MEAN_Y2016.kriged.raster,
     MEAN_Y2017.kriged.raster, MEAN_Y2018.kriged.raster,
     MEAN_Y2019.kriged.raster, MEAN_Y2020.kriged.raster,
     MEAN_Y2021.kriged.raster, file = "05_CoefficientRaster/MEAN_raster.RData")

save(kriging.cv.mean.dataset, file = "04_Results/krigingMeanResult.RData")
write.csv(kriging.cv.mean.dataset, file = "08_Tables/krigingCVMean.csv")
