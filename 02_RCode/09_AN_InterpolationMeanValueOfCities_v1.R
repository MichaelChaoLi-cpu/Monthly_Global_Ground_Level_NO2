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

# output: MEAN_raster.RData
# MEAN_raster.RData: "MEAN_mg_m2_troposphere_no2.kriged.raster" interpolation of troposheric no2 mean value
#                                                               in the usedDataset based on ordinary kriging (OK)
# MEAN_raster.RData: "MEAN_ndvi.kriged.raster" interpolation of ndvi mean value (OK)
# MEAN_raster.RData: "MEAN_temp.kriged.raster" interpolation of temperature mean value (OK)
# MEAN_raster.RData: "MEAN_NTL.kriged.raster" interpolation of night time light mean value (OK)
# MEAN_raster.RData: "MEAN_PBLH.kriged.raster" interpolation of PBLH mean value (OK)
# MEAN_raster.RData: "MEAN_precipitation.kriged.raster" interpolation of precipitation mean value (OK)
# MEAN_raster.RData: "MEAN_ter_pressure.kriged.raster" interpolation of air pressure mean value (OK)

# output: krigingMeanResult.RData
# krigingMeanResult.RData: "mean_error" mean error
# krigingMeanResult.RData: "MSPE" mean square error
# krigingMeanResult.RData: "MSNE" Mean square normalized error
# krigingMeanResult.RData: "CoOP" correlation observed and predicted
# krigingMeanResult.RData: "CoPR" correlation predicted and residual
# krigingMeanResult.RData: "R2" centered R2

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
  no2_measured_mg.m3 ~ mg_m2_troposphere_no2 + ter_pressure + temp +
  ndvi +  precipitation + NTL + PBLH + 
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

### no2_measured_mg.m3
no2_measured_mg.m3_emp_OK <- gstat::variogram(
  no2_measured_mg.m3 ~ 1, meanValueSpatialDataFrame
)
plot(no2_measured_mg.m3_emp_OK)
dat.fit  <- fit.variogram(no2_measured_mg.m3_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(no2_measured_mg.m3~1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
mean_error <- mean(out$residual)
MSPE <- mean(out$residual^2)
Mean_square_normalized_error <- mean(out$zscore^2)
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
correlation_predicted_residuals <- cor(out$observed - out$residual, out$residual)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed)^2 )
cv.line <- c("no2_measured_mg.m3.mean", mean_error, MSPE, Mean_square_normalized_error,
             correlation_observed_predicted, correlation_predicted_residuals, R2_interpolation)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
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
out <- 
  krige.cv(mg_m2_troposphere_no2~1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
mean_error <- mean(out$residual)
MSPE <- mean(out$residual^2)
Mean_square_normalized_error <- mean(out$zscore^2)
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
correlation_predicted_residuals <- cor(out$observed - out$residual, out$residual)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed)^2 )
cv.line <- c("mg_m2_troposphere_no2.mean", mean_error, MSPE, Mean_square_normalized_error,
             correlation_observed_predicted, correlation_predicted_residuals, R2_interpolation)
print(cv.line)
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
out <- 
  krige.cv(ter_pressure~1, meanValueSpatialDataFrame, meanValueSpatialDataFrame@data, 
           model = dat.fit)
mean_error <- mean(out$residual)
MSPE <- mean(out$residual^2)
Mean_square_normalized_error <- mean(out$zscore^2)
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
correlation_predicted_residuals <- cor(out$observed - out$residual, out$residual)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed)^2 )
cv.line <- c("ter_pressure.mean", mean_error, MSPE, Mean_square_normalized_error,
             correlation_observed_predicted, correlation_predicted_residuals, R2_interpolation)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
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
mean_error <- mean(out$residual)
MSPE <- mean(out$residual^2)
Mean_square_normalized_error <- mean(out$zscore^2)
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
correlation_predicted_residuals <- cor(out$observed - out$residual, out$residual)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed)^2 )
cv.line <- c("temp.mean", mean_error, MSPE, Mean_square_normalized_error,
             correlation_observed_predicted, correlation_predicted_residuals, R2_interpolation)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
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
1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$temp)^2) / 
  sum( ( meanValueSpatialDataFrame$temp - mean(meanValueSpatialDataFrame$temp) )^2 )
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
mean_error <- mean(out$residual)
MSPE <- mean(out$residual^2)
Mean_square_normalized_error <- mean(out$zscore^2)
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
correlation_predicted_residuals <- cor(out$observed - out$residual, out$residual)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed)^2 )
cv.line <- c("ndvi.mean", mean_error, MSPE, Mean_square_normalized_error,
             correlation_observed_predicted, correlation_predicted_residuals, R2_interpolation)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
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
mean_error <- mean(out$residual)
MSPE <- mean(out$residual^2)
Mean_square_normalized_error <- mean(out$zscore^2)
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
correlation_predicted_residuals <- cor(out$observed - out$residual, out$residual)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed)^2 )
cv.line <- c("precipitation.mean", mean_error, MSPE, Mean_square_normalized_error,
             correlation_observed_predicted, correlation_predicted_residuals, R2_interpolation)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
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
1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$precipitation)^2) / 
  sum( ( meanValueSpatialDataFrame$precipitation - mean(meanValueSpatialDataFrame$precipitation) )^2 )
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$precipitation)
precipitation.kriged.raster@data <- precipitation.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_precipitation.kriged.raster <- raster(precipitation.kriged.raster)
rm(precipitation_emp_OK, dat.fit, precipitation.kriged, 
   precipitation.kriged.raster)

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
mean_error <- mean(out$residual)
MSPE <- mean(out$residual^2)
Mean_square_normalized_error <- mean(out$zscore^2)
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
correlation_predicted_residuals <- cor(out$observed - out$residual, out$residual)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed)^2 )
cv.line <- c("NTL.mean", mean_error, MSPE, Mean_square_normalized_error,
             correlation_observed_predicted, correlation_predicted_residuals, R2_interpolation)
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
mean_error <- mean(out$residual)
MSPE <- mean(out$residual^2)
Mean_square_normalized_error <- mean(out$zscore^2)
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
correlation_predicted_residuals <- cor(out$observed - out$residual, out$residual)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed)^2 )
cv.line <- c("PBLH.mean", mean_error, MSPE, Mean_square_normalized_error,
             correlation_observed_predicted, correlation_predicted_residuals, R2_interpolation)
print(cv.line)
kriging.cv.mean.dataset <- rbind(kriging.cv.mean.dataset, cv.line)
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
1 - sum( (meanValueSpatialDataFrame$predict.value - meanValueSpatialDataFrame$PBLH)^2) / 
  sum( ( meanValueSpatialDataFrame$PBLH - mean(meanValueSpatialDataFrame$PBLH) )^2 )
cor.test(meanValueSpatialDataFrame$predict.value, meanValueSpatialDataFrame$PBLH)
PBLH.kriged.raster@data <- PBLH.kriged.raster@data %>% dplyr::select(var1.pred)
MEAN_PBLH.kriged.raster <- raster(PBLH.kriged.raster)
rm(PBLH_emp_OK, dat.fit, PBLH.kriged, 
   PBLH.kriged.raster)

colnames(kriging.cv.mean.dataset) <- c("Variable", "mean_error", "MSPE", "MSNE",
                                         "CoOP", "CoPR", "R2")

save(MEAN_no2_measured_mg.m3.kriged.raster,
     MEAN_mg_m2_troposphere_no2.kriged.raster,
     MEAN_ndvi.kriged.raster, MEAN_temp.kriged.raster, MEAN_NTL.kriged.raster,
     MEAN_PBLH.kriged.raster, MEAN_precipitation.kriged.raster, 
     MEAN_ter_pressure.kriged.raster, file = "05_CoefficientRaster/MEAN_raster.RData")

save(kriging.cv.mean.dataset, file = "04_Results/krigingMeanResult.RData")