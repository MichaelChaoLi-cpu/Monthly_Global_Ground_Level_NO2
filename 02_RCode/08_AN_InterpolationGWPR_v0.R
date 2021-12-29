# Author: M.L.

# input: GWPR_FEM_CV_F_result.Rdata
# GWPR_FEM_CV_F_result.Rdata: "GWPR.FEM.CV.F.result" GWPR result with 2.25 fixed distance bandwidth 
#                                                    based on the FEM

# output: COEF_raster.RData
# COEF_raster.RData: "CO_dayTimeTemperature.kriged.raster" interpolation of day time temperature coefficient
#                                                          based on ordinary kriging 
# COEF_raster.RData: "CO_humidity.kriged.raster" interpolation of humidity coefficient
#                                                based on ordinary kriging
# COEF_raster.RData: "CO_mg_m2_troposphere_no2.kriged.raster" interpolation of troposheric no2 coefficient
#                                                             based on ordinary kriging
# COEF_raster.RData: "CO_ndvi.kriged.raster" interpolation of ndvi coefficient based on ordinary kriging (OK)
# COEF_raster.RData: "CO_nightTimeTemperature.kriged.raster" interpolation of night temperature coefficient (OK)
# COEF_raster.RData: "CO_NTL.kriged.raster" interpolation of night time light coefficient (OK)
# COEF_raster.RData: "CO_PBLH.kriged.raster" interpolation of PBLH coefficient (OK)
# COEF_raster.RData: "CO_precipitation.kriged.raster" interpolation of precipitation coefficient (OK)
# COEF_raster.RData: "CO_speedwind.kriged.raster" interpolation of wind speed coefficient (OK)
# COEF_raster.RData: "CO_ter_pressure.kriged.raster" interpolation of air pressure coefficient (OK)
# COEF_raster.RData: "CO_Y2016.kriged.raster" interpolation of 2016 dummy variable coefficient (OK)
# COEF_raster.RData: "CO_Y2017.kriged.raster" interpolation of 2017 dummy variable coefficient (OK)
# COEF_raster.RData: "CO_Y2018.kriged.raster" interpolation of 2018 dummy variable coefficient (OK)
# COEF_raster.RData: "CO_Y2019.kriged.raster" interpolation of 2019 dummy variable coefficient (OK)
# COEF_raster.RData: "CO_Y2020.kriged.raster" interpolation of 2020 dummy variable coefficient (OK)
# COEF_raster.RData: "CO_Y2021.kriged.raster" interpolation of 2021 dummy variable coefficient (OK)



# end

library(tidyverse)
library(gstat)
library(sp) 
library(raster)
library(dplyr)
library(tmap)
library(automap)

setwd("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub")
load("04_Results/GWPR_FEM_CV_F_result.Rdata")
# make the result into point
GWPR.point.dataset <- GWPR.FEM.CV.F.result$SDF

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
# check normality
shapiro.test(GWPR.point.dataset$mg_m2_troposphere_no2)

bubble(GWPR.point.dataset, "mg_m2_troposphere_no2", fill = FALSE, maxsize = 2, identify = FALSE)

coords <- addcoord(nx,xmin,xsize,ny,ymin,ysize)

# the tutorial of this interpolation
# https://github.com/GeostatsGuy/geostatsr/blob/master/kriging_demo_Rnotebook.ipynb
GWPR.point.dataset <- spTransform(GWPR.point.dataset, CRS(proj))
coords@proj4string <- CRS(proj)
summary(coords) 
class(coords)

mg_m2_troposphere_no2.idw.2.0 = idw(mg_m2_troposphere_no2~1, idp = 2.0, GWPR.point.dataset, coords)
class(mg_m2_troposphere_no2.idw.2.0)                        # check the inverse distance object
cuts = c(-0.1, -0.06, -0.03, 0, 0.03, 0.06, 0.09, 0.12)

spplot(mg_m2_troposphere_no2.idw.2.0["var1.pred"],main = "mg_m2_troposphere_no2 Inverse Distance (p=2.0)",
       key.space = "right", cuts = cuts, xlab = "Degree", ylab = "Degree")


# Fine-tuning the interpolation
# the following procedure are learned from https://mgimond.github.io/Spatial/interpolation-in-r.html
IDW.out <- vector(length = length(GWPR.point.dataset))
for (i in 1:length(GWPR.point.dataset)) {
  IDW.out[i] <- idw(mg_m2_troposphere_no2~1, GWPR.point.dataset[-i,],
                    GWPR.point.dataset[i,], idp = 2.0)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ GWPR.point.dataset$mg_m2_troposphere_no2, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ GWPR.point.dataset$mg_m2_troposphere_no2), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
cor.test(IDW.out, GWPR.point.dataset$mg_m2_troposphere_no2)
# RMSE
sqrt( sum((IDW.out - GWPR.point.dataset$mg_m2_troposphere_no2)^2) / length(GWPR.point.dataset))
# R2
1 - sum( (IDW.out - GWPR.point.dataset$mg_m2_troposphere_no2)^2) / 
  sum( ( GWPR.point.dataset$mg_m2_troposphere_no2)^2 )


## 
#-------------------------------------------kriging----------------------------------
## https://swilke-geoscience.net/post/2020-09-10-kriging_with_r/kriging/
### mg_m2_troposphere_no2
mg_m2_troposphere_no2_emp_OK <- gstat::variogram(
  mg_m2_troposphere_no2 ~ 1, GWPR.point.dataset
)
plot(mg_m2_troposphere_no2_emp_OK)
# Compute the variogram model by passing the nugget, sill and range values
# to fit.variogram() via the vgm() function.
dat.fit  <- fit.variogram(mg_m2_troposphere_no2_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
mg_m2_troposphere_no2.kriged <- 
  krige(mg_m2_troposphere_no2~1, GWPR.point.dataset, coords, 
        model = dat.fit)
mg_m2_troposphere_no2.kriged@data$Ttest <- mg_m2_troposphere_no2.kriged@data$var1.pred / 
  mg_m2_troposphere_no2.kriged@data$var1.var
mg_m2_troposphere_no2.kriged.raster <- as(mg_m2_troposphere_no2.kriged, 'SpatialPixelsDataFrame')
mg_m2_troposphere_no2.kriged.raster <- as(mg_m2_troposphere_no2.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, mg_m2_troposphere_no2.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$mg_m2_troposphere_no2)^2) / 
  sum( ( GWPR.point.dataset$mg_m2_troposphere_no2 - mean(GWPR.point.dataset$mg_m2_troposphere_no2) )^2 )
mg_m2_troposphere_no2.kriged.raster@data <- mg_m2_troposphere_no2.kriged.raster@data %>% dplyr::select(var1.pred)
CO_mg_m2_troposphere_no2.kriged.raster <- raster(mg_m2_troposphere_no2.kriged.raster)

### ter_pressure
ter_pressure_emp_OK <- gstat::variogram(
  ter_pressure ~ 1, GWPR.point.dataset
)
plot(ter_pressure_emp_OK)
dat.fit  <- fit.variogram(ter_pressure_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
ter_pressure.kriged <- 
  krige(ter_pressure~1, GWPR.point.dataset, coords, 
        model = dat.fit)
ter_pressure.kriged@data$Ttest <- ter_pressure.kriged@data$var1.pred / 
  ter_pressure.kriged@data$var1.var
ter_pressure.kriged.raster <- as(ter_pressure.kriged, 'SpatialPixelsDataFrame')
ter_pressure.kriged.raster <- as(ter_pressure.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, ter_pressure.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$ter_pressure)^2) / 
  sum( ( GWPR.point.dataset$ter_pressure - mean(GWPR.point.dataset$ter_pressure) )^2 )
ter_pressure.kriged.raster@data <- ter_pressure.kriged.raster@data %>% dplyr::select(var1.pred)
CO_ter_pressure.kriged.raster <- raster(ter_pressure.kriged.raster)

### dayTimeTemperature
dayTimeTemperature_emp_OK <- gstat::variogram(
  dayTimeTemperature ~ 1, GWPR.point.dataset
)
plot(dayTimeTemperature_emp_OK)
dat.fit  <- fit.variogram(dayTimeTemperature_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
dayTimeTemperature.kriged <- 
  krige(dayTimeTemperature~1, GWPR.point.dataset, coords, 
        model = dat.fit)
dayTimeTemperature.kriged@data$Ttest <- dayTimeTemperature.kriged@data$var1.pred / 
  dayTimeTemperature.kriged@data$var1.var
dayTimeTemperature.kriged.raster <- as(dayTimeTemperature.kriged, 'SpatialPixelsDataFrame')
dayTimeTemperature.kriged.raster <- as(dayTimeTemperature.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, dayTimeTemperature.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$dayTimeTemperature)^2) / 
  sum( ( GWPR.point.dataset$dayTimeTemperature - mean(GWPR.point.dataset$dayTimeTemperature) )^2 )
dayTimeTemperature.kriged.raster@data <- dayTimeTemperature.kriged.raster@data %>% dplyr::select(var1.pred)
CO_dayTimeTemperature.kriged.raster <- raster(dayTimeTemperature.kriged.raster)
rm(dayTimeTemperature_emp_OK, dat.fit, dayTimeTemperature.kriged, 
   dayTimeTemperature.kriged.raster)

### nightTimeTemperature
nightTimeTemperature_emp_OK <- gstat::variogram(
  nightTimeTemperature ~ 1, GWPR.point.dataset
)
plot(nightTimeTemperature_emp_OK)
dat.fit  <- fit.variogram(nightTimeTemperature_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
nightTimeTemperature.kriged <- 
  krige(nightTimeTemperature~1, GWPR.point.dataset, coords, 
        model = dat.fit)
nightTimeTemperature.kriged@data$Ttest <- nightTimeTemperature.kriged@data$var1.pred / 
  nightTimeTemperature.kriged@data$var1.var
nightTimeTemperature.kriged.raster <- as(nightTimeTemperature.kriged, 'SpatialPixelsDataFrame')
nightTimeTemperature.kriged.raster <- as(nightTimeTemperature.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, nightTimeTemperature.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$nightTimeTemperature)^2) / 
  sum( ( GWPR.point.dataset$nightTimeTemperature - mean(GWPR.point.dataset$nightTimeTemperature) )^2 )
nightTimeTemperature.kriged.raster@data <- nightTimeTemperature.kriged.raster@data %>% dplyr::select(var1.pred)
CO_nightTimeTemperature.kriged.raster <- raster(nightTimeTemperature.kriged.raster)
rm(nightTimeTemperature_emp_OK, dat.fit, nightTimeTemperature.kriged, 
   nightTimeTemperature.kriged.raster)

### ndvi
ndvi_emp_OK <- gstat::variogram(
  ndvi ~ 1, GWPR.point.dataset
)
plot(ndvi_emp_OK)
dat.fit  <- fit.variogram(ndvi_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
ndvi.kriged <- 
  krige(ndvi~1, GWPR.point.dataset, coords, 
        model = dat.fit)
ndvi.kriged@data$Ttest <- ndvi.kriged@data$var1.pred / 
  ndvi.kriged@data$var1.var
ndvi.kriged.raster <- as(ndvi.kriged, 'SpatialPixelsDataFrame')
ndvi.kriged.raster <- as(ndvi.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, ndvi.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$ndvi)^2) / 
  sum( ( GWPR.point.dataset$ndvi - mean(GWPR.point.dataset$ndvi) )^2 )
ndvi.kriged.raster@data <- ndvi.kriged.raster@data %>% dplyr::select(var1.pred)
CO_ndvi.kriged.raster <- raster(ndvi.kriged.raster)
rm(ndvi_emp_OK, dat.fit, ndvi.kriged, 
   ndvi.kriged.raster)

### humidity
humidity_emp_OK <- gstat::variogram(
  humidity ~ 1, GWPR.point.dataset
)
plot(humidity_emp_OK)
dat.fit  <- fit.variogram(humidity_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
humidity.kriged <- 
  krige(humidity~1, GWPR.point.dataset, coords, 
        model = dat.fit)
humidity.kriged@data$Ttest <- humidity.kriged@data$var1.pred / 
  humidity.kriged@data$var1.var
humidity.kriged.raster <- as(humidity.kriged, 'SpatialPixelsDataFrame')
humidity.kriged.raster <- as(humidity.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, humidity.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$humidity)^2) / 
  sum( ( GWPR.point.dataset$humidity - mean(GWPR.point.dataset$humidity) )^2 )
humidity.kriged.raster@data <- humidity.kriged.raster@data %>% dplyr::select(var1.pred)
CO_humidity.kriged.raster <- raster(humidity.kriged.raster)
rm(humidity_emp_OK, dat.fit, humidity.kriged, 
   humidity.kriged.raster)

### precipitation
precipitation_emp_OK <- gstat::variogram(
  precipitation ~ 1, GWPR.point.dataset
)
plot(precipitation_emp_OK)
dat.fit  <- fit.variogram(precipitation_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
precipitation.kriged <- 
  krige(precipitation~1, GWPR.point.dataset, coords, 
        model = dat.fit)
precipitation.kriged@data$Ttest <- precipitation.kriged@data$var1.pred / 
  precipitation.kriged@data$var1.var
precipitation.kriged.raster <- as(precipitation.kriged, 'SpatialPixelsDataFrame')
precipitation.kriged.raster <- as(precipitation.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, precipitation.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$precipitation)^2) / 
  sum( ( GWPR.point.dataset$precipitation - mean(GWPR.point.dataset$precipitation) )^2 )
precipitation.kriged.raster@data <- precipitation.kriged.raster@data %>% dplyr::select(var1.pred)
CO_precipitation.kriged.raster <- raster(precipitation.kriged.raster)
rm(precipitation_emp_OK, dat.fit, precipitation.kriged, 
   precipitation.kriged.raster)

### NTL
NTL_emp_OK <- gstat::variogram(
  NTL ~ 1, GWPR.point.dataset
)
plot(NTL_emp_OK)
dat.fit  <- fit.variogram(NTL_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
NTL.kriged <- 
  krige(NTL~1, GWPR.point.dataset, coords, 
        model = dat.fit)
NTL.kriged@data$Ttest <- NTL.kriged@data$var1.pred / 
  NTL.kriged@data$var1.var
NTL.kriged.raster <- as(NTL.kriged, 'SpatialPixelsDataFrame')
NTL.kriged.raster <- as(NTL.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, NTL.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$NTL)^2) / 
  sum( ( GWPR.point.dataset$NTL - mean(GWPR.point.dataset$NTL) )^2 )
NTL.kriged.raster@data <- NTL.kriged.raster@data %>% dplyr::select(var1.pred)
CO_NTL.kriged.raster <- raster(NTL.kriged.raster)
rm(NTL_emp_OK, dat.fit, NTL.kriged, 
   NTL.kriged.raster)

### speedwind
speedwind_emp_OK <- gstat::variogram(
  speedwind ~ 1, GWPR.point.dataset
)
plot(speedwind_emp_OK)
dat.fit  <- fit.variogram(speedwind_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
speedwind.kriged <- 
  krige(speedwind~1, GWPR.point.dataset, coords, 
        model = dat.fit)
speedwind.kriged@data$Ttest <- speedwind.kriged@data$var1.pred / 
  speedwind.kriged@data$var1.var
speedwind.kriged.raster <- as(speedwind.kriged, 'SpatialPixelsDataFrame')
speedwind.kriged.raster <- as(speedwind.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, speedwind.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$speedwind)^2) / 
  sum( ( GWPR.point.dataset$speedwind - mean(GWPR.point.dataset$speedwind) )^2 )
speedwind.kriged.raster@data <- speedwind.kriged.raster@data %>% dplyr::select(var1.pred)
CO_speedwind.kriged.raster <- raster(speedwind.kriged.raster)
rm(speedwind_emp_OK, dat.fit, speedwind.kriged, 
   speedwind.kriged.raster)

### PBLH
PBLH_emp_OK <- gstat::variogram(
  PBLH ~ 1, GWPR.point.dataset
)
plot(PBLH_emp_OK)
dat.fit  <- fit.variogram(PBLH_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
PBLH.kriged <- 
  krige(PBLH~1, GWPR.point.dataset, coords, 
        model = dat.fit)
PBLH.kriged@data$Ttest <- PBLH.kriged@data$var1.pred / 
  PBLH.kriged@data$var1.var
PBLH.kriged.raster <- as(PBLH.kriged, 'SpatialPixelsDataFrame')
PBLH.kriged.raster <- as(PBLH.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, PBLH.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$PBLH)^2) / 
  sum( ( GWPR.point.dataset$PBLH - mean(GWPR.point.dataset$PBLH) )^2 )
PBLH.kriged.raster@data <- PBLH.kriged.raster@data %>% dplyr::select(var1.pred)
CO_PBLH.kriged.raster <- raster(PBLH.kriged.raster)
rm(PBLH_emp_OK, dat.fit, PBLH.kriged, 
   PBLH.kriged.raster)

### Y2016
Y2016_emp_OK <- gstat::variogram(
  Y2016 ~ 1, GWPR.point.dataset
)
plot(Y2016_emp_OK)
dat.fit  <- fit.variogram(Y2016_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
Y2016.kriged <- 
  krige(Y2016~1, GWPR.point.dataset, coords, 
        model = dat.fit)
Y2016.kriged@data$Ttest <- Y2016.kriged@data$var1.pred / 
  Y2016.kriged@data$var1.var
Y2016.kriged.raster <- as(Y2016.kriged, 'SpatialPixelsDataFrame')
Y2016.kriged.raster <- as(Y2016.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, Y2016.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2016)^2) / 
  sum( ( GWPR.point.dataset$Y2016 - mean(GWPR.point.dataset$Y2016) )^2 )
Y2016.kriged.raster@data <- Y2016.kriged.raster@data %>% dplyr::select(var1.pred)
CO_Y2016.kriged.raster <- raster(Y2016.kriged.raster)
rm(Y2016_emp_OK, dat.fit, Y2016.kriged, 
   Y2016.kriged.raster)

### Y2017
Y2017_emp_OK <- gstat::variogram(
  Y2017 ~ 1, GWPR.point.dataset
)
plot(Y2017_emp_OK)
dat.fit  <- fit.variogram(Y2017_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
Y2017.kriged <- 
  krige(Y2017~1, GWPR.point.dataset, coords, 
        model = dat.fit)
Y2017.kriged@data$Ttest <- Y2017.kriged@data$var1.pred / 
  Y2017.kriged@data$var1.var
Y2017.kriged.raster <- as(Y2017.kriged, 'SpatialPixelsDataFrame')
Y2017.kriged.raster <- as(Y2017.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, Y2017.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2017)^2) / 
  sum( ( GWPR.point.dataset$Y2017 - mean(GWPR.point.dataset$Y2017) )^2 )
Y2017.kriged.raster@data <- Y2017.kriged.raster@data %>% dplyr::select(var1.pred)
CO_Y2017.kriged.raster <- raster(Y2017.kriged.raster)
rm(Y2017_emp_OK, dat.fit, Y2017.kriged, 
   Y2017.kriged.raster)

### Y2018
Y2018_emp_OK <- gstat::variogram(
  Y2018 ~ 1, GWPR.point.dataset
)
plot(Y2018_emp_OK)
dat.fit  <- fit.variogram(Y2018_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
Y2018.kriged <- 
  krige(Y2018~1, GWPR.point.dataset, coords, 
        model = dat.fit)
Y2018.kriged@data$Ttest <- Y2018.kriged@data$var1.pred / 
  Y2018.kriged@data$var1.var
Y2018.kriged.raster <- as(Y2018.kriged, 'SpatialPixelsDataFrame')
Y2018.kriged.raster <- as(Y2018.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, Y2018.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2018)^2) / 
  sum( ( GWPR.point.dataset$Y2018 - mean(GWPR.point.dataset$Y2018) )^2 )
Y2018.kriged.raster@data <- Y2018.kriged.raster@data %>% dplyr::select(var1.pred)
CO_Y2018.kriged.raster <- raster(Y2018.kriged.raster)
rm(Y2018_emp_OK, dat.fit, Y2018.kriged, 
   Y2018.kriged.raster)

### Y2019
Y2019_emp_OK <- gstat::variogram(
  Y2019 ~ 1, GWPR.point.dataset
)
plot(Y2019_emp_OK)
dat.fit  <- fit.variogram(Y2019_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
Y2019.kriged <- 
  krige(Y2019~1, GWPR.point.dataset, coords, 
        model = dat.fit)
Y2019.kriged@data$Ttest <- Y2019.kriged@data$var1.pred / 
  Y2019.kriged@data$var1.var
Y2019.kriged.raster <- as(Y2019.kriged, 'SpatialPixelsDataFrame')
Y2019.kriged.raster <- as(Y2019.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, Y2019.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2019)^2) / 
  sum( ( GWPR.point.dataset$Y2019 - mean(GWPR.point.dataset$Y2019) )^2 )
Y2019.kriged.raster@data <- Y2019.kriged.raster@data %>% dplyr::select(var1.pred)
CO_Y2019.kriged.raster <- raster(Y2019.kriged.raster)
rm(Y2019_emp_OK, dat.fit, Y2019.kriged, 
   Y2019.kriged.raster)

### Y2020
Y2020_emp_OK <- gstat::variogram(
  Y2020 ~ 1, GWPR.point.dataset
)
plot(Y2020_emp_OK)
dat.fit  <- fit.variogram(Y2020_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
Y2020.kriged <- 
  krige(Y2020~1, GWPR.point.dataset, coords, 
        model = dat.fit)
Y2020.kriged@data$Ttest <- Y2020.kriged@data$var1.pred / 
  Y2020.kriged@data$var1.var
Y2020.kriged.raster <- as(Y2020.kriged, 'SpatialPixelsDataFrame')
Y2020.kriged.raster <- as(Y2020.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, Y2020.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2020)^2) / 
  sum( ( GWPR.point.dataset$Y2020 - mean(GWPR.point.dataset$Y2020) )^2 )
Y2020.kriged.raster@data <- Y2020.kriged.raster@data %>% dplyr::select(var1.pred)
CO_Y2020.kriged.raster <- raster(Y2020.kriged.raster)
rm(Y2020_emp_OK, dat.fit, Y2020.kriged, 
   Y2020.kriged.raster)

### Y2021
Y2021_emp_OK <- gstat::variogram(
  Y2021 ~ 1, GWPR.point.dataset
)
plot(Y2021_emp_OK)
dat.fit  <- fit.variogram(Y2021_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
Y2021.kriged <- 
  krige(Y2021~1, GWPR.point.dataset, coords, 
        model = dat.fit)
Y2021.kriged@data$Ttest <- Y2021.kriged@data$var1.pred / 
  Y2021.kriged@data$var1.var
Y2021.kriged.raster <- as(Y2021.kriged, 'SpatialPixelsDataFrame')
Y2021.kriged.raster <- as(Y2021.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, Y2021.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2021)^2) / 
  sum( ( GWPR.point.dataset$Y2021 - mean(GWPR.point.dataset$Y2021) )^2 )
Y2021.kriged.raster@data <- Y2021.kriged.raster@data %>% dplyr::select(var1.pred)
CO_Y2021.kriged.raster <- raster(Y2021.kriged.raster)
rm(Y2021_emp_OK, dat.fit, Y2021.kriged, 
   Y2021.kriged.raster)

save(CO_dayTimeTemperature.kriged.raster, CO_humidity.kriged.raster, CO_mg_m2_troposphere_no2.kriged.raster,
     CO_ndvi.kriged.raster, CO_nightTimeTemperature.kriged.raster, CO_NTL.kriged.raster,
     CO_PBLH.kriged.raster, CO_precipitation.kriged.raster, CO_speedwind.kriged.raster,
     CO_ter_pressure.kriged.raster, CO_Y2016.kriged.raster, CO_Y2017.kriged.raster,
     CO_Y2018.kriged.raster, CO_Y2019.kriged.raster, CO_Y2020.kriged.raster,
     CO_Y2021.kriged.raster, file = "05_CoefficientRaster/COEF_raster.RData")
