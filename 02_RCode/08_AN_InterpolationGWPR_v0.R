# Author: M.L.

# input: GWPR_FEM_CV_F_result.Rdata
# GWPR_FEM_CV_F_result.Rdata: "GWPR.FEM.CV.F.result" GWPR result with 2.25 fixed distance bandwidth 
#                                                    based on the FEM

# output: COEF_raster.RData
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

# output: krigingCVResult.RData
# krigingCVResult.RData: "N" number of observation
# krigingCVResult.RData: "R2" R2
# krigingCVResult.RData: "RMSE" root mean square error
# krigingCVResult.RData: "MAE" mean absolute error
# krigingCVResult.RData: "r" correlation coefficient
# krigingCVResult.RData: "Intercept" intercept of regression
# krigingCVResult.RData: "Slope" slope of regression
# krigingCVResult.RData: "onPoint.R2" R2 after interpolation

# Note: in this version, "NTL" is dropped. (COEF_raster.RData: "CO_NTL.kriged.raster" 
#                                           interpolation of night time light coefficient (OK))
# Note: we begin to test the results with adaptive bw, now we set 7

# end

library(tidyverse)
library(gstat)
library(sp) 
library(raster)
library(dplyr)
library(tmap)
library(automap)

setwd("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub")
load("04_Results/GWPR_FEM_CV_A_result.Rdata")
# make the result into point
GWPR.point.dataset <- GWPR.FEM.CV.A.result$SDF

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
shapiro.test(GWPR.point.dataset$ug_m2_troposphere_no2)

bubble(GWPR.point.dataset, "ug_m2_troposphere_no2", fill = FALSE, maxsize = 2, identify = FALSE)

coords <- addcoord(nx,xmin,xsize,ny,ymin,ysize)

# the tutorial of this interpolation
# https://github.com/GeostatsGuy/geostatsr/blob/master/kriging_demo_Rnotebook.ipynb
GWPR.point.dataset <- spTransform(GWPR.point.dataset, CRS(proj))
coords@proj4string <- CRS(proj)
summary(coords) 
class(coords)

kriging.cv.result.dataset <- data.frame(Doubles = double(), Ints = integer(),
                                        Factors = factor(), Logicals = logical(),
                                        Characters = character(), stringsAsFactors = FALSE)

#IDW test
RUN <- F
if(RUN) {
  ug_m2_troposphere_no2.idw.2.0 = idw(ug_m2_troposphere_no2~1, idp = 2.0, GWPR.point.dataset, coords)
  class(ug_m2_troposphere_no2.idw.2.0)                        # check the inverse distance object
  cuts = c(-0.1, -0.06, -0.03, 0, 0.03, 0.06, 0.09, 0.12)
  
  spplot(ug_m2_troposphere_no2.idw.2.0["var1.pred"],main = "ug_m2_troposphere_no2 Inverse Distance (p=2.0)",
         key.space = "right", cuts = cuts, xlab = "Degree", ylab = "Degree")
  
  
  # Fine-tuning the interpolation
  # the following procedure are learned from https://mgimond.github.io/Spatial/interpolation-in-r.html
  IDW.out <- vector(length = length(GWPR.point.dataset))
  for (i in 1:length(GWPR.point.dataset)) {
    IDW.out[i] <- idw(ug_m2_troposphere_no2~1, GWPR.point.dataset[-i,],
                      GWPR.point.dataset[i,], idp = 2.0)$var1.pred
  }
  
  # Plot the differences
  OP <- par(pty="s", mar=c(4,3,0,0))
  plot(IDW.out ~ GWPR.point.dataset$ug_m2_troposphere_no2, asp=1, xlab="Observed", ylab="Predicted", pch=16,
       col=rgb(0,0,0,0.5))
  abline(lm(IDW.out ~ GWPR.point.dataset$ug_m2_troposphere_no2), col="red", lw=2,lty=2)
  abline(0,1)
  par(OP)
  cor.test(IDW.out, GWPR.point.dataset$ug_m2_troposphere_no2)
  # RMSE
  sqrt( sum((IDW.out - GWPR.point.dataset$ug_m2_troposphere_no2)^2) / length(GWPR.point.dataset))
  # R2
  1 - sum( (IDW.out - GWPR.point.dataset$ug_m2_troposphere_no2)^2) / 
    sum( ( GWPR.point.dataset$ug_m2_troposphere_no2)^2 )
}

## 
#-------------------------------------------kriging----------------------------------
## https://swilke-geoscience.net/post/2020-09-10-kriging_with_r/kriging/
### ug_m2_troposphere_no2
ug_m2_troposphere_no2_emp_OK <- gstat::variogram(
  ug_m2_troposphere_no2 ~ 1, GWPR.point.dataset
)
plot(ug_m2_troposphere_no2_emp_OK)
# Compute the variogram model by passing the nugget, sill and range values
# to fit.variogram() via the vgm() function.
dat.fit  <- fit.variogram(ug_m2_troposphere_no2_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
#v_mod_ok = autofitVariogram(ug_m2_troposphere_no2 ~ 1, GWPR.point.dataset)
#plot(v_mod_ok)
out <- 
  krige.cv(ug_m2_troposphere_no2~1, GWPR.point.dataset, GWPR.point.dataset@data, 
        model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

ug_m2_troposphere_no2.kriged <- 
  krige(ug_m2_troposphere_no2~1, GWPR.point.dataset, coords, 
        model = dat.fit)
ug_m2_troposphere_no2.kriged@data$Ttest <- ug_m2_troposphere_no2.kriged@data$var1.pred / 
  ug_m2_troposphere_no2.kriged@data$var1.var
ug_m2_troposphere_no2.kriged.raster <- as(ug_m2_troposphere_no2.kriged, 'SpatialPixelsDataFrame')
ug_m2_troposphere_no2.kriged.raster <- as(ug_m2_troposphere_no2.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, ug_m2_troposphere_no2.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$ug_m2_troposphere_no2)^2) / 
  sum( ( GWPR.point.dataset$ug_m2_troposphere_no2 - mean(GWPR.point.dataset$ug_m2_troposphere_no2) )^2 )
cv.line <- c("ug_m2_troposphere_no2", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
ug_m2_troposphere_no2.kriged.raster@data <- ug_m2_troposphere_no2.kriged.raster@data %>% dplyr::select(var1.pred)
CO_ug_m2_troposphere_no2.kriged.raster <- raster(ug_m2_troposphere_no2.kriged.raster)

### ter_pressure
ter_pressure_emp_OK <- gstat::variogram(
  ter_pressure ~ 1, GWPR.point.dataset
)
plot(ter_pressure_emp_OK)
dat.fit  <- fit.variogram(ter_pressure_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(ter_pressure~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

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
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$ter_pressure)^2) / 
  sum( ( GWPR.point.dataset$ter_pressure - mean(GWPR.point.dataset$ter_pressure) )^2 )
cv.line <- c("ter_pressure",
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
ter_pressure.kriged.raster@data <- ter_pressure.kriged.raster@data %>% dplyr::select(var1.pred)
CO_ter_pressure.kriged.raster <- raster(ter_pressure.kriged.raster)

### temp
temp_emp_OK <- gstat::variogram(
  temp ~ 1, GWPR.point.dataset
)
plot(temp_emp_OK)
dat.fit  <- fit.variogram(temp_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(temp~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

temp.kriged <- 
  krige(temp~1, GWPR.point.dataset, coords, 
        model = dat.fit)
temp.kriged@data$Ttest <- temp.kriged@data$var1.pred / 
  temp.kriged@data$var1.var
temp.kriged.raster <- as(temp.kriged, 'SpatialPixelsDataFrame')
temp.kriged.raster <- as(temp.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, temp.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$temp)^2) / 
  sum( ( GWPR.point.dataset$temp - mean(GWPR.point.dataset$temp) )^2 )
cv.line <- c("temp",
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
temp.kriged.raster@data <- temp.kriged.raster@data %>% dplyr::select(var1.pred)
CO_temp.kriged.raster <- raster(temp.kriged.raster)
rm(temp_emp_OK, dat.fit, temp.kriged, 
   temp.kriged.raster)

### ndvi
ndvi_emp_OK <- gstat::variogram(
  ndvi ~ 1, GWPR.point.dataset
)
plot(ndvi_emp_OK)
dat.fit  <- fit.variogram(ndvi_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(ndvi~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

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
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$ndvi)^2) / 
  sum( ( GWPR.point.dataset$ndvi - mean(GWPR.point.dataset$ndvi) )^2 )
cv.line <- c("ndvi", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
ndvi.kriged.raster@data <- ndvi.kriged.raster@data %>% dplyr::select(var1.pred)
CO_ndvi.kriged.raster <- raster(ndvi.kriged.raster)
rm(ndvi_emp_OK, dat.fit, ndvi.kriged, 
   ndvi.kriged.raster)

### precipitation
precipitation_emp_OK <- gstat::variogram(
  precipitation ~ 1, GWPR.point.dataset
)
plot(precipitation_emp_OK)
dat.fit  <- fit.variogram(precipitation_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(precipitation~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

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
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$precipitation)^2) / 
  sum( ( GWPR.point.dataset$precipitation - mean(GWPR.point.dataset$precipitation) )^2 )
cv.line <- c("precipitation", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
precipitation.kriged.raster@data <- precipitation.kriged.raster@data %>% dplyr::select(var1.pred)
CO_precipitation.kriged.raster <- raster(precipitation.kriged.raster)
rm(precipitation_emp_OK, dat.fit, precipitation.kriged, 
   precipitation.kriged.raster)

run <- F
if(run){
  ### NTL low accuracy of interpolation. R2 is rough 0.20 
  NTL_emp_OK <- gstat::variogram(
    NTL ~ 1, GWPR.point.dataset
  )
  plot(NTL_emp_OK)
  dat.fit  <- fit.variogram(NTL_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                            vgm(model = "Sph"))
  NTL.kriged <- 
    krige(NTL~1, GWPR.point.dataset, coords, 
          model = dat.fit)
  out <- 
    krige.cv(NTL~1, GWPR.point.dataset, GWPR.point.dataset@data, 
             model = dat.fit)
  N <- nrow(GWPR.point.dataset@data)
  R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
  MRSE <- sqrt(mean(out$residual^2))
  MAE <- mean(abs(out$residual))
  correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
  out$pred <- out$observed - out$residual
  coef.reg <- coef(lm(observed ~ pred, data = out))
  
  cv.line <- c("NTL", 
               N, R2_interpolation, MRSE, MAE, 
               correlation_observed_predicted, coef.reg, onpoint.R2)
  kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
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
}

### PBLH
PBLH_emp_OK <- gstat::variogram(
  PBLH ~ 1, GWPR.point.dataset
)
plot(PBLH_emp_OK)
dat.fit  <- fit.variogram(PBLH_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(model = "Sph"))
out <- 
  krige.cv(PBLH~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))
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
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$PBLH)^2) / 
  sum( ( GWPR.point.dataset$PBLH - mean(GWPR.point.dataset$PBLH) )^2 )
cv.line <- c("PBLH",
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
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
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

out <- 
  krige.cv(Y2016~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
Y2016.kriged@data$Ttest <- Y2016.kriged@data$var1.pred / 
  Y2016.kriged@data$var1.var
Y2016.kriged.raster <- as(Y2016.kriged, 'SpatialPixelsDataFrame')
Y2016.kriged.raster <- as(Y2016.kriged.raster, "SpatialGridDataFrame")
predict.value <- over(GWPR.point.dataset, Y2016.kriged.raster)
GWPR.point.dataset$predict.value <- predict.value$var1.pred
# R2
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2016)^2) / 
  sum( ( GWPR.point.dataset$Y2016 - mean(GWPR.point.dataset$Y2016) )^2 )
cv.line <- c("Y2016", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
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
out <- 
  krige.cv(Y2017~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

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
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2017)^2) / 
  sum( ( GWPR.point.dataset$Y2017 - mean(GWPR.point.dataset$Y2017) )^2 )
cv.line <- c("Y2017",              
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
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
out <- 
  krige.cv(Y2018~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

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
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2018)^2) / 
  sum( ( GWPR.point.dataset$Y2018 - mean(GWPR.point.dataset$Y2018) )^2 )
cv.line <- c("Y2018", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
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
out <- 
  krige.cv(Y2019~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

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
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2019)^2) / 
  sum( ( GWPR.point.dataset$Y2019 - mean(GWPR.point.dataset$Y2019) )^2 )
cv.line <- c("Y2019", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
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
out <- 
  krige.cv(Y2020~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

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
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2020)^2) / 
  sum( ( GWPR.point.dataset$Y2020 - mean(GWPR.point.dataset$Y2020) )^2 )
cv.line <- c("Y2020", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
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
out <- 
  krige.cv(Y2021~1, GWPR.point.dataset, GWPR.point.dataset@data, 
           model = dat.fit)
N <- nrow(GWPR.point.dataset@data)
R2_interpolation <- 1 - sum( out$residual^2) / sum( (out$observed - mean(out$observed) )^2 )
MRSE <- sqrt(mean(out$residual^2))
MAE <- mean(abs(out$residual))
correlation_observed_predicted <- cor(out$observed, out$observed - out$residual)
out$pred <- out$observed - out$residual
coef.reg <- coef(lm(observed ~ pred, data = out))

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
onpoint.R2 <- 1 - sum( (GWPR.point.dataset$predict.value - GWPR.point.dataset$Y2021)^2) / 
  sum( ( GWPR.point.dataset$Y2021 - mean(GWPR.point.dataset$Y2021) )^2 )
cv.line <- c("Y2021", 
             N, R2_interpolation, MRSE, MAE, 
             correlation_observed_predicted, coef.reg, onpoint.R2)
print(cv.line)
kriging.cv.result.dataset <- rbind(kriging.cv.result.dataset, cv.line)
Y2021.kriged.raster@data <- Y2021.kriged.raster@data %>% dplyr::select(var1.pred)
CO_Y2021.kriged.raster <- raster(Y2021.kriged.raster)
rm(Y2021_emp_OK, dat.fit, Y2021.kriged, 
   Y2021.kriged.raster)


colnames(kriging.cv.result.dataset) <- c("Variable", "N", "R2", "MRSE", "MAE", 
                                         "r", "Intercept", "Slope", "onpoint.R2")

save(CO_ug_m2_troposphere_no2.kriged.raster,
     CO_ndvi.kriged.raster, CO_temp.kriged.raster, 
     CO_PBLH.kriged.raster, CO_precipitation.kriged.raster, 
     CO_ter_pressure.kriged.raster, CO_Y2016.kriged.raster, CO_Y2017.kriged.raster,
     CO_Y2018.kriged.raster, CO_Y2019.kriged.raster, CO_Y2020.kriged.raster,
     CO_Y2021.kriged.raster, file = "05_CoefficientRaster/COEF_raster.RData")

save(kriging.cv.result.dataset, file = "04_Results/krigingCVResult.RData")
write.csv(kriging.cv.result.dataset, file = "08_Tables/krigingCVResult.csv")
