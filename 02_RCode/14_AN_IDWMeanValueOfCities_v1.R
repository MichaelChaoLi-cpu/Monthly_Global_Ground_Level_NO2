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

# output: idw.MEAN_raster.RData
# idw.MEAN_raster.RData: "MEAN_ug_m2_troposphere_no2.idw.raster" interpolation of troposheric no2 mean value
#                                                               in the usedDataset (IDW)
# idw.MEAN_raster.RData: "MEAN_ndvi.idw.raster" interpolation of ndvi mean value (IDW)
# idw.MEAN_raster.RData: "MEAN_temp.idw.raster" interpolation of temperature mean value (IDW)
# idw.MEAN_raster.RData: "MEAN_NTL.idw.raster" interpolation of night time light mean value (IDW)
# idw.MEAN_raster.RData: "MEAN_PBLH.idw.raster" interpolation of PBLH mean value (IDW)
# idw.MEAN_raster.RData: "MEAN_precipitation.idw.raster" interpolation of precipitation mean value (IDW)
# idw.MEAN_raster.RData: "MEAN_ter_pressure.idw.raster" interpolation of air pressure mean value (IDW)

# output: idwMeanResult.RData
# idwMeanResult.RData: "N" number of observation
# idwMeanResult.RData: "R2" R2
# idwMeanResult.RData: "RMSE" root mean square error
# idwMeanResult.RData: "MAE" mean absolute error
# idwMeanResult.RData: "r" correlation coefficient
# idwMeanResult.RData: "Intercept" intercept of regression
# idwMeanResult.RData: "Slope" slope of regression

# end

library(tidyverse)
library(gstat)
library(sp) 
library(raster)
library(dplyr)
library(tmap)

IDW.cv <- function(formula.use, GWPR.point.dataset){
  IDW.out <- vector(length = length(GWPR.point.dataset))
  for (i in 1:length(GWPR.point.dataset)) {
    IDW.out[i] <- idw(formula.use, GWPR.point.dataset[-i,],
                      GWPR.point.dataset[i,], idp = 2.0, nmax = 7)$var1.pred
  }
  IDW.out <- IDW.out %>% as.data.frame()
  aim.data <- GWPR.point.dataset@data %>% dplyr::select(all.vars(formula.use))
  IDW.out <- cbind(IDW.out, aim.data)
  colnames(IDW.out) <- c("cv.Predict", "observed")
  IDW.out$residual <- IDW.out$cv.Predict - IDW.out$observed
  cor.coef <- cor.test(IDW.out$cv.Predict, IDW.out$observed)
  cor.coef <- cor.coef$estimate
  mae <- mean(abs(IDW.out$residual))
  rmse <- sqrt(mean(IDW.out$residual ^ 2) )
  r2 <- 1 - (sum(IDW.out$residual ^ 2) / sum(  (IDW.out$observed - mean(IDW.out$observed))^2 ))
  N <- length(GWPR.point.dataset)
  test.reg <- lm(observed ~ cv.Predict, IDW.out)
  coef <- coef(test.reg)
  names <- all.vars(formula.use)
  result <- c(names, N, r2, rmse, mae, cor.coef, coef)
  return(result)
}

idw.raster <- function(formula.use, GWPR.point.dataset = GWPR.point.dataset,
                       coords = coords){
  idw.result <-
    gstat::idw(formula.use, GWPR.point.dataset,
               newdata=coords, nmax = 7, idp = 2.0)
  idw.result.raster <- as(idw.result, 'SpatialPixelsDataFrame')
  idw.result.raster <- as(idw.result.raster, "SpatialGridDataFrame")
  idw.result.raster@data <- idw.result.raster@data %>% dplyr::select(var1.pred)
  raster::raster(idw.result.raster) ->
    out.raster
  return(out.raster)
}

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

# leave-one-out cv
formula <- no2_measured_ug.m3 ~ ug_m2_troposphere_no2 + 
  ter_pressure + temp + ndvi + precipitation + PBLH + 
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021

IDW.cv.mean.dataset <- data.frame(Doubles = double(), Ints = integer(),
                                    Factors = factor(), Logicals = logical(),
                                    Characters = character(), stringsAsFactors = FALSE)
line <- IDW.cv(no2_measured_ug.m3~1, meanValueSpatialDataFrame)
IDW.cv.mean.dataset <- rbind(IDW.cv.mean.dataset, line)

indep.variables <- all.vars(formula)[2:length(all.vars(formula))]
for (indep.variable in indep.variables){
  aim.formula <- paste0(indep.variable, "~ 1") %>% as.formula()
  line <- IDW.cv(aim.formula, meanValueSpatialDataFrame)
  IDW.cv.mean.dataset <- rbind(IDW.cv.mean.dataset, line)
}
colnames(IDW.cv.mean.dataset) <- c("Variable", "N", "R2", "RMSE", "MAE", "r",
                                     "Intercept","Slope") 
save(IDW.cv.mean.dataset, file = "04_Results/idwMeanResult.RData")

# idw interpolation mean
MEAN_no2_measured_ug.m3.idw.raster <- idw.raster(no2_measured_ug.m3 ~ 1, meanValueSpatialDataFrame, coords)
MEAN_ug_m2_troposphere_no2.idw.raster <- idw.raster(ug_m2_troposphere_no2 ~ 1, meanValueSpatialDataFrame, coords)
MEAN_ndvi.idw.raster <- idw.raster(ndvi ~ 1, meanValueSpatialDataFrame, coords)
MEAN_temp.idw.raster  <- idw.raster(temp ~ 1, meanValueSpatialDataFrame, coords)
MEAN_PBLH.idw.raster <- idw.raster(PBLH ~ 1, meanValueSpatialDataFrame, coords)
MEAN_precipitation.idw.raster <- idw.raster(precipitation ~ 1, meanValueSpatialDataFrame, coords)
MEAN_ter_pressure.idw.raster <- idw.raster(ter_pressure ~ 1, meanValueSpatialDataFrame, coords)
MEAN_Y2016.idw.raster <- idw.raster(Y2016 ~ 1, meanValueSpatialDataFrame, coords)
MEAN_Y2017.idw.raster <- idw.raster(Y2017 ~ 1, meanValueSpatialDataFrame, coords)
MEAN_Y2018.idw.raster <- idw.raster(Y2018 ~ 1, meanValueSpatialDataFrame, coords)
MEAN_Y2019.idw.raster <- idw.raster(Y2019 ~ 1, meanValueSpatialDataFrame, coords)
MEAN_Y2020.idw.raster <- idw.raster(Y2020 ~ 1, meanValueSpatialDataFrame, coords)
MEAN_Y2021.idw.raster <- idw.raster(Y2021 ~ 1, meanValueSpatialDataFrame, coords)


save(MEAN_no2_measured_ug.m3.idw.raster,
     MEAN_ug_m2_troposphere_no2.idw.raster,
     MEAN_ndvi.idw.raster, MEAN_temp.idw.raster, 
     MEAN_PBLH.idw.raster, MEAN_precipitation.idw.raster, 
     MEAN_ter_pressure.idw.raster, MEAN_Y2016.idw.raster,
     MEAN_Y2017.idw.raster, MEAN_Y2018.idw.raster,
     MEAN_Y2019.idw.raster, MEAN_Y2020.idw.raster,
     MEAN_Y2021.idw.raster, file = "05_CoefficientRaster/idw.MEAN_raster.RData")

save(IDW.cv.mean.dataset, file = "04_Results/idwMeanResult.RData")

