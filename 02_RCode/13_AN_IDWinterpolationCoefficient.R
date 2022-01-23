# Author: M.L.

# input: GWPR_FEM_CV_F_result.Rdata
# GWPR_FEM_CV_F_result.Rdata: "GWPR.FEM.CV.F.result" GWPR result with 2.25 fixed distance bandwidth 
#                                                    based on the FEM

# output: COEF_raster.RData
# COEF_raster.RData: "CO_ug_m2_troposphere_no2.idw.raster" interpolation of troposheric no2 coefficient
#                                                             based on IDW
# COEF_raster.RData: "CO_ndvi.idw.raster" interpolation of ndvi coefficient (IDW)
# COEF_raster.RData: "CO_temp.idw.raster" interpolation of night temperature coefficient (IDW)
# COEF_raster.RData: "CO_PBLH.idw.raster" interpolation of PBLH coefficient (IDW)
# COEF_raster.RData: "CO_precipitation.idw.raster" interpolation of precipitation coefficient (IDW)
# COEF_raster.RData: "CO_ter_pressure.idw.raster" interpolation of air pressure coefficient (IDW)
# COEF_raster.RData: "CO_Y2016.idw.raster" interpolation of 2016 dummy variable coefficient (IDW)
# COEF_raster.RData: "CO_Y2017.idw.raster" interpolation of 2017 dummy variable coefficient (IDW)
# COEF_raster.RData: "CO_Y2018.idw.raster" interpolation of 2018 dummy variable coefficient (IDW)
# COEF_raster.RData: "CO_Y2019.idw.raster" interpolation of 2019 dummy variable coefficient (IDW)
# COEF_raster.RData: "CO_Y2020.idw.raster" interpolation of 2020 dummy variable coefficient (IDW)
# COEF_raster.RData: "CO_Y2021.idw.raster" interpolation of 2021 dummy variable coefficient (IDW)

# output: idwCVResult.RData
# idwCVResult.RData: "N" number of observation
# idwCVResult.RData: "R2" R2
# idwCVResult.RData: "RMSE" root mean square error
# idwCVResult.RData: "MAE" mean absolute error
# idwCVResult.RData: "r" correlation coefficient
# idwCVResult.RData: "Intercept" intercept of regression
# idwCVResult.RData: "Slope" slope of regression

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

coords <- addcoord(nx,xmin,xsize,ny,ymin,ysize)

# the tutorial of this interpolation
# https://github.com/GeostatsGuy/geostatsr/blob/master/kriging_demo_Rnotebook.ipynb
GWPR.point.dataset <- spTransform(GWPR.point.dataset, CRS(proj))
coords@proj4string <- CRS(proj)
summary(coords) 
class(coords)


# leave-one-out cv
formula <- no2_measured_ug.m3 ~ ug_m2_troposphere_no2 + 
  ter_pressure + temp + ndvi + precipitation + PBLH + 
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021

IDW.cv.result.dataset <- data.frame(Doubles = double(), Ints = integer(),
                                        Factors = factor(), Logicals = logical(),
                                        Characters = character(), stringsAsFactors = FALSE)

indep.variables <- all.vars(formula)[2:length(all.vars(formula))]
for (indep.variable in indep.variables){
  aim.formula <- paste0(indep.variable, "~ 1") %>% as.formula()
  line <- IDW.cv(aim.formula, GWPR.point.dataset)
  IDW.cv.result.dataset <- rbind(IDW.cv.result.dataset, line)
}
colnames(IDW.cv.result.dataset) <- c("Variable", "N", "R2", "RMSE", "MAE", "r",
                                     "Intercept","Slope") 
save(IDW.cv.result.dataset, file = "04_Results/IDWCVResult.RData")

# idw interpolation coefficient
CO_ug_m2_troposphere_no2.idw.raster <- idw.raster(ug_m2_troposphere_no2 ~ 1, GWPR.point.dataset, coords)
CO_ndvi.idw.raster <- idw.raster(ndvi ~ 1, GWPR.point.dataset, coords)
CO_temp.idw.raster  <- idw.raster(temp ~ 1, GWPR.point.dataset, coords)
CO_PBLH.idw.raster <- idw.raster(PBLH ~ 1, GWPR.point.dataset, coords)
CO_precipitation.idw.raster <- idw.raster(precipitation ~ 1, GWPR.point.dataset, coords)
CO_ter_pressure.idw.raster <- idw.raster(ter_pressure ~ 1, GWPR.point.dataset, coords)
CO_Y2016.idw.raster <- idw.raster(Y2016 ~ 1, GWPR.point.dataset, coords)
CO_Y2017.idw.raster <- idw.raster(Y2017 ~ 1, GWPR.point.dataset, coords)
CO_Y2018.idw.raster <- idw.raster(Y2018 ~ 1, GWPR.point.dataset, coords)
CO_Y2019.idw.raster <- idw.raster(Y2019 ~ 1, GWPR.point.dataset, coords)
CO_Y2020.idw.raster <- idw.raster(Y2020 ~ 1, GWPR.point.dataset, coords)
CO_Y2021.idw.raster <- idw.raster(Y2021 ~ 1, GWPR.point.dataset, coords)

save(CO_ug_m2_troposphere_no2.idw.raster,
     CO_ndvi.idw.raster, CO_temp.idw.raster, 
     CO_PBLH.idw.raster, CO_precipitation.idw.raster, 
     CO_ter_pressure.idw.raster, CO_Y2016.idw.raster, CO_Y2017.idw.raster,
     CO_Y2018.idw.raster, CO_Y2019.idw.raster, CO_Y2020.idw.raster,
     CO_Y2021.idw.raster, file = "05_CoefficientRaster/idw_COEF_raster.RData")

