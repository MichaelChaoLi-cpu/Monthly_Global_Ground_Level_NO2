# Author: M.L.

# end

library(gstat)
library(sp) 
library(raster)
library(dplyr)

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
cuts = c(-1., -0.5, 0, .005,.007,.009,.011,.013,.015,.017,.019,.021,.023)

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
# RMSE
sqrt( sum((IDW.out - GWPR.point.dataset$mg_m2_troposphere_no2)^2) / length(GWPR.point.dataset))
# R2
1 - sum( (IDW.out - GWPR.point.dataset$mg_m2_troposphere_no2)^2) / 
  sum( ( GWPR.point.dataset$mg_m2_troposphere_no2 - mean(GWPR.point.dataset$mg_m2_troposphere_no2))^2 )


## 
#-------------------------------------------kriging----------------------------------
## https://swilke-geoscience.net/post/2020-09-10-kriging_with_r/kriging/
mg_m2_troposphere_no2_emp_OK <- gstat::variogram(
  mg_m2_troposphere_no2 ~ 1, GWPR.point.dataset, cutoff = 30, width = 0.25
)

# Compute the variogram model by passing the nugget, sill and range values
# to fit.variogram() via the vgm() function.
dat.fit  <- fit.variogram(mg_m2_troposphere_no2_emp_OK, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill = 2.25, model="Sph", range = 10, nugget = 0))

plot(mg_m2_troposphere_no2_emp_OK$dist, mg_m2_troposphere_no2_emp_OK$gamma, 
     main="mg_m2_troposphere_no2_emp_OK Variogram",xlab="  Lag Distance (Degree) ",
     ylab=" Semivariogram ", pch=16, col = "red", ylim=c(0,max(mg_m2_troposphere_no2_emp_OK$gamma)*1.1))

mg_m2_troposphere_no2.kriged <- 
  krige(mg_m2_troposphere_no2~1, GWPR.point.dataset, coords, 
        model = dat.fit)
