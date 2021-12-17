# Author: M.L.

# end

library(gstat)
library(sp) 
library(raster)


# make the result into point
GWPR.point.dataset <- GWPR.FEM.CV.F.result$SDF

# set the basic raster
h_rst <- raster(nrows = 720, ncols = 1440,
                xmn = -180,ymx = 90,
                xmx = 180,ymn = -90, crs = CRS(proj))

# functionalized preprocess
# note: we have too many layer therefore we have to functionalized
P = GWPR.point.dataset # P is the points
keyParameter.fomula = mg_m2_troposphere_no2 ~ 1
idp.power = 2

# interpolation
fit_IDW <- gstat::gstat(
  formula = keyParameter.fomula,
  data = GWPR.point.dataset,
  set = list(idp = idp.power)
)
interp_IDW <- interpolate(h_rst,fit_IDW)

plot(interp_IDW) #result

# Fine-tuning the interpolation
# the following procedure are learned from https://mgimond.github.io/Spatial/interpolation-in-r.html
IDW.out <- vector(length = length(P))
for (i in 1:length(P)) {
  IDW.out[i] <- idw(keyParameter.fomula, P[-i,], P[i,], idp=idp.power)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ P$mg_m2_troposphere_no2, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ P$mg_m2_troposphere_no2), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
# RMSE
sqrt( sum((IDW.out - P$mg_m2_troposphere_no2)^2) / length(P))
# R2
1 - sum( (IDW.out - P$mg_m2_troposphere_no2)^2) / sum( ( P$mg_m2_troposphere_no2 - mean(P$mg_m2_troposphere_no2))^2 )


## 
#-------------------------------------------kriging----------------------------------
## https://swilke-geoscience.net/post/2020-09-10-kriging_with_r/kriging/
v_emp_OK <- gstat::variogram(
  mg_m2_troposphere_no2 ~ 1, P
)

plot(v_emp_OK)
v_mod_OK <- automap::autofitVariogram(mg_m2_troposphere_no2 ~ 1, P)$var_model
plot(automap::autofitVariogram(mg_m2_troposphere_no2 ~ 1, P))

gridRaster <- as(h_rst, 'SpatialGridDataFrame')

OK <- krige(
  mg_m2_troposphere_no2 ~ 1,  
  P,
  gridRaster,
  model = v_mod_OK  
)




P$X <- coordinates(P)[,1]
P$Y <- coordinates(P)[,2]
# Define the 1st order polynomial equation
f.1 <- as.formula(mg_m2_troposphere_no2 ~ X + Y) 

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)
dat.krg <- krige( f.1, P, grd, dat.fit)

# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
r.m <- mask(r, W)
