# Author: M.L.

# output: trendenceMonthGroundLevel.RData
# trendenceMonthGroundLevel.RData: "month.slope" the slope of regression between prediction and month

# output: meanOfRasterNo2.RData
# meanOfRasterNo2.RData: "ave.value" the mean of 77 month predictions 

# end


library(tidyverse)
library(dplyr)
library(sp) 
library(raster)
library(tmap)

# make the raster -180 -60 180 90
nx = 600                                       # number of cells in the x direction
ny = 1440                                     # number of cells in the y direction
xmin = -179.875                                     # x coordinate of lower, left cell center 
ymin = -59.875                                     # y coordinate of lower, left cell center 
xsize = 0.25                                   # extent of cells in x direction
ysize = 0.25                                   # extent of cells in y direction
proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

addcoord <- function(nx,xmin,xsize,ny,ymin,ysize,proj) { # Michael Pyrcz, March, 2018                      
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
  coords.df$id = 1:nrow(coords.df)
  xy <- coords.df[,c(1,2)]
  coords.df <- SpatialPointsDataFrame(coords = xy, data = coords.df %>% dplyr::select("id"),
                                      proj4string = CRS(proj))
  return (coords.df)
  
}

coords <- addcoord(nx,xmin,xsize,ny,ymin,ysize,proj)

predict_raster_folder <- "D:/10_Article/11_PredictRaster/01_Test0104/"
raster.list <- list.files(predict_raster_folder)

month.count <- 1
while (month.count < length(raster.list) + 1) {
  filename = raster.list[month.count]
  test_tiff <- raster::raster(paste0(predict_raster_folder, filename))
  data_ext <- raster::extract(test_tiff, coords)
  coords@data <- cbind(coords@data, data_ext)
  colnames(coords@data)[ncol(coords@data)] <- paste0("Mon_", as.character(month.count))
  month.count <- month.count + 1
}

month.grid.dataset <- coords@data %>% as.data.frame()
month.grid.dataset$NACount <- rowSums(is.na(month.grid.dataset))
month.grid.dataset <- month.grid.dataset %>%
  filter(NACount < 82)
month.grid.dataset <- month.grid.dataset %>%
  dplyr::select(-NACount)

test.coeff.grid <- data.frame(Doubles=double(),
                                    Ints=integer(),
                                    Factors=factor(),
                                    Logicals=logical(),
                                    Characters=character(),
                                    stringsAsFactors=FALSE)
line.num = 1
while (line.num < nrow(month.grid.dataset) + 1){
  test.single.grid <- month.grid.dataset[line.num,]
  test.single.grid <- as.matrix(test.single.grid)
  id <- test.single.grid[1,1]
  test.single.grid <- test.single.grid[,2:83]
  test.single.grid <- test.single.grid %>% as.data.frame()
  colnames(test.single.grid)[1] <- "data" 
  test.single.grid$month <- 1:82
  
  test.reg <- lm(data ~ month, test.single.grid)
  coeff <- summary(test.reg)
  if (nrow(coeff$coefficients) == 2){
    line <- c(id, coeff$coefficients[1,], coeff$coefficients[2,])
  }
  test.coeff.grid <- rbind(test.coeff.grid, line)
  line.num <- line.num + 1
}
colnames(test.coeff.grid) <- c("id", "Intercept", "I.Std.Error","I.t.value","I.Pr",
                               "month.slope", "M.Std.Error","M.t.value","M.Pr")
xy <- coordinates(coords)
id.xy <- cbind(coords@data, xy) 
id.xy <- id.xy %>% dplyr::select(id, X, Y)

test.coeff.grid$slope.sig <- NA
test.coeff.grid <- test.coeff.grid %>% 
  mutate(slope.sig = ifelse((M.Pr < 0.1), T, F))
test.coeff.grid.sp <- left_join(test.coeff.grid, id.xy, by = "id")
xy <- test.coeff.grid.sp[,c("X", "Y")]
test.coeff.grid.sp <- SpatialPointsDataFrame(coords = xy, data = test.coeff.grid.sp,
                                    proj4string = CRS(proj))

test.coeff.grid.raster <- as(test.coeff.grid.sp, 'SpatialPixelsDataFrame')
test.coeff.grid.raster <- as(test.coeff.grid.raster, "SpatialGridDataFrame")
test.coeff.grid.raster.output <- test.coeff.grid.raster
test.coeff.grid.raster.output@data <- test.coeff.grid.raster.output@data %>%
  mutate(month.slope = ifelse(slope.sig == F, NA, month.slope))
test.coeff.grid.raster.output@data <- test.coeff.grid.raster.output@data %>% dplyr::select("month.slope")
test.coeff.grid.raster.output <- raster(test.coeff.grid.raster.output)

save(test.coeff.grid.raster.output, file = "04_Results/trendenceMonthGroundLevel.RData")

month.grid.dataset$ave.value <- rowMeans(month.grid.dataset[,2:83], na.rm = T)
month.grid.dataset.mean <- month.grid.dataset %>%
  dplyr::select("id", "ave.value")

xy <- coordinates(coords)
id.xy <- cbind(coords@data, xy) 
id.xy <- id.xy %>% dplyr::select(id, X, Y)
id.xy <- left_join(id.xy, month.grid.dataset.mean, by = 'id')
test.mean.grid.sp <- SpatialPointsDataFrame(coords = xy, data = id.xy, proj4string = CRS(proj))
test.mean.grid.raster <- as(test.mean.grid.sp, 'SpatialPixelsDataFrame')
test.mean.grid.raster <- as(test.mean.grid.raster, "SpatialGridDataFrame")
test.mean.grid.raster.output <- test.mean.grid.raster
test.mean.grid.raster.output@data <- test.mean.grid.raster.output@data %>% dplyr::select("ave.value")
test.mean.grid.raster.output <- raster(test.mean.grid.raster.output)

save(test.mean.grid.raster.output, file = "04_Results/meanOfRasterNo2.RData")
