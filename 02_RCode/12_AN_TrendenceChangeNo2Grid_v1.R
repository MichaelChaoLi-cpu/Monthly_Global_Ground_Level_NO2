# Author: M.L.

# end


library(tidyverse)
library(dplyr)
library(sp) 
library(raster)

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
  filter(NACount < 77)
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
  test.single.grid <- test.single.grid[,2:78]
  test.single.grid <- test.single.grid %>% as.data.frame()
  colnames(test.single.grid)[1] <- "data" 
  test.single.grid$month <- 1:77
  
  test.reg <- lm(data ~ month, test.single.grid)
  coeff = coefficients(test.reg)
  line <- c(id, coeff)
  test.coeff.grid <- rbind(test.coeff.grid, line)
  line.num <- line.num + 1
}
colnames(test.coeff.grid) <- c("id", "Intercept", "month.slope")