# Author: M.L.

# end

library(tmap)
library(sp)
library(ggplot2)
library(grid)
library(dplyr)
library(MASS)
library(rgdal)
library(rgeos)
library("rnaturalearth")
library(ggpubr)
library(gridExtra)
library("viridisLite")
library("viridis") 

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

set.seed(1)
dat <- data.frame(
  x = c(
    rnorm(1e4, mean = 0, sd = 0.1),
    rnorm(1e3, mean = 0, sd = 0.1)
  ),
  y = c(
    rnorm(1e4, mean = 0, sd = 0.1),
    rnorm(1e3, mean = 0.1, sd = 0.2)
  )
)

world <- ne_countries(scale = "medium", returnclass = "sf")
plot(world)

load("03_Rawdata/cityLocationSpatialPoint.Rdata")
proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
cityLocationSpatialPoint@data$code <- 1
load("03_Rawdata/usedDataset.RData")

# tm set
title_size = .0001
legend_title_size = 1
margin = 0
brk_ols_res = c(-2, -0.5, -1.0, -0.5, 0, 0.5, 1, 0.5, 2)
tmap_mode('plot')

(a_city_location <- 
  tm_shape(world) +
  tm_polygons(col = 'level', lwd = 0.01, alpha = .8, pal = "grey94", legend.show = F) +
  tm_shape(cityLocationSpatialPoint) +
  tm_dots(col = 'code', pal = "red2", border.alpha = 0, shape = 17,
          legend.show = F, size = 0.1) +
  tm_grid(projection = proj, alpha = .25) + 
  tm_layout(
    inner.margins = c(margin, margin, margin, margin),
    title.size = title_size, 
    legend.position = c("left", "top"),
    legend.title.size = legend_title_size,
    legend.text.size = legend_title_size * 0.75
  ) + 
  tm_scale_bar())

a_city_location %>%
  tmap_save(filename = "07_Figure/01_point_location.jpg" ,width = 210, height = 120, units = 'mm', dpi = 1000)

#-------------descriptive statistics--------------
Mean <- round(mean(usedDataset$no2_measured_mg.m3), 2)
SD <- round(sd(usedDataset$no2_measured_mg.m3), 2)
N = nrow(usedDataset)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("a",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(a <- ggplot(usedDataset) +
  aes(x = no2_measured_mg.m3) +
  xlim(0, 1.5) +
  geom_histogram(colour = "black", fill = "white") +
  xlab("Ground-level NO2 (mg/m3)") + 
  ylab("Frequency") +
  annotation_custom(grob) +
  annotation_custom(grob_add))

Mean <- round(mean(usedDataset$mg_m2_troposphere_no2), 2)
SD <- round(sd(usedDataset$mg_m2_troposphere_no2), 2)
N = nrow(usedDataset)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("b",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(b <- ggplot(usedDataset) +
  aes(x = mg_m2_troposphere_no2) +
  geom_histogram(colour = "black", fill = "white") +
  xlim(0, 20) +
  xlab("OMI TrCA NO2 (mg/m2)") + 
  ylab("Frequency") +
  annotation_custom(grob) +
  annotation_custom(grob_add))

Mean <- round(mean(usedDataset$ter_pressure), 2)
SD <- round(sd(usedDataset$ter_pressure), 2)
N = nrow(usedDataset)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.61,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("c",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(c <- ggplot(usedDataset) +
    aes(x = ter_pressure) +
    geom_histogram(colour = "black", fill = "white") +
    #xlim(0, 20) +
    xlab("Terrain Atmospheric Pressure (hPa)") + 
    ylab("Frequency") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

Mean <- round(mean(usedDataset$temp), 2)
SD <- round(sd(usedDataset$temp), 2)
N = nrow(usedDataset)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("d",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(d <- ggplot(usedDataset) +
    aes(x = temp) +
    geom_histogram(colour = "black", fill = "white") +
    #xlim(0, 20) +
    xlab("Temperature (Celsius Degree)") + 
    ylab("Frequency") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

Mean <- round(mean(usedDataset$ndvi), 2)
SD <- round(sd(usedDataset$ndvi), 2)
N = nrow(usedDataset)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("e",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(e <- ggplot(usedDataset) +
    aes(x = ndvi) +
    geom_histogram(colour = "black", fill = "white") +
    #xlim(0, 20) +
    xlab("NDVI") + 
    ylab("Frequency") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

Mean <- round(mean(usedDataset$precipitation), 2)
SD <- round(sd(usedDataset$precipitation), 2)
N = nrow(usedDataset)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("f",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(f <- ggplot(usedDataset) +
    aes(x = precipitation) +
    geom_histogram(colour = "black", fill = "white") +
    xlim(0, 0.8) +
    xlab("Precipitation kg/(m2*h)") + 
    ylab("Frequency") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

Mean <- round(mean(usedDataset$PBLH), 2)
SD <- round(sd(usedDataset$PBLH), 2)
N = nrow(usedDataset)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("g",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(g <- ggplot(usedDataset) +
    aes(x = PBLH) +
    geom_histogram(colour = "black", fill = "white") +
    xlim(0, 2500) +
    xlab("PBLH (m)") + 
    ylab("Frequency") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

jpeg(file="07_Figure\\descriptive_stat.jpeg", width = 297, height = 315, units = "mm", quality = 300, res = 300)
grid.arrange(a, b, c,
             d, e, f,
             g,
             nrow = 3)
dev.off()

formula <- no2_measured_mg.m3 ~ mg_m2_troposphere_no2 + 
  ter_pressure + 
  temp +
  ndvi + precipitation + PBLH + 
  #humidity + UVAerosolIndex + ozone + speedwind + NTL +
  #cloudfraction + cloudpressure + # add this two variables effect are limited, only increase 0.2% R2
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021

#-------------trends of XY------------------------
(plot1 <- ggscatter(usedDataset %>% data.frame(), x = "mg_m2_troposphere_no2", y = "no2_measured_mg.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "OMI TrCA NO2", ylab = "Ground-level NO2",           
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21
) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

(plot2 <- ggscatter(usedDataset, x = "ter_pressure", y = "no2_measured_mg.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "Terrian Atmospheric Pressure", ylab = "Ground-level NO2",
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

(plot3 <- ggscatter(usedDataset, x = "temp", y = "no2_measured_mg.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "Temperature", ylab = "Ground-level NO2",
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

(plot4 <- ggscatter(usedDataset, x = "ndvi", y = "no2_measured_mg.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "NDVI", ylab = "Ground-level NO2",
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

(plot5 <- ggscatter(usedDataset, x = "precipitation", y = "no2_measured_mg.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "Precipitation", ylab = "Ground-level NO2",
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

(plot6 <- ggscatter(usedDataset, x = "PBLH", y = "no2_measured_mg.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "PBLH", ylab = "Ground-level NO2",
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

jpeg(file="07_Figure/cor_line1.jpeg", width = 297, height = 210, units = "mm", quality = 300, res = 300)
grid.arrange(plot1, plot2, plot3, 
             plot4, plot5, plot6,
             nrow = 2)
dev.off()
#-------------trends of XY------------------------

load("04_Results/dataPrediction.Rdata")
#---------------------GWPR-------------------------------------
data.predict$Density <- get_density(data.predict$predictNo2, data.predict$no2_measured_mg.m3.ori, n = 100)
reg <- lm(predictNo2 ~ no2_measured_mg.m3.ori, data = data.predict )
coeff = coefficients(reg)
eq = paste0("y = ", round(coeff[2],3), "x + ", round(coeff[1],3))
grob <- grobTree(textGrob(eq,
                          x = 0.05,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 12)))
corre <- cor(data.predict$predictNo2, data.predict$no2_measured_mg.m3.ori)
corr.text <- paste0("r = ", round(corre,4))
grob.corr <- grobTree(textGrob(corr.text,
                               x = 0.05,  y = 0.85, hjust = 0,
                               gp = gpar(col = "black", fontsize = 12)))
N <- length(data.predict$no2_measured_mg.m3.ori)
N.text <- paste0("N = ", N)
grob.N <- grobTree(textGrob(N.text,
                            x = 0.05,  y = 0.80, hjust = 0,
                            gp = gpar(col = "black", fontsize = 12)))
grob_add <- grobTree(textGrob("GWPR",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
gwpr.cv <- ggplot(data.predict ) +
  geom_point(aes(x = no2_measured_mg.m3.ori, y = predictNo2, color = Density)) +
  scale_color_viridis() + 
  scale_x_continuous(name = "Measured NO2 (mg/m3)", limits = c(0, 1.5)) +
  scale_y_continuous(name = "Predicted NO2 (mg/m3)", limits = c(0, 1.5)) +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype = "dashed", size = 0.5) + 
  geom_abline(intercept = coeff[1], slope = coeff[2], color="blue", 
              size= 0.5) + 
  annotation_custom(grob) + 
  annotation_custom(grob.corr) +
  annotation_custom(grob.N) +
  annotation_custom(grob_add)

jpeg(file="07_Figure/gwpr.cv.jpeg", width = 210, height = 210, units = "mm", quality = 300, res = 300)
gwpr.cv 
dev.off()

load("04_Results/FinalRasterCrossValidation.Rdata")
#---------------------Raster-------------------------------------
testDataset$Density <- get_density(testDataset$predict_no2, testDataset$no2_measured_mg.m3, n = 100)
reg <- lm(predict_no2 ~ no2_measured_mg.m3, data = testDataset )
coeff = coefficients(reg)
eq = paste0("y = ", round(coeff[2],3), "x + ", round(coeff[1],3))
grob <- grobTree(textGrob(eq,
                          x = 0.05,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 12)))
corre <- cor(testDataset$predict_no2, testDataset$no2_measured_mg.m3)
corr.text <- paste0("r = ", round(corre,4))
grob.corr <- grobTree(textGrob(corr.text,
                               x = 0.05,  y = 0.85, hjust = 0,
                               gp = gpar(col = "black", fontsize = 12)))
N <- length(testDataset$no2_measured_mg.m3)
N.text <- paste0("N = ", N)
grob.N <- grobTree(textGrob(N.text,
                            x = 0.05,  y = 0.80, hjust = 0,
                            gp = gpar(col = "black", fontsize = 12)))
grob_add <- grobTree(textGrob("Raster",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
raster.cv <- ggplot(testDataset) +
  geom_point(aes(x = no2_measured_mg.m3, y = predict_no2, color = Density)) +
  scale_color_viridis() + 
  scale_x_continuous(name = "Measured NO2 (mg/m3)", limits = c(0, 1.5)) +
  scale_y_continuous(name = "Predicted NO2 (mg/m3)", limits = c(0, 1.5)) +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype = "dashed", size = 0.5) + 
  geom_abline(intercept = coeff[1], slope = coeff[2], color="blue", 
              size= 0.5) + 
  annotation_custom(grob) + 
  annotation_custom(grob.corr) +
  annotation_custom(grob.N) +
  annotation_custom(grob_add)

jpeg(file="07_Figure/gwpr.Raster.cv.jpeg", width = 210, height = 210, units = "mm", quality = 300, res = 300)
raster.cv
dev.off()
