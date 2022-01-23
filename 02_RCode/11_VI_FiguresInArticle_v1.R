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
library(raster)
library(stringr)
library(magick)
library(cowplot)
library(patchwork)
library(plotrix)

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
Mean <- round(mean(usedDataset$no2_measured_ug.m3), 2)
SD <- round(sd(usedDataset$no2_measured_ug.m3), 2)
N = nrow(usedDataset)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("a",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(a <- ggplot(usedDataset) +
  aes(x = no2_measured_ug.m3) +
  xlim(0, 100) +
  geom_histogram(colour = "black", fill = "white") +
  xlab("Ground-level NO2 (ug/m3)") + 
  ylab("Frequency") +
  annotation_custom(grob) +
  annotation_custom(grob_add))

Mean <- round(mean(usedDataset$ug_m2_troposphere_no2), 2)
SD <- round(sd(usedDataset$ug_m2_troposphere_no2), 2)
N = nrow(usedDataset)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("b",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(b <- ggplot(usedDataset) +
  aes(x = ug_m2_troposphere_no2) +
  geom_histogram(colour = "black", fill = "white") +
  xlim(0, 10000) +
  xlab("OMI TrCA NO2 (ug/m2)") + 
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

formula <- no2_measured_ug.m3 ~ ug_m2_troposphere_no2 + 
  ter_pressure + 
  temp +
  ndvi + precipitation + PBLH + 
  #humidity + UVAerosolIndex + ozone + speedwind + NTL +
  #cloudfraction + cloudpressure + # add this two variables effect are limited, only increase 0.2% R2
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021

#-------------trends of XY------------------------
(plot1 <- ggscatter(usedDataset %>% data.frame(), x = "ug_m2_troposphere_no2", y = "no2_measured_ug.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "OMI TrCA NO2", ylab = "Ground-level NO2",           
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21
) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

(plot2 <- ggscatter(usedDataset, x = "ter_pressure", y = "no2_measured_ug.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "Terrain Atmospheric Pressure", ylab = "Ground-level NO2",
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

(plot3 <- ggscatter(usedDataset, x = "temp", y = "no2_measured_ug.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "Temperature", ylab = "Ground-level NO2",
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

(plot4 <- ggscatter(usedDataset, x = "ndvi", y = "no2_measured_ug.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "NDVI", ylab = "Ground-level NO2",
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

(plot5 <- ggscatter(usedDataset, x = "precipitation", y = "no2_measured_ug.m3", size = 1,
                   add = "reg.line", conf.int = TRUE,
                   cor.coef = F, cor.method = "pearson",
                   xlab = "Precipitation", ylab = "Ground-level NO2",
                   add.params = list(color = "blue", fill = "lightskyblue1"),
                   color = "grey76", shape = 21) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01))

(plot6 <- ggscatter(usedDataset, x = "PBLH", y = "no2_measured_ug.m3", size = 1,
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
data.predict$Density <- get_density(data.predict$predictNo2, data.predict$no2_measured_ug.m3.ori, n = 1000)
reg <- lm(predictNo2 ~ no2_measured_ug.m3.ori, data = data.predict )
coeff = coefficients(reg)
eq = paste0("y = ", round(coeff[2],3), "x + ", round(coeff[1],3))
grob <- grobTree(textGrob(eq,
                          x = 0.05,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 12)))
corre <- cor(data.predict$predictNo2, data.predict$no2_measured_ug.m3.ori)
corr.text <- paste0("r = ", round(corre,4))
grob.corr <- grobTree(textGrob(corr.text,
                               x = 0.05,  y = 0.85, hjust = 0,
                               gp = gpar(col = "black", fontsize = 12)))
N <- length(data.predict$no2_measured_ug.m3.ori)
N.text <- paste0("N = ", N)
grob.N <- grobTree(textGrob(N.text,
                            x = 0.05,  y = 0.80, hjust = 0,
                            gp = gpar(col = "black", fontsize = 12)))
grob_add <- grobTree(textGrob("GWPR",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(gwpr.cv <- ggplot(data.predict ) +
  geom_point(aes(x = no2_measured_ug.m3.ori, y = predictNo2, color = Density)) +
  scale_color_viridis() + 
  scale_x_continuous(name = "Measured NO2 (ug/m3)", limits = c(0, 100)) +
  scale_y_continuous(name = "Predicted NO2 (ug/m3)", limits = c(0, 100)) +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype = "dashed", size = 0.5) + 
  geom_abline(intercept = coeff[1], slope = coeff[2], color="blue", 
              size= 0.5) + 
  annotation_custom(grob) + 
  annotation_custom(grob.corr) +
  annotation_custom(grob.N) +
  annotation_custom(grob_add))

jpeg(file="07_Figure/gwpr.cv.jpeg", width = 210, height = 210, units = "mm", quality = 300, res = 300)
gwpr.cv 
dev.off()

load("04_Results/kriged.FinalRasterCrossValidation.Rdata")
#---------------------Raster-------------------------------------
testDataset$Density <- get_density(testDataset$predict_no2, testDataset$no2_measured_ug.m3, n = 1000)
reg <- lm(predict_no2 ~ no2_measured_ug.m3, data = testDataset )
coeff = coefficients(reg)
eq = paste0("y = ", round(coeff[2],3), "x + ", round(coeff[1],3))
grob <- grobTree(textGrob(eq,
                          x = 0.05,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 12)))
corre <- cor(testDataset$predict_no2, testDataset$no2_measured_ug.m3)
corr.text <- paste0("r = ", round(corre,4))
grob.corr <- grobTree(textGrob(corr.text,
                               x = 0.05,  y = 0.85, hjust = 0,
                               gp = gpar(col = "black", fontsize = 12)))
N <- length(testDataset$no2_measured_ug.m3)
N.text <- paste0("N = ", N)
grob.N <- grobTree(textGrob(N.text,
                            x = 0.05,  y = 0.80, hjust = 0,
                            gp = gpar(col = "black", fontsize = 12)))
grob_add <- grobTree(textGrob("Raster",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(raster.cv <- ggplot(testDataset) +
  geom_point(aes(x = no2_measured_ug.m3, y = predict_no2, color = Density)) +
  scale_color_viridis() + 
  scale_x_continuous(name = "Measured NO2 (ug/m3)", limits = c(0, 100)) +
  scale_y_continuous(name = "Predicted NO2 (ug/m3)", limits = c(0, 100)) +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype = "dashed", size = 0.5) + 
  geom_abline(intercept = coeff[1], slope = coeff[2], color="blue", 
              size= 0.5) + 
  annotation_custom(grob) + 
  annotation_custom(grob.corr) +
  annotation_custom(grob.N) +
  annotation_custom(grob_add))

jpeg(file="07_Figure/gwpr.Raster.cv.jpeg", width = 210, height = 210, units = "mm", quality = 300, res = 300)
raster.cv
dev.off()

load("04_Results/trendenceMonthGroundLevel.RData")
pal <- colorRampPalette(c("blue", "white", "red"))
brks = c(-0.10, -0.08, -0.06, -0.04, -0.02,
         0, 0.02, 0.04, 0.06, 0.08, 0.10)
labels_brks = c("-10", "-8", "-6", "-4", "-2",
                "0", "2", "4", "6", "8", "10")
title_size = .0001
legend_title_size = 1
margin = 0

slope.tmap <- tm_shape(test.coeff.grid.raster.output) +
  tm_raster("month.slope", palette = pal(11), breaks = brks, 
            style = 'cont', legend.is.portrait = F, title = "The Slope of Monthly Chage\n[ *0.01 ug/(m3 * month)]",
            labels = labels_brks) +
  tm_grid(alpha = .25) + 
  tm_layout(
    inner.margins = c(margin, margin, margin, margin),
    title.size = title_size, 
    legend.position = c("left", "bottom"),
    legend.title.size = legend_title_size,
    legend.text.size = legend_title_size * 0.75
  ) + 
  tm_scale_bar()

#slope.tmap
slope.tmap %>%
  tmap_save(filename = "07_Figure/monthSlope.jpg", width = 300, height = 140, units = 'mm', dpi = 1000)

pal <- colorRampPalette(c("blue","green","yellow","red"))
load("04_Results/meanOfRasterNo2.RData")
#test.mean.grid.raster.output %>% plot()
brks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
labels_brks = brks %>% as.character()
mean.tmap <- tm_shape(test.mean.grid.raster.output) +
  tm_raster("ave.value", palette = pal(11), breaks = brks, 
            style = 'cont', legend.is.portrait = F, title = "The Mean of Monthly Concentration (ug/m3)",
            labels = labels_brks) +
  tm_grid(alpha = .25) + 
  tm_layout(
    inner.margins = c(margin, margin, margin, margin),
    title.size = title_size, 
    legend.position = c("left", "bottom"),
    legend.title.size = legend_title_size,
    legend.text.size = legend_title_size * 0.75
  ) + 
  tm_scale_bar()
mean.tmap %>%
  tmap_save(filename = "07_Figure/meanConcentration.jpg", width = 300, height = 140, units = 'mm', dpi = 1000)

predict_raster_folder <- "D:/10_Article/11_PredictRaster/01_Test0104/"
raster.filelist <- list.files(predict_raster_folder)
#brks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80)
brks = c(0, 10, 20, 30, 40, 50,  60, 70, 80, 90, 100)
labels_brks = brks %>% as.character()

for (file in raster.filelist){
  predict.tiff <- raster(paste0(predict_raster_folder, file))
  filename <- substr(file, 1, 6)
  tiff.tmap <- tm_shape(predict.tiff) +
    tm_raster(paste0("X", filename), palette = pal(17), breaks = brks, 
              style = 'cont', legend.is.portrait = F, 
              title = paste0("Monthly Concentration in ", filename, " (ug/m3)"),
              labels = labels_brks) +
    tm_grid(alpha = .25) + 
    tm_layout(
      inner.margins = c(margin, margin, margin, margin),
      title.size = title_size, 
      legend.position = c("left", "bottom"),
      legend.title.size = legend_title_size,
      legend.text.size = legend_title_size * 0.75
    ) + 
    tm_scale_bar()
  tiff.tmap %>%
    tmap_save(filename = paste0("07_Figure/month/", filename, ".jpg"),
              width = 300, height = 140, units = 'mm', dpi = 1000)
}

for (file in raster.filelist){
  predict.tiff <- raster(paste0(predict_raster_folder, file))
  filename <- substr(file, 1, 6)
  tiff.tmap <- tm_shape(predict.tiff) +
    tm_raster(paste0("X", filename), palette = pal(11), breaks = brks, 
              style = 'cont', legend.is.portrait = F, 
              title = paste0("Monthly Concentration in ", filename, " (ug/m3)"),
              labels = labels_brks) +
    tm_grid(alpha = .25) + 
    tm_layout(
      inner.margins = c(margin, margin, margin, margin),
      title.size = title_size, 
      legend.position = c("left", "bottom"),
      legend.title.size = legend_title_size,
      legend.text.size = legend_title_size * 0.75
    ) + 
    tm_scale_bar()
  tiff.tmap %>%
    tmap_save(filename = paste0("07_Figure/lowRes/", filename, ".jpg"),
              width = 300, height = 140, units = 'mm', dpi = 200)
}

jpg.list <- list.files("07_Figure/lowRes/")
frames <- paste0("07_Figure/lowRes/", jpg.list)
m <- image_read(frames)
m <- image_animate(m, fps = 2)
image_write(m, paste0("06_Animate/","ani.gif"))

load("04_Results/GWPR_BW_setp_list.Rdata")
GWPR.FEM.bandwidth.step.list <- GWPR.FEM.bandwidth.step.list %>% as.data.frame()
plot1 <- ggplot(GWPR.FEM.bandwidth.step.list, aes(x = BandwidthVector, y = ScoreVector)) +
  geom_point() +
  scale_x_continuous(name = "Fixed Distance Bandwidth (Arc Degree)") +
  scale_y_continuous(name = "Mean Square Prediction Error") +
  theme_bw()


load("04_Results/GWPR_Adaptive_BW_setp_list.Rdata")
GWPR.FEM.Adaptive.bandwidth.step.list <- GWPR.FEM.Adaptive.bandwidth.step.list %>% as.data.frame()
GWPR.FEM.Adaptive.bandwidth.step.list <- GWPR.FEM.Adaptive.bandwidth.step.list[3:98,]
plot2 <- ggplot(GWPR.FEM.Adaptive.bandwidth.step.list, aes(x = BandwidthVector, y = ScoreVector)) +
  geom_point() +
  scale_x_continuous(name = "Adaptive Distance Bandwidth") +
  scale_y_continuous(name = "Mean Square Prediction Error") +
  theme_bw()


jpeg(file="07_Figure/bwselection.jpeg", width = 297, height = 105, units = "mm", quality = 300, res = 300)
plot_grid(plot1, plot2, 
          labels = c("A", "B"), nrow = 1, ncol = 2)
dev.off()

# time series of data
load("04_Results/GWPR_FEM_CV_A_result.Rdata")
residual.GWPR <- GWPR.FEM.CV.A.result$GWPR.residuals 

load("03_Rawdata/usedDataset.RData")
merge.tropono2 <- usedDataset %>% dplyr::select(ug_m2_troposphere_no2,
                                                CityCode, period, year, month) %>%
  rename("id" = "CityCode")
residual.GWPR <- left_join(residual.GWPR, merge.tropono2, by = c("id", "period"))
residual.GWPR$ug_m2_troposphere_no2 <- residual.GWPR$ug_m2_troposphere_no2/160
monthly.GWPR <- aggregate(residual.GWPR %>% dplyr::select(y, yhat, ug_m2_troposphere_no2), 
                          list(residual.GWPR$period), FUN = mean)
monthly.GWPR$Group.2 <- monthly.GWPR$Group.1 %>% as.factor()
monthly.pre <- monthly.GWPR %>% dplyr::select(-Group.1 )
monthly.pre <- monthly.pre %>% pivot_longer(cols = c("y", "yhat", "ug_m2_troposphere_no2"))
monthly.GWPR.stderr <- aggregate(residual.GWPR %>% dplyr::select(y, yhat, ug_m2_troposphere_no2), 
                                 list(residual.GWPR$period), FUN = std.error)
monthly.GWPR.stderr$Group.2 <- monthly.GWPR.stderr$Group.1 %>% as.factor()
monthly.pre.stderr <- monthly.GWPR.stderr %>% dplyr::select(-Group.1 )
monthly.pre.stderr <- monthly.pre.stderr %>% pivot_longer(cols = c("y", "yhat", "ug_m2_troposphere_no2"))
monthly.pre.stderr <- monthly.pre.stderr %>% rename("SE" = "value")
monthly.pre <- left_join(monthly.pre, monthly.pre.stderr, by = c("Group.2", "name"))

(p3 <- ggplot(monthly.pre) +
    geom_col(aes(x = Group.2, y = value, fill = name), alpha = 0.7, stat = "identity",
             position = position_dodge(1), width = 0.8) +
    geom_errorbar(aes(x = Group.2, ymin = value - SE * 1.96, ymax = value + SE * 1.96, color = name), 
                  position = position_dodge(1), width = 0.5, size = 0.2) +
    scale_y_continuous(name = "Ground-Level NO2 Concentration",
                       sec.axis = sec_axis(~.*160, name="Tropospheric NO2 Concentration")) +
    scale_fill_manual(values = c("red", "blue", "gold1"), name =NULL, 
                      labels = c("Tropospheric NO2 Concentration", 
                                 "Measured Ground-Level NO2 Concentration",
                                 "Predicted Ground-Level NO2 Concentration")) +
    scale_color_manual(values = c("red4", "blue3", "salmon4")) +
    scale_x_discrete(name = NULL)+ 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = c(0.85, 0.88)) +
    guides(color = FALSE)
  )
jpeg(file="07_Figure/monthlymean.jpeg", width = 297, height = 150, units = "mm", quality = 300, res = 300)
p3
dev.off()
