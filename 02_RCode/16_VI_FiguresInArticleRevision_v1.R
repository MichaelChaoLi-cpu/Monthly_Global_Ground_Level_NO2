# Author: M.L.

# end

library(dplyr)
library(tidyverse)
library(moments)

load("03_Rawdata/usedDataset.RData")
formula <- no2_measured_ug.m3 ~ ug_m2_troposphere_no2 + 
  ter_pressure + temp + ndvi + precipitation + PBLH 

usedDataset.tranformed <- usedDataset %>% dplyr::select(CityCode, all.vars(formula))

usedDataset.tranformed.mean <- usedDataset.tranformed %>%
  aggregate(by = list(usedDataset.tranformed$CityCode), FUN = mean)
usedDataset.tranformed.mean <- usedDataset.tranformed.mean %>% dplyr::select(-Group.1)

colnames(usedDataset.tranformed.mean) <- paste0(colnames(usedDataset.tranformed.mean), "_m")
colnames(usedDataset.tranformed.mean)[1] <- "CityCode"

usedDataset.tranformed <- left_join(usedDataset.tranformed, usedDataset.tranformed.mean, by = "CityCode")
usedDataset.tranformed <- usedDataset.tranformed %>%
  mutate(no2_measured_ug.m3_t = no2_measured_ug.m3 - no2_measured_ug.m3_m,
         ug_m2_troposphere_no2_t = ug_m2_troposphere_no2 - ug_m2_troposphere_no2_m,
         ter_pressure_t = ter_pressure - ter_pressure_m,
         temp_t = temp - temp_m,
         ndvi_t = ndvi - ndvi_m,
         precipitation_t = precipitation - precipitation_m,
         PBLH_t = PBLH - PBLH_m)

#-------------descriptive statistics--------------
Mean <- round(mean(usedDataset.tranformed$no2_measured_ug.m3_t), 2)
SD <- round(sd(usedDataset.tranformed$no2_measured_ug.m3_t), 2)
N = nrow(usedDataset.tranformed)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("a",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(a <- ggplot(usedDataset.tranformed) +
    aes(x = no2_measured_ug.m3_t) +
    xlim(-50, 50) +
    geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(usedDataset.tranformed$no2_measured_ug.m3_t),
                              sd = sd(usedDataset.tranformed$no2_measured_ug.m3_t)),
                  col = 'red', size = 2) +
    xlab("Ground-level NO2 (ug/m3)") + 
    ylab("Density") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

Mean <- round(mean(usedDataset.tranformed$ug_m2_troposphere_no2_t), 2)
SD <- round(sd(usedDataset.tranformed$ug_m2_troposphere_no2_t), 2)
N = nrow(usedDataset.tranformed)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.72,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("b",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(b <- ggplot(usedDataset.tranformed) +
    aes(x = ug_m2_troposphere_no2_t) +
    geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(usedDataset.tranformed$ug_m2_troposphere_no2_t),
                              sd = sd(usedDataset.tranformed$ug_m2_troposphere_no2_t)),
                  col = 'red', size = 2) +
    xlim(-10000, 10000) +
    xlab("OMI TrCA NO2 (ug/m2)") + 
    ylab("Density") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

Mean <- round(mean(usedDataset.tranformed$ter_pressure_t), 2)
SD <- round(sd(usedDataset.tranformed$ter_pressure_t), 2)
N = nrow(usedDataset.tranformed)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("c",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(c <- ggplot(usedDataset.tranformed) +
    aes(x = ter_pressure_t) +
    geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(usedDataset.tranformed$ter_pressure_t),
                              sd = sd(usedDataset.tranformed$ter_pressure_t)),
                  col = 'red', size = 2) +
    xlim(-20, 20) +
    xlab("Terrain Atmospheric Pressure (hPa)") + 
    ylab("Density") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

Mean <- round(mean(usedDataset.tranformed$temp_t), 2)
SD <- round(sd(usedDataset.tranformed$temp_t), 2)
N = nrow(usedDataset.tranformed)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("d",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(d <- ggplot(usedDataset.tranformed) +
    aes(x = temp_t) +
    geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(usedDataset.tranformed$temp_t),
                              sd = sd(usedDataset.tranformed$temp_t)),
                  col = 'red', size = 2) +
    xlim(-30, 30) +
    xlab("Temperature (Celsius Degree)") + 
    ylab("Density") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

Mean <- round(mean(usedDataset.tranformed$ndvi_t), 2)
SD <- round(sd(usedDataset.tranformed$ndvi_t), 2)
N = nrow(usedDataset.tranformed)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("e",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(e <- ggplot(usedDataset.tranformed) +
    aes(x = ndvi_t) +
    geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(usedDataset.tranformed$ndvi_t),
                              sd = sd(usedDataset.tranformed$ndvi_t)),
                  col = 'red', size = 2) +
    xlim(-0.4, 0.4) +
    xlab("NDVI") + 
    ylab("Density") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

Mean <- round(mean(usedDataset.tranformed$precipitation_t), 2)
SD <- round(sd(usedDataset.tranformed$precipitation_t), 2)
N = nrow(usedDataset.tranformed)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("f",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(f <- ggplot(usedDataset.tranformed) +
    aes(x = precipitation_t) +
    geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(usedDataset.tranformed$precipitation_t),
                              sd = sd(usedDataset.tranformed$precipitation_t)),
                  col = 'red', size = 2) +
    xlim(-0.4, 0.4) +
    xlab("Precipitation kg/(m2*h)") + 
    ylab("Denisty") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

Mean <- round(mean(usedDataset.tranformed$PBLH_t), 2)
SD <- round(sd(usedDataset.tranformed$PBLH_t), 2)
N = nrow(usedDataset.tranformed)
grob <- grobTree(textGrob(paste0("Mean = ", Mean, "\nStd.dev = ", SD,"\nN = ", N),
                          x = 0.75,  y = 0.90, hjust = 0,
                          gp = gpar(col = "black", fontsize = 8)))
grob_add <- grobTree(textGrob("g",
                              x = 0.02,  y = 0.95, hjust = 0,
                              gp = gpar(col = "black", fontsize = 18)))
(g <- ggplot(usedDataset.tranformed) +
    aes(x = PBLH_t) +
    geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(usedDataset.tranformed$PBLH_t),
                              sd = sd(usedDataset.tranformed$PBLH_t)),
                  col = 'red', size = 2) +
    #xlim(0, 2500) +
    xlab("PBLH (m)") + 
    ylab("Denisty") +
    annotation_custom(grob) +
    annotation_custom(grob_add))

jpeg(file="07_Figure\\descriptive_stat_transform.jpeg", width = 297, height = 315, units = "mm", quality = 300, res = 300)
grid.arrange(a, b, c,
             d, e, f,
             g,
             nrow = 3)
dev.off()
#-------------descriptive statistics--------------

#### check skewness
skewness(usedDataset.tranformed$no2_measured_ug.m3_t)
skewness(usedDataset.tranformed$ug_m2_troposphere_no2_t)
skewness(usedDataset.tranformed$ter_pressure)
skewness(usedDataset.tranformed$temp_t)
skewness(usedDataset.tranformed$ndvi_t)
skewness(usedDataset.tranformed$precipitation_t)
skewness(usedDataset.tranformed$PBLH_t)

kurtosis(usedDataset.tranformed$no2_measured_ug.m3_t)
kurtosis(usedDataset.tranformed$ug_m2_troposphere_no2_t)
kurtosis(usedDataset.tranformed$ter_pressure)
kurtosis(usedDataset.tranformed$temp_t)
kurtosis(usedDataset.tranformed$ndvi_t)
kurtosis(usedDataset.tranformed$precipitation_t)
kurtosis(usedDataset.tranformed$PBLH_t)

write.csv(usedDataset.tranformed, file = "03_Rawdata/00_transformed_total_dataset.csv")
#### check skewness

#-------------trends of XY------------------------
(plot1 <- ggscatter(usedDataset.tranformed, x = "ug_m2_troposphere_no2_t", y = "no2_measured_ug.m3_t", size = 1,
                    add = "reg.line", conf.int = TRUE,
                    cor.coef = F, cor.method = "pearson",
                    xlab = "OMI TrCA NO2", ylab = "Ground-level NO2",           
                    add.params = list(color = "blue", fill = "lightskyblue1"),
                    color = "grey76", shape = 21, ylim = c(-25, 50), xlim = c(-10000, 25000)
) +
  stat_cor( p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = 0.03, label.y.npc = 0.17) +
  annotate("text", x = 23000, y = 46, label = 'bold("a")', parse = TRUE, size = 5)
)


(plot2 <- ggscatter(usedDataset.tranformed, x = "ter_pressure_t", y = "no2_measured_ug.m3_t", size = 1,
                    add = "reg.line", conf.int = TRUE,
                    cor.coef = F, cor.method = "pearson",
                    xlab = "Terrain Atmospheric Pressure", ylab = "Ground-level NO2",
                    add.params = list(color = "blue", fill = "lightskyblue1"),
                    color = "grey76", shape = 21, ylim = c(-50, 50), xlim = c(-15, 15)) +
    stat_cor( p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = 0.23, label.y.npc = 0.17) +
    annotate("text", x = 13, y = 46, label = 'bold("b")', parse = TRUE, size = 5)
  )

(plot3 <- ggscatter(usedDataset.tranformed, x = "temp_t", y = "no2_measured_ug.m3_t", size = 1,
                    add = "reg.line", conf.int = TRUE,
                    cor.coef = F, cor.method = "pearson",
                    xlab = "Temperature", ylab = "Ground-level NO2",
                    add.params = list(color = "blue", fill = "lightskyblue1"),
                    color = "grey76", shape = 21,ylim = c(-50, 50)) +
    stat_cor( p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = "left", label.y.npc = 0.17) +
    annotate("text", x = 21, y = 46, label = 'bold("c")', parse = TRUE, size = 5)
  )

(plot4 <- ggscatter(usedDataset.tranformed, x = "ndvi_t", y = "no2_measured_ug.m3_t", size = 1,
                    add = "reg.line", conf.int = TRUE,
                    cor.coef = F, cor.method = "pearson",
                    xlab = "NDVI", ylab = "Ground-level NO2",
                    add.params = list(color = "blue", fill = "lightskyblue1"),
                    color = "grey76", shape = 21, ylim = c(-50, 50)) +
    stat_cor( p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = "left", label.y.npc = 0.17) +
    annotate("text", x = 0.4, y = 46, label = 'bold("d")', parse = TRUE, size = 5)
  )

(plot5 <- ggscatter(usedDataset.tranformed, x = "precipitation_t", y = "no2_measured_ug.m3_t", size = 1,
                    add = "reg.line", conf.int = TRUE,
                    cor.coef = F, cor.method = "pearson",
                    xlab = "Precipitation", ylab = "Ground-level NO2",
                    add.params = list(color = "blue", fill = "lightskyblue1"),
                    color = "grey76", shape = 21, ylim = c(-50, 50)) +
    stat_cor( p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = "left", label.y.npc = 0.17) +
    annotate("text", x = 0.75, y = 46, label = 'bold("e")', parse = TRUE, size = 5)
  )

(plot6 <- ggscatter(usedDataset.tranformed, x = "PBLH_t", y = "no2_measured_ug.m3_t", size = 1,
                    add = "reg.line", conf.int = TRUE,
                    cor.coef = F, cor.method = "pearson",
                    xlab = "PBLH", ylab = "Ground-level NO2",
                    add.params = list(color = "blue", fill = "lightskyblue1"),
                    color = "grey76", shape = 21, ylim = c(-50, 50)) +
    stat_cor( p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = "left", label.y.npc = 0.17) +
    annotate("text", x = 1200, y = 46, label = 'bold("f")', parse = TRUE, size = 5)
  )

jpeg(file="07_Figure/cor_line1_transform.jpeg", width = 297, height = 210, units = "mm", quality = 300, res = 300)
grid.arrange(plot1, plot2, plot3, 
             plot4, plot5, plot6,
             nrow = 2)
dev.off()
#-------------trends of XY------------------------