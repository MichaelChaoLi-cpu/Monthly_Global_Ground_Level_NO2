# Author: M.L.

# end

library(dplyr)
library(tidyverse)
library(moments)
library(raster)

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
skewness(usedDataset.tranformed$ter_pressure_t)
skewness(usedDataset.tranformed$temp_t)
skewness(usedDataset.tranformed$ndvi_t)
skewness(usedDataset.tranformed$precipitation_t)
skewness(usedDataset.tranformed$PBLH_t)

test.dataset <- usedDataset.tranformed %>% 
  filter(no2_measured_ug.m3_t < 10.6 * 10)
skewness(test.dataset$no2_measured_ug.m3_t)
kurtosis(test.dataset$no2_measured_ug.m3_t)
test.dataset <- usedDataset.tranformed %>% 
  filter(no2_measured_ug.m3_t > 10.6 * 10)
test.dataset
test.dataset <- usedDataset.tranformed %>% 
  filter(ug_m2_troposphere_no2_t < 2751.1 * 10)
skewness(test.dataset$ug_m2_troposphere_no2_t)
test.dataset <- usedDataset.tranformed %>% 
  filter(ug_m2_troposphere_no2_t > 2751.1 * 7)

kurtosis(usedDataset.tranformed$no2_measured_ug.m3_t)
kurtosis(usedDataset.tranformed$ug_m2_troposphere_no2_t)
kurtosis(usedDataset.tranformed$ter_pressure_t)
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

#-------------change to line plot of monthlymean--------------
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

(p3.point <- ggplot(monthly.pre) +
    #geom_line(aes(x = Group.2, y = value, group = name, color = name),  stat = "identity",
    #          size = 1) +
    geom_errorbar(aes(x = Group.2, ymin = value - SE * 1.96, ymax = value + SE * 1.96), color = "gray88",
                  width = 0, size = 3, alpha = 0.8) +
    geom_point(aes(x = Group.2, y = value, group = name, color = name),  stat = "identity",
               size = 2, alpha = 0.8)+
    geom_errorbar(aes(x = Group.2, ymin = value - SE * 1.96, ymax = value + SE * 1.96, color = name), 
                  width = 0, size = 0.8, alpha = 0.6) +
    scale_y_continuous(name = "Ground-Level NO2 Concentration",
                       sec.axis = sec_axis(~.*160, name="Tropospheric NO2 Concentration")) +
    scale_color_manual(values = c("red", "blue", "gold1"), name =NULL, 
                       labels = c("Tropospheric NO2 Concentration", 
                                  "Measured Ground-Level NO2 Concentration",
                                  "Predicted Ground-Level NO2 Concentration")) +
    #scale_color_manual(values = c("red4", "blue3", "salmon4")) +
    scale_x_discrete(name = NULL)+ 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = c(0.85, 0.88))
)
jpeg(file="07_Figure/monthlymean_point_only.jpeg", width = 297, height = 180, units = "mm", quality = 300, res = 300)
p3.point
dev.off()

(p3.line <- ggplot(monthly.pre) +
    geom_line(aes(x = Group.2, y = value, group = name, color = name),  stat = "identity",
              size = 0.4) +
    geom_errorbar(aes(x = Group.2, ymin = value - SE * 1.96, ymax = value + SE * 1.96), color = "gray88",
                  width = 0, size = 3, alpha = 0.8) +
    geom_point(aes(x = Group.2, y = value, group = name, color = name),  stat = "identity",
               size = 2, alpha = 0.8) +
    geom_errorbar(aes(x = Group.2, ymin = value - SE * 1.96, ymax = value + SE * 1.96, color = name), 
                  width = 0, size = 0.8, alpha = 0.6) +
    scale_y_continuous(name = "Ground-Level NO2 Concentration",
                       sec.axis = sec_axis(~.*160, name="Tropospheric NO2 Concentration")) +
    scale_color_manual(values = c("red", "blue", "gold1"), name =NULL, 
                       labels = c("Tropospheric NO2 Concentration", 
                                  "Measured Ground-Level NO2 Concentration",
                                  "Predicted Ground-Level NO2 Concentration")) +
    #scale_color_manual(values = c("red4", "blue3", "salmon4")) +
    scale_x_discrete(name = NULL)+ 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = c(0.85, 0.88))
)
jpeg(file="07_Figure/monthlymean_line.jpeg", width = 297, height = 180, units = "mm", quality = 300, res = 300)
p3.line
dev.off()
#-------------change to line plot of monthlymean--------------

#### month trend raster and city global
rasterFolder <- "C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/10_DataArchive/"
rasterFileList <- list.files(rasterFolder)
mean_value_df <- data.frame(
  Doubles=double(),Ints=integer(), Factors=factor(), Logicals=logical(),
  Characters=character(), stringsAsFactors=FALSE)
for (rasterFile in rasterFileList){
  test_tiff <- raster::raster(paste0(rasterFolder, rasterFile))
  mean_value <- cellStats(test_tiff, stat='mean', na.rm = T)
  date <- rasterFile %>% substr(1, 6) %>% as.numeric()
  line <- c(mean_value, date)
  mean_value_df <- rbind(mean_value_df, line)
}

colnames(mean_value_df) <- c("mean_tiff", "date")
mean_value_df <- mean_value_df %>% arrange(date)

yhat_df <- monthly.pre %>% filter(name == 'yhat') %>%
  dplyr::select(Group.2, value)
yhat_df$Group.2 <- yhat_df$Group.2 %>% as.character() %>% as.numeric()
colnames(yhat_df) <- c("date", "mean_city")

mean_value_df <- left_join(mean_value_df, yhat_df, by = "date")
mean_value_df$month_order <- 1:82
date_vector <- mean_value_df$date
mean_value_df <- mean_value_df %>% 
  dplyr::select(-date)
mean_value_df <- mean_value_df %>% 
  pivot_longer(!month_order, names_to = "data_source", values_to = "value")
mean_value_df <- mean_value_df %>% na.omit()

(p4 <- ggplot(mean_value_df) +
    geom_point(aes(x = month_order, y = value, group = data_source, color = data_source),  stat = "identity",
               size = 2, alpha = 0.8) +
    geom_smooth(aes(x = month_order, y = value, group = data_source, color = data_source), method = lm) +
    scale_color_manual(values = c("red", "blue", "gold1"), name =NULL, 
                       labels = c("Mean of Predictions of the Cities", 
                                  "Mean of Predictions of All the Grids")) +
    scale_y_continuous(name = "NO2 Concentration") +
    scale_x_continuous(name = "Date", breaks = seq(1, 82, by = 1), labels = date_vector,
                       minor_breaks = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = c(0.85, 0.88))
  )
jpeg(file="07_Figure/meantiff_city_trend.jpeg", width = 297, height = 180, units = "mm", quality = 300, res = 300)
p4
dev.off()

lm(value ~ month_order, data = mean_value_df %>% filter(data_source == "mean_city")) %>%
  summary()
lm(value ~ month_order, data = mean_value_df %>% filter(data_source == "mean_tiff")) %>%
  summary()
