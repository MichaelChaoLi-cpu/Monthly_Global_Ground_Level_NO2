# Author: M.L.

# input: mergedDataset.Rdata
# mergedDataset.Rdata "no2" monthly average no2 concentration ppm.
# mergedDataset.Rdata "mg_m2_total_no2" monthly average total column amount of no2 (mg/m2)
# mergedDataset.Rdata "mg_m2_troposphere_no2" monthly average tropospheric column amount of no2 (mg/m2)
# mergedDataset.Rdata "ter_pressure" monthly average terrain surface pressure (hpa)
# mergedDataset.Rdata "dayTimeTemperature" monthly average day time temperature (C)
# mergedDataset.Rdata "nightTimeTemperature" monthly average night time temperature (C)
# mergedDataset.Rdata "ndvi" NDVI -1 to 1
# mergedDataset.Rdata "humidity" g/kg means 1 gram water in the 1 kg air.
# mergedDataset.Rdata "precipitation" the precipitation unit is kg/(m2 * h)
# mergedDataset.Rdata "NTL" nighttime light
# mergedDataset.Rdata "CityCode" identity index.
# mergedDataset.Rdata "period" year * 100 + month, time index. 

# input CityLocationOfficial.csv
# CityLocationOfficial.csv: "Country", "City", "Latitude", "Longitude"

# end

library(tidyverse)
library(dplyr)
library(plm)
library(GWPR.light)
library(tmap)
library(sp)

load("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/03_Rawdata/mergedDataset.Rdata")

na.test <- mergedDataset %>% na.omit()
na.test$count <- 1
na.test <- aggregate(na.test$count, by = list(na.test$City, na.test$Country), FUN=sum)
colnames(na.test) <- c("City", "Country", "RecordCount")
na.test <- na.test %>% filter(RecordCount > 2)

usedDataset <- left_join(mergedDataset, na.test, by = c("City", "Country"))
usedDataset <- usedDataset %>% filter(!is.na(RecordCount))
usedDataset <- usedDataset %>% dplyr::select(
  no2, mg_m2_total_no2, ter_pressure, dayTimeTemperature, nightTimeTemperature, ndvi,
  humidity, precipitation, NTL, CityCode, period, City, Country
) %>% na.omit()
# preprocessing of the panel data set not we take the total column as the dependent variable

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

pdata <- pdata.frame(usedDataset, index = c("CityCode", "period"))
formula <- no2 ~ mg_m2_total_no2 + ter_pressure + dayTimeTemperature + nightTimeTemperature + ndvi +
  humidity + precipitation + NTL
ols <- plm(formula, pdata, model = "pooling")
summary(ols)
fem <- plm(formula, pdata, model = "within")
summary(fem)
rem <- plm(formula, pdata, model = "random")
summary(rem)
pFtest(fem, ols)
phtest(fem, rem)
plmtest(ols, type = c("bp"))
# base f test and hausman test fem is preferred global model.
# test linear model 

# we exiamine from the GWPR based on fem 
GWPR.FEM.bandwidth <- bw.GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                              SDF = cityLocationSpatialPoint, adaptive = F, p = 2, bigdata = F,
                              upperratio = 0.10, effect = "individual", model = "within", approach = "CV",
                              kernel = "bisquare",doParallel = T, cluster.number = 4)

GWPR.plmtest.Fixed.result <-
  GWPR.plmtest(formula = formula, data = usedDataset, index = c("CityCode", "period"),
               SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
               p = 2, kernel = "bisquare", longlat = F)
tm_shape(GWPR.plmtest.Fixed.result$SDF) +
  tm_dots(col = "p.value", breaks = c(0, 0.05, 1))
  ### this indicate that OLS is better than REM in most samples

GWPR.pFtest.Fixed.result <- 
  GWPR.pFtest(formula = formula, data = usedDataset, index = c("CityCode", "period"),
              SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
              p = 2, effect = "individual", kernel = "bisquare", longlat = F)
tm_shape(GWPR.pFtest.Fixed.result$SDF) +
  tm_dots(col = "p.value", breaks = c(0, 0.05, 1))
  ### this indicate that FEM is better than OLS in most samples

GWPR.phtest.Fixed.result <- 
  GWPR.phtest(formula = formula, data = usedDataset, index = c("CityCode", "period"),
              SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
              p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
              random.method = "amemiya")
tm_shape(GWPR.phtest.Fixed.result$SDF) +
  tm_dots(col = "p.value", breaks = c(0, 0.05, 1))
  ### this indicate that FEM is better than REM in most samples

GWPR.CV.F.result <- GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                         SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
                         p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                         model = "within")


# let us test pooled regression
GWPR.OLS.bandwidth <- bw.GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                              SDF = cityLocationSpatialPoint, adaptive = F, p = 2, bigdata = T,
                              upperratio = 0.10, effect = "individual", model = "pooling", approach = "CV",
                              kernel = "bisquare")
GWPR.OLS.CV.F.result <- GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                         SDF = cityLocationSpatialPoint, bw = GWPR.OLS.bandwidth, adaptive = F,
                         p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                         model = "pooling")

# let us test REM
GWPR.REM.bandwidth <- bw.GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                              SDF = cityLocationSpatialPoint, adaptive = F, p = 2, bigdata = T,
                              upperratio = 0.40, effect = "individual", model = "random", approach = "CV",
                              kernel = "bisquare", random.method = "swar")
