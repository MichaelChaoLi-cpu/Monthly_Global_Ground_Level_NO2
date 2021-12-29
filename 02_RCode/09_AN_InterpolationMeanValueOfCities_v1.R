# end

library(tidyverse)
library(gstat)
library(sp) 
library(raster)
library(dplyr)
library(tmap)
library(automap)

### make the mean raster
load("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/03_Rawdata/mergedDataset.Rdata")

na.test <- mergedDataset %>% na.omit()
na.test$count <- 1
na.test <- aggregate(na.test$count, by = list(na.test$City, na.test$Country), FUN=sum)
colnames(na.test) <- c("City", "Country", "RecordCount")
na.test <- na.test %>% filter(RecordCount > 17) #freedom is 8

usedDataset <- left_join(mergedDataset, na.test, by = c("City", "Country"))
usedDataset <- usedDataset %>% filter(!is.na(RecordCount))
usedDataset <- usedDataset %>% dplyr::select(
  no2_measured_mg.m3,
  no2, mg_m2_total_no2, mg_m2_troposphere_no2,
  #mg_m2_total_no2_lag, mg_m2_troposphere_no2_lag,
  ter_pressure, dayTimeTemperature, nightTimeTemperature, ndvi,
  humidity, precipitation, NTL, speedwind, PBLH, 
  #UVAerosolIndex, ozone, cloudfraction, cloudpressure,
  CityCode, City, Country, 
  month, year, Date, Y2016, Y2017, Y2018, Y2019, Y2020, Y2021
) %>% na.omit()
usedDataset$humidity <- usedDataset$humidity %>% as.numeric()
usedDataset$precipitation <- usedDataset$precipitation %>% as.numeric()
usedDataset$year <- usedDataset$year %>% as.character() %>% as.numeric()
usedDataset$month <- usedDataset$month %>% as.character() %>% as.numeric()
usedDataset$period <- usedDataset$year * 100 + usedDataset$month
##### this is a strange value
usedDataset <- usedDataset %>% filter(CityCode != 499)
# preprocessing of the panel data set not we take the total column as the dependent variable

formula.CV.FEM <-
  no2_measured_mg.m3 ~ mg_m2_troposphere_no2 + ter_pressure + dayTimeTemperature + nightTimeTemperature +
  ndvi + humidity + precipitation + NTL + speedwind + PBLH + 
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021
rawCrossValidationDataset <- usedDataset %>% 
  dplyr::select("CityCode", "period", all.vars(formula.CV.FEM))
meanValueOfVariables <- stats::aggregate(rawCrossValidationDataset[,all.vars(formula.CV.FEM)],
                                         by = list(rawCrossValidationDataset$CityCode), mean)
colnames(meanValueOfVariables)[1] <- "CityCode"

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