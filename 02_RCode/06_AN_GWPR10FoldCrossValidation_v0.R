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
# mergedDataset.Rdata "speedwind" the wind speed unit is m/s
# mergedDataset.Rdata "NTL" nighttime light
# mergedDataset.Rdata "PBLR" planetary boundary layer height unit is m.
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
usedDataset$period <- usedDataset$year * 100 + usedDataset$month
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
# get the city points 

#FEM Cross Validation, fixed bw 2.25 
formula.CV.FEM <-
  no2_measured_mg.m3 ~ mg_m2_troposphere_no2 + ter_pressure + dayTimeTemperature + nightTimeTemperature +
  ndvi + humidity + precipitation + NTL + speedwind + PBLH + 
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021 + 0

rawCrossValidationDataset <- usedDataset %>% 
  dplyr::select("CityCode", "period", all.vars(formula.CV.FEM))
meanValueOfVariables <- stats::aggregate(rawCrossValidationDataset[,all.vars(formula.CV.FEM)],
                                         by = list(rawCrossValidationDataset$CityCode), mean)
colnames(meanValueOfVariables)[1] <- "CityCode"
meanValueOfVariablesCity <- meanValueOfVariables
meanValueOfVariables <- dplyr::left_join(dplyr::select(rawCrossValidationDataset, "CityCode", "period"),
                                         meanValueOfVariables, by = "CityCode")
meanValueOfVariables <- meanValueOfVariables %>% arrange("CityCode", "period")
rawCrossValidationDataset <- rawCrossValidationDataset %>% arrange("CityCode", "period")

#### get FEM Transformation Dataset
femTransformationDataset <- (dplyr::select(rawCrossValidationDataset, -"CityCode", -"period")) - 
  (dplyr::select(meanValueOfVariables, -"CityCode", -"period"))
femTransformationDataset$CityCode <- rawCrossValidationDataset$CityCode
femTransformationDataset$period <- rawCrossValidationDataset$period

# Randomly order dataset
set.seed(42)

rows <- sample(nrow(femTransformationDataset))
femTransformationDataset <- femTransformationDataset[rows,]

singleFoldNumber <- floor(nrow(femTransformationDataset)/10)
foldNumberth <- 1

CV.result.table <- data.frame(Doubles=double(),
                              Ints=integer(),
                              Factors=factor(),
                              Logicals=logical(),
                              Characters=character(),
                              stringsAsFactors=FALSE)

while (foldNumberth < 11){
  if (foldNumberth == 10){
    rows.test <- rows[((foldNumberth-1)*singleFoldNumber+1):nrow(femTransformationDataset)]
  } else {
    rows.test <- rows[((foldNumberth-1)*singleFoldNumber+1):(foldNumberth*singleFoldNumber)]
  }
    
  test <- femTransformationDataset[rows.test,]
  train <- femTransformationDataset[-rows.test,]
  
  trainCode <- train %>%
    dplyr::select(CityCode) %>% distinct()
  
  trainCityLocation <- left_join(trainCode, cityLocation, by = "CityCode")
  xy <- trainCityLocation %>% dplyr::select(Longitude, Latitude)
  trainCityLocationSpatialPoint <- SpatialPointsDataFrame(coords = xy, data = cityLocation[,c(1, 2, 3, 4, 5)],
                                                          proj4string = CRS(proj))
  rm(xy)
  # get the train city points 
  
  GWPR.FEM.bandwidth = 2.25 ###
  GWPR.FEM.CV.F.result.CV1 <- GWPR.user(formula = formula.CV.FEM, data = train, index = c("CityCode", "period"),
                               SDF = trainCityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
                               p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                               model = "pooling")
  CVtrain.R2 <- GWPR.FEM.CV.F.result.CV1$R2
  coef.CV1 <- GWPR.FEM.CV.F.result.CV1$SDF@data
  coef.CV1 <- coef.CV1[,1:17]
  colnames(coef.CV1) <- paste0(colnames(coef.CV1), "_Coef")
  colnames(coef.CV1)[1] <- "CityCode"
  colnames(meanValueOfVariablesCity) <- paste0(colnames(meanValueOfVariablesCity), "_mean")
  colnames(meanValueOfVariablesCity)[1] <- "CityCode"
  
  test.predict <- left_join(test, coef.CV1, by = "CityCode")
  test.predict <- left_join(test.predict, meanValueOfVariablesCity, by = "CityCode")
  test.predict <- test.predict %>%
    mutate(predictNo2 = mg_m2_troposphere_no2_Coef * (mg_m2_troposphere_no2) + 
             ter_pressure_Coef * (ter_pressure) + 
             dayTimeTemperature_Coef * (dayTimeTemperature) +
             nightTimeTemperature_Coef * (nightTimeTemperature) +
             ndvi_Coef * (ndvi) +
             humidity_Coef * (humidity) +
             precipitation_Coef * (precipitation) +
             NTL_Coef * (NTL) + speedwind_Coef * (speedwind) +
             PBLH_Coef * (PBLH) +
             Y2016_Coef * (Y2016) + Y2017_Coef * (Y2017) +
             Y2018_Coef * (Y2018) + Y2019_Coef * (Y2019) +
             Y2020_Coef * (Y2020) + Y2021_Coef * (Y2021) + test.predict$no2_measured_mg.m3_mean
           )
  test.predict$no2_measured_mg.m3.ori <- test.predict$no2_measured_mg.m3 + test.predict$no2_measured_mg.m3_mean
  #ss.tot <- sum((test.predict$no2_measured_mg.m3.ori - mean(test.predict$no2_measured_mg.m3.ori))^2)
  ss.tot <- sum((test.predict$no2_measured_mg.m3.ori)^2)
  ss.res <- sum((test.predict$no2_measured_mg.m3.ori - test.predict$predictNo2)^2)
  CVtest.R2 <- 1 - ss.res/ss.tot
  result <- c(foldNumberth, CVtrain.R2, CVtest.R2)
  print(result)
  CV.result.table <- rbind(CV.result.table, result)
  foldNumberth <- foldNumberth + 1
}

colnames(CV.result.table) <- c("foldNumber", "CVtrain.R2", "CVtest.R2")