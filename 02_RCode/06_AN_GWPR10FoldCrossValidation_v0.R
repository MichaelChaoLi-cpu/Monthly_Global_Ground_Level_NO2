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

# input: CityLocationOfficial.csv
# CityLocationOfficial.csv: "Country", "City", "Latitude", "Longitude"

# output: femCrossValidation.Rdata
# femCrossValidation.Rdata: "foldNumber", "CVtrain.R2", "CVtest.R2" 

# end

library(tidyverse)
library(dplyr)
library(plm)
library(GWPR.light)
library(tmap)
library(sp)

load("03_Rawdata/usedDataset.Rdata")
load("04_Results/GWPR_FEM_CV_A_result.Rdata")
load("04_Results/GWPR_OLS_CV_A_result.Rdata")

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

#FEM Cross Validation, fixed bw 10
formula.CV.FEM <-
  no2_measured_mg.m3 ~ mg_m2_troposphere_no2 + ter_pressure + temp +
  ndvi + precipitation +  PBLH +
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
source("02_RCode/05_AF_GWPRRevisedForCrossValidation_v1.R")
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
  meanValueOfVariables.use <- meanValueOfVariables
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
  
  GWPR.FEM.bandwidth = GWPR.FEM.CV.F.result$GW.arguments$bw ###
  GWPR.FEM.CV.F.result.CV1 <- GWPR.user(formula = formula.CV.FEM, data = train, index = c("CityCode", "period"),
                               SDF = trainCityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
                               p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                               model = "pooling")
  #CVtrain.R2 <- GWPR.FEM.CV.F.result.CV1$R2
  coef.CV1 <- GWPR.FEM.CV.F.result.CV1$SDF@data
  coef.CV1 <- coef.CV1[,1:17]
  colnames(coef.CV1) <- paste0(colnames(coef.CV1), "_Coef")
  colnames(coef.CV1)[1] <- "CityCode"
  colnames(meanValueOfVariables.use) <- paste0(colnames(meanValueOfVariables.use), "_mean")
  colnames(meanValueOfVariables.use)[1] <- "CityCode"
  
  train.predict <- left_join(train, coef.CV1, by = "CityCode")
  train.predict <- left_join(train.predict, meanValueOfVariables.use, by = "CityCode")
  train.predict <- train.predict %>%
    mutate(predictNo2 = mg_m2_troposphere_no2_Coef * (mg_m2_troposphere_no2) + 
             ter_pressure_Coef * (ter_pressure) + 
             temp_Coef * temp +
             ndvi_Coef * (ndvi) +
             precipitation_Coef * (precipitation) +
             PBLH_Coef * (PBLH) +
             Y2016_Coef * (Y2016) + Y2017_Coef * (Y2017) +
             Y2018_Coef * (Y2018) + Y2019_Coef * (Y2019) +
             Y2020_Coef * (Y2020) + Y2021_Coef * (Y2021) + no2_measured_mg.m3_mean
    )
  train.predict$no2_measured_mg.m3.ori <- train.predict$no2_measured_mg.m3 + train.predict$no2_measured_mg.m3_mean
  #ss.tot <- sum((train.predict$no2_measured_mg.m3.ori - mean(train.predict$no2_measured_mg.m3.ori))^2)
  ss.tot <- sum((train.predict$no2_measured_mg.m3.ori - mean(test.predict$no2_measured_mg.m3))^2)
  ss.res <- sum((train.predict$no2_measured_mg.m3.ori - train.predict$predictNo2)^2)
  CVtrain.R2 <- 1 - ss.res/ss.tot
  
  test.predict <- left_join(test, coef.CV1, by = "CityCode")
  test.predict <- left_join(test.predict, meanValueOfVariables.use, by = "CityCode")
  test.predict <- test.predict %>%
    mutate(predictNo2 = mg_m2_troposphere_no2_Coef * (mg_m2_troposphere_no2) + 
             ter_pressure_Coef * (ter_pressure) + 
             temp_Coef * temp +
             ndvi_Coef * (ndvi) +
             precipitation_Coef * (precipitation) +
             PBLH_Coef * (PBLH) +
             Y2016_Coef * (Y2016) + Y2017_Coef * (Y2017) +
             Y2018_Coef * (Y2018) + Y2019_Coef * (Y2019) +
             Y2020_Coef * (Y2020) + Y2021_Coef * (Y2021) + no2_measured_mg.m3_mean
           )
  test.predict$no2_measured_mg.m3.ori <- test.predict$no2_measured_mg.m3 + test.predict$no2_measured_mg.m3_mean
  #ss.tot <- sum((test.predict$no2_measured_mg.m3.ori - mean(test.predict$no2_measured_mg.m3.ori))^2)
  ss.tot <- sum((test.predict$no2_measured_mg.m3.ori - mean(test.predict$no2_measured_mg.m3))^2)
  ss.res <- sum((test.predict$no2_measured_mg.m3.ori - test.predict$predictNo2)^2)
  CVtest.R2 <- 1 - ss.res/ss.tot
  result <- c(foldNumberth, CVtrain.R2, CVtest.R2)
  print(result)
  CV.result.table <- rbind(CV.result.table, result)
  foldNumberth <- foldNumberth + 1
}
colnames(CV.result.table) <- c("foldNumber", "CVtrain.R2", "CVtest.R2")
save(CV.result.table, file = "04_Results/femCrossValidation.Rdata")

#PoM Cross Validation, fixed bw 2.25 
formula.CV.PoM <-
  no2_measured_mg.m3 ~ mg_m2_troposphere_no2 + ter_pressure + temp +
  ndvi + precipitation +  PBLH +
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021
rawCrossValidationDataset <- usedDataset %>% 
  dplyr::select("CityCode", "period", all.vars(formula.CV.PoM))
rawCrossValidationDataset <- rawCrossValidationDataset %>% arrange("CityCode", "period")
pomTransformationDataset <- rawCrossValidationDataset
pomTransformationDataset <- pomTransformationDataset[rows,]

singleFoldNumber <- floor(nrow(pomTransformationDataset)/10)
foldNumberth <- 1

CV.result.pom.table <- data.frame(Doubles=double(),
                                  Ints=integer(),
                                  Factors=factor(),
                                  Logicals=logical(),
                                  Characters=character(),
                                  stringsAsFactors=FALSE)
while (foldNumberth < 11){
  if (foldNumberth == 10){
    rows.test <- rows[((foldNumberth-1)*singleFoldNumber+1):nrow(pomTransformationDataset)]
  } else {
    rows.test <- rows[((foldNumberth-1)*singleFoldNumber+1):(foldNumberth*singleFoldNumber)]
  }
  
  test <- pomTransformationDataset[rows.test,]
  train <- pomTransformationDataset[-rows.test,]
  
  trainCode <- train %>%
    dplyr::select(CityCode) %>% distinct()
  
  trainCityLocation <- left_join(trainCode, cityLocation, by = "CityCode")
  xy <- trainCityLocation %>% dplyr::select(Longitude, Latitude)
  trainCityLocationSpatialPoint <- SpatialPointsDataFrame(coords = xy, data = cityLocation[,c(1, 2, 3, 4, 5)],
                                                          proj4string = CRS(proj))
  rm(xy)
  # get the train city points 
  
  GWPR.PoM.bandwidth = GWPR.FEM.CV.F.result$GW.arguments$bw ###
  GWPR.PoM.CV.F.result.CV1 <- GWPR(formula = formula.CV.PoM, data = train, index = c("CityCode", "period"),
                                   SDF = trainCityLocationSpatialPoint, bw = GWPR.PoM.bandwidth, adaptive = F,
                                   p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                                   model = "pooling")
  CVtrain.R2 <- GWPR.PoM.CV.F.result.CV1$R2
  coef.CV1 <- GWPR.PoM.CV.F.result.CV1$SDF@data
  coef.CV1 <- coef.CV1[,1:18]
  colnames(coef.CV1) <- paste0(colnames(coef.CV1), "_Coef")
  colnames(coef.CV1)[1] <- "CityCode"
  
  test.predict <- left_join(test, coef.CV1, by = "CityCode")
  test.predict <- test.predict %>%
    mutate(predictNo2 = mg_m2_troposphere_no2_Coef * (mg_m2_troposphere_no2) + 
             ter_pressure_Coef * (ter_pressure) + 
             temp_Coef * temp +
             ndvi_Coef * (ndvi) +
             precipitation_Coef * (precipitation) +
             PBLH_Coef * (PBLH) +
             Y2016_Coef * (Y2016) + Y2017_Coef * (Y2017) +
             Y2018_Coef * (Y2018) + Y2019_Coef * (Y2019) +
             Y2020_Coef * (Y2020) + Y2021_Coef * (Y2021) +
             Intercept_Coef
    )
  ss.tot <- sum((test.predict$no2_measured_mg.m3 - mean(test.predict$no2_measured_mg.m3))^2)
  ss.res <- sum((test.predict$no2_measured_mg.m3 - test.predict$predictNo2)^2)
  CVtest.R2 <- 1 - ss.res/ss.tot
  result <- c(foldNumberth, CVtrain.R2, CVtest.R2)
  print(result)
  CV.result.pom.table <- rbind(CV.result.pom.table, result)
  foldNumberth <- foldNumberth + 1
}
colnames(CV.result.pom.table) <- c("foldNumber", "CVtrain.R2", "CVtest.R2")
save(CV.result.pom.table, file = "04_Results/pomCrossValidation.Rdata")


#### Adaptive bw FEM
set.seed(42)

rows <- sample(nrow(femTransformationDataset))
femTransformationDataset <- femTransformationDataset[rows,]

singleFoldNumber <- floor(nrow(femTransformationDataset)/10)
foldNumberth <- 1

CV.A.result.table <- data.frame(Doubles=double(),
                              Ints=integer(),
                              Factors=factor(),
                              Logicals=logical(),
                              Characters=character(),
                              stringsAsFactors=FALSE)

while (foldNumberth < 11){
  meanValueOfVariables.use <- meanValueOfVariables
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
  
  GWPR.FEM.bandwidth = GWPR.FEM.CV.A.result$GW.arguments$bw ###
  GWPR.FEM.CV.F.result.CV1 <- GWPR.user(formula = formula.CV.FEM, data = train, index = c("CityCode", "period"),
                                        SDF = trainCityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = T,
                                        p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                                        model = "pooling")
  #CVtrain.R2 <- GWPR.FEM.CV.F.result.CV1$R2
  coef.CV1 <- GWPR.FEM.CV.F.result.CV1$SDF@data
  coef.CV1 <- coef.CV1[,1:17]
  colnames(coef.CV1) <- paste0(colnames(coef.CV1), "_Coef")
  colnames(coef.CV1)[1] <- "CityCode"
  colnames(meanValueOfVariables.use) <- paste0(colnames(meanValueOfVariables.use), "_mean")
  colnames(meanValueOfVariables.use)[1] <- "CityCode"
  
  train.predict <- left_join(train, coef.CV1, by = "CityCode")
  train.predict <- left_join(train.predict, meanValueOfVariables.use, by = "CityCode")
  train.predict <- train.predict %>%
    mutate(predictNo2 = mg_m2_troposphere_no2_Coef * (mg_m2_troposphere_no2) + 
             ter_pressure_Coef * (ter_pressure) + 
             temp_Coef * temp +
             ndvi_Coef * (ndvi) +
             precipitation_Coef * (precipitation) +
             PBLH_Coef * (PBLH) +
             Y2016_Coef * (Y2016) + Y2017_Coef * (Y2017) +
             Y2018_Coef * (Y2018) + Y2019_Coef * (Y2019) +
             Y2020_Coef * (Y2020) + Y2021_Coef * (Y2021) + no2_measured_mg.m3_mean
    )
  train.predict$no2_measured_mg.m3.ori <- train.predict$no2_measured_mg.m3 + train.predict$no2_measured_mg.m3_mean
  #ss.tot <- sum((train.predict$no2_measured_mg.m3.ori - mean(train.predict$no2_measured_mg.m3.ori))^2)
  ss.tot <- sum((train.predict$no2_measured_mg.m3.ori - mean(test.predict$no2_measured_mg.m3))^2)
  ss.res <- sum((train.predict$no2_measured_mg.m3.ori - train.predict$predictNo2)^2)
  CVtrain.R2 <- 1 - ss.res/ss.tot
  
  test.predict <- left_join(test, coef.CV1, by = "CityCode")
  test.predict <- left_join(test.predict, meanValueOfVariables.use, by = "CityCode")
  test.predict <- test.predict %>%
    mutate(predictNo2 = mg_m2_troposphere_no2_Coef * (mg_m2_troposphere_no2) + 
             ter_pressure_Coef * (ter_pressure) + 
             temp_Coef * temp +
             ndvi_Coef * (ndvi) +
             precipitation_Coef * (precipitation) +
             PBLH_Coef * (PBLH) +
             Y2016_Coef * (Y2016) + Y2017_Coef * (Y2017) +
             Y2018_Coef * (Y2018) + Y2019_Coef * (Y2019) +
             Y2020_Coef * (Y2020) + Y2021_Coef * (Y2021) + no2_measured_mg.m3_mean
    )
  test.predict$no2_measured_mg.m3.ori <- test.predict$no2_measured_mg.m3 + test.predict$no2_measured_mg.m3_mean
  #ss.tot <- sum((test.predict$no2_measured_mg.m3.ori - mean(test.predict$no2_measured_mg.m3.ori))^2)
  ss.tot <- sum((test.predict$no2_measured_mg.m3.ori - mean(test.predict$no2_measured_mg.m3))^2)
  ss.res <- sum((test.predict$no2_measured_mg.m3.ori - test.predict$predictNo2)^2)
  CVtest.R2 <- 1 - ss.res/ss.tot
  result <- c(foldNumberth, CVtrain.R2, CVtest.R2)
  print(result)
  CV.A.result.table <- rbind(CV.A.result.table, result)
  foldNumberth <- foldNumberth + 1
}
colnames(CV.A.result.table) <- c("foldNumber", "CVtrain.R2", "CVtest.R2")
save(CV.A.result.table, file = "04_Results/AdaptivefemCrossValidation.Rdata")

#PoM Cross Validation, adaptive bw 7 
formula.CV.PoM <-
  no2_measured_mg.m3 ~ mg_m2_troposphere_no2 + ter_pressure + temp +
  ndvi + precipitation +  PBLH +
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021
rawCrossValidationDataset <- usedDataset %>% 
  dplyr::select("CityCode", "period", all.vars(formula.CV.PoM))
rawCrossValidationDataset <- rawCrossValidationDataset %>% arrange("CityCode", "period")
pomTransformationDataset <- rawCrossValidationDataset
pomTransformationDataset <- pomTransformationDataset[rows,]

singleFoldNumber <- floor(nrow(pomTransformationDataset)/10)
foldNumberth <- 1

CV.A.result.pom.table <- data.frame(Doubles=double(),
                                  Ints=integer(),
                                  Factors=factor(),
                                  Logicals=logical(),
                                  Characters=character(),
                                  stringsAsFactors=FALSE)
while (foldNumberth < 11){
  if (foldNumberth == 10){
    rows.test <- rows[((foldNumberth-1)*singleFoldNumber+1):nrow(pomTransformationDataset)]
  } else {
    rows.test <- rows[((foldNumberth-1)*singleFoldNumber+1):(foldNumberth*singleFoldNumber)]
  }
  
  test <- pomTransformationDataset[rows.test,]
  train <- pomTransformationDataset[-rows.test,]
  
  trainCode <- train %>%
    dplyr::select(CityCode) %>% distinct()
  
  trainCityLocation <- left_join(trainCode, cityLocation, by = "CityCode")
  xy <- trainCityLocation %>% dplyr::select(Longitude, Latitude)
  trainCityLocationSpatialPoint <- SpatialPointsDataFrame(coords = xy, data = cityLocation[,c(1, 2, 3, 4, 5)],
                                                          proj4string = CRS(proj))
  rm(xy)
  # get the train city points 
  
  GWPR.PoM.bandwidth = GWPR.OLS.CV.A.result$GW.arguments$bw ###
  GWPR.PoM.CV.F.result.CV1 <- GWPR(formula = formula.CV.PoM, data = train, index = c("CityCode", "period"),
                                   SDF = trainCityLocationSpatialPoint, bw = GWPR.PoM.bandwidth, adaptive = T,
                                   p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                                   model = "pooling")
  CVtrain.R2 <- GWPR.PoM.CV.F.result.CV1$R2
  coef.CV1 <- GWPR.PoM.CV.F.result.CV1$SDF@data
  coef.CV1 <- coef.CV1[,1:18]
  colnames(coef.CV1) <- paste0(colnames(coef.CV1), "_Coef")
  colnames(coef.CV1)[1] <- "CityCode"
  
  test.predict <- left_join(test, coef.CV1, by = "CityCode")
  test.predict <- test.predict %>%
    mutate(predictNo2 = mg_m2_troposphere_no2_Coef * (mg_m2_troposphere_no2) + 
             ter_pressure_Coef * (ter_pressure) + 
             temp_Coef * temp +
             ndvi_Coef * (ndvi) +
             precipitation_Coef * (precipitation) +
             PBLH_Coef * (PBLH) +
             Y2016_Coef * (Y2016) + Y2017_Coef * (Y2017) +
             Y2018_Coef * (Y2018) + Y2019_Coef * (Y2019) +
             Y2020_Coef * (Y2020) + Y2021_Coef * (Y2021) +
             Intercept_Coef
    )
  ss.tot <- sum((test.predict$no2_measured_mg.m3 - mean(test.predict$no2_measured_mg.m3))^2)
  ss.res <- sum((test.predict$no2_measured_mg.m3 - test.predict$predictNo2)^2)
  CVtest.R2 <- 1 - ss.res/ss.tot
  result <- c(foldNumberth, CVtrain.R2, CVtest.R2)
  print(result)
  CV.A.result.pom.table <- rbind(CV.A.result.pom.table, result)
  foldNumberth <- foldNumberth + 1
}
colnames(CV.A.result.pom.table) <- c("foldNumber", "CVtrain.R2", "CVtest.R2")
save(CV.A.result.pom.table, file = "04_Results/AdaptivePomCrossValidation.Rdata")
