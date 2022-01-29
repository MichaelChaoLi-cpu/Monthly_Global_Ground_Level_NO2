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
# mergedDataset.Rdata "PBLH" planetary boundary layer height unit is m.
# mergedDataset.Rdata "CityCode" identity index.
# mergedDataset.Rdata "period" year * 100 + month, time index. 
# mergedDataset.Rdata "ug_m2_total_no2" monthly average total column amount of no2 (ug/m2)
# mergedDataset.Rdata "ug_m2_troposphere_no2" monthly average tropospheric column amount of no2 (ug/m2)
# mergedDataset.Rdata "no2_measured_ug.m3" (ug/m3)

# input CityLocationOfficial.csv
# CityLocationOfficial.csv: "Country", "City", "Latitude", "Longitude"

# output: usedDataset.RData
# usedDataset.RData: "no2_measured_mg.m3", "no2", "mg_m2_total_no2", "mg_m2_troposphere_no2",
#                    "ter_pressure", "dayTimeTemperature", "nightTimeTemperature", "ndvi",
#                    "precipitation", "NTL", "PBLH", "CityCode", "City", "Country","month",
#                    "year", "Date", "Y2016", "Y2017", "Y2018", "Y2019", "Y2020", "Y2021",
#                    "period"
# usedDataset.RData: "temp" the average temperature of "dayTimeTemperature" and "nightTimeTemperature"

# output: GWPR_BW_setp_list.Rdata
# GWPR_BW_setp_list.Rdata: "BandwidthVector" from 0.25 to 50, step length is 0.25.
# GWPR_BW_setp_list.Rdata: "ScoreVector" CV score.

# output: GWPR_FEM_CV_F_result.Rdata
# note: this is the result of GWPR based on FEM. R2 is 0.8167. Fixed bandwidth is 2.25 arc degrees.
#       But islands appear.

# output: GWPR_OLS_CV_F_result.Rdata
# note: this is the result of GWPR based on OLS. R2 is 0.7436. Fixed bandwidth is 2.25 arc degrees.
#       But islands appear

# output: GWPR_FEM_CV_A_result.Rdata
# note: this is the result of GWPR based on FEM. R2 is 0.7445. Adaptive bandwidth is 7.

# output: GWPR_OLS_CV_A_result.Rdata
# note: this is the result of GWPR based on OLS. R2 is 0.7004. Adaptive bandwidth is 7.

# Note: "mg_m2_troposphere_no2" is with high accuracy.
# Note: in this version, we drop "NTL", because low accuary of interpolation in both coefficient and 
#       mean value. After dropping this variable the R2 only reduces 0.1%. So, we decide to do so. Bandwidth
#       is 10 arc degrees.
# Note: we test 2.25 degrees, this time.

# end

library(tidyverse)
library(dplyr)
library(plm)
library(GWPR.light)
library(tmap)
library(sp)
library(doParallel)
library(foreach)

load("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/03_Rawdata/mergedDataset.Rdata")

na.test <- mergedDataset %>% dplyr::select(
  no2_measured_ug.m3,
  no2, ug_m2_total_no2, ug_m2_troposphere_no2,
  #ug_m2_total_no2_lag, ug_m2_troposphere_no2_lag,
  ter_pressure, dayTimeTemperature, nightTimeTemperature, ndvi,
  precipitation, PBLH, #speedwind, humidity,  NTL,
  #UVAerosolIndex, ozone, cloudfraction, cloudpressure,
  CityCode, City, Country, 
  month, year, Date, Y2016, Y2017, Y2018, Y2019, Y2020, Y2021
) %>% na.omit()
na.test$count <- 1
na.test <- aggregate(na.test$count, by = list(na.test$City, na.test$Country), FUN=sum)
colnames(na.test) <- c("City", "Country", "RecordCount")
na.test <- na.test %>% filter(RecordCount > 5) #freedom 

usedDataset <- left_join(mergedDataset, na.test, by = c("City", "Country"))
usedDataset <- usedDataset %>% filter(!is.na(RecordCount))
usedDataset <- usedDataset %>% dplyr::select(
  no2_measured_ug.m3,
  no2, ug_m2_total_no2, ug_m2_troposphere_no2,
  #ug_m2_total_no2_lag, ug_m2_troposphere_no2_lag,
  ter_pressure, dayTimeTemperature, nightTimeTemperature, ndvi,
  precipitation, PBLH, #speedwind, humidity,  NTL,
  #UVAerosolIndex, ozone, cloudfraction, cloudpressure,
  CityCode, City, Country, 
  month, year, Date, Y2016, Y2017, Y2018, Y2019, Y2020, Y2021
) %>% na.omit()
#usedDataset$humidity <- usedDataset$humidity %>% as.numeric()
usedDataset$precipitation <- usedDataset$precipitation %>% as.numeric()
usedDataset$year <- usedDataset$year %>% as.character() %>% as.numeric()
usedDataset$month <- usedDataset$month %>% as.character() %>% as.numeric()
usedDataset$period <- usedDataset$year * 100 + usedDataset$month
usedDataset$temp <- (usedDataset$dayTimeTemperature + usedDataset$nightTimeTemperature) / 2
usedDataset.ori <- usedDataset
#this is a strange city, island
usedDataset <- usedDataset %>% 
  filter(CityCode != 499, CityCode != 4) ###CityCode != 291, CityCode != 163,
# preprocessing of the panel data set not we take the total column as the dependent variable

cityLocation <- read.csv("D:/10_Article/01_RawData/12_LocationJson/CityLocationOfficial.csv",
                         encoding="UTF-8") %>%
  dplyr::select(X0, X1, X2, X3)
colnames(cityLocation) <- c("Latitude", "Longitude", "City", "Country")
cityLocation <- cityLocation %>% 
  mutate(City = ifelse(City == "Washington, D.C.", "Washington D.C.", City))
cityNameCode <- usedDataset %>% dplyr::select(CityCode, City, Country) %>% distinct()
cityLocation <- left_join(cityNameCode, cityLocation, by = c("City", "Country"))
cityLocation.out <- cityLocation[order(cityLocation$Country, cityLocation$City),]
cityLocation.out %>%
  readr::write_csv(
    file = "C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/08_Tables/cityLocationTable.csv")
rm(cityNameCode)

proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
xy <- cityLocation %>% dplyr::select(Longitude, Latitude)
cityLocationSpatialPoint <- SpatialPointsDataFrame(coords = xy, data = cityLocation[,c(1, 2, 3, 4, 5)],
                                                   proj4string = CRS(proj))
rm(xy)
# get the city points 

# distance matrix
drop_island <- F # this time we test adaptive bandwidth
if (drop_island) {
  coord.point <- coordinates(cityLocationSpatialPoint)
  dist.matrix <- GWmodel::gw.dist(dp.locat = coord.point, focus=0, p=2, theta=0, longlat=F)
  dist.matrix <- dist.matrix %>% as.data.frame()
  for (i in 1:nrow(dist.matrix)){
    dist.matrix[i, i] <- Inf
  }
  dist.matrix$min <- rje::rowMins(dist.matrix) 
  dist.matrix$CityCode <- cityLocationSpatialPoint$CityCode
  dist.matrix <- dist.matrix %>%
    dplyr::select(min, CityCode)
  dist.matrix <- dist.matrix %>%
    filter(min < 2.25)
# distance matrix

##### this is a strange value
  usedDataset <- left_join(usedDataset, dist.matrix, by = "CityCode")
  usedDataset <- usedDataset %>%
    filter(!is.na(min))
}

##### re-obtain the point
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
plot(cityLocationSpatialPoint)
# get the city points 
save(usedDataset, file = "03_Rawdata/usedDataset.RData")
save(cityLocationSpatialPoint, file = "03_Rawdata/cityLocationSpatialPoint.Rdata")

na.test <- usedDataset %>% na.omit()
na.test$count <- 1
na.test <- aggregate(na.test$count, by = list(na.test$City, na.test$Country), FUN=sum)
colnames(na.test) <- c("City", "Country", "RecordCount")

pdata <- pdata.frame(usedDataset, index = c("CityCode", "Date"))
formula <- no2_measured_ug.m3 ~ ug_m2_troposphere_no2 + 
  ter_pressure + 
  temp +
  ndvi + precipitation + PBLH + 
  #humidity + UVAerosolIndex + ozone + speedwind + NTL +
  #cloudfraction + cloudpressure + # add this two variables effect are limited, only increase 0.2% R2
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021
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


setwd("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/")
source("02_RCode/07_AF_GWPRBandwidthStepSelection_v1.R")
# we exiamine from the GWPR based on fem 
GWPR.FEM.bandwidth <- # this is about fixed bandwidth
  bw.GWPR.step.selection(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                              SDF = cityLocationSpatialPoint, adaptive = F, p = 2, bigdata = F,
                              upperratio = 0.10, effect = "individual", model = "within", approach = "CV",
                              kernel = "bisquare",doParallel = T, cluster.number = 6, gradientIncrecement = T,
                              GI.step = 0.25, GI.upper = 20, GI.lower = 0.25)
GWPR.FEM.bandwidth.step.list <- GWPR.FEM.bandwidth
plot(GWPR.FEM.bandwidth.step.list[,1], GWPR.FEM.bandwidth.step.list[,2])
save(GWPR.FEM.bandwidth.step.list,
     file = "C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/04_Results/GWPR_BW_setp_list.Rdata")

GWPR.FEM.bandwidth.golden <- 
  bw.GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
          SDF = cityLocationSpatialPoint, adaptive = F, p = 2, bigdata = T,
          upperratio = 0.15, effect = "individual", model = "within", approach = "CV",
          kernel = "bisquare",doParallel = T, cluster.number = 6)
GWPR.FEM.bandwidth = 2.25 ###

################################ this is GWPR based on FEM
GWPR.FEM.CV.F.result <- GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                             SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = T,
                             p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                             model = "within")
GWPR.FEM.CV.F.result$SDF@data %>% view()
save(GWPR.FEM.CV.F.result, file = "C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/04_Results/GWPR_FEM_CV_F_result.Rdata")

# let us test pooled regression
GWPR.OLS.CV.F.result <- GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                             SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
                             p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                             model = "pooling")
save(GWPR.OLS.CV.F.result, file = "C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/04_Results/GWPR_OLS_CV_F_result.Rdata")

run <- F
if (run){
  # let us test REM #### random effect fail
  GWPR.REM.CV.F.result <- GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                               SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
                               p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                               model = "random", random.method = "amemiya")
  #note: the REM requires very high freedom. therefore, the bandwidth should be super large, since the points are
  #      not evenly distributed.
  
  formula.total <- no2_measured_ug.m3 ~ ug_m2_total_no2 + 
    ter_pressure + ndvi +  precipitation + NTL + PBLH +
    #UVAerosolIndex + ozone + dayTimeTemperature + nightTimeTemperature + humidity + speedwind +
    #cloudfraction + cloudpressure + # add this two variables effect are limited, only increase 0.2% R2
    Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021
  GWPR.FEM.CV.F.result <- GWPR(formula = formula.total, data = usedDataset, index = c("CityCode", "period"),
                               SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
                               p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                               model = "within")
}


GWPR.FEM.Adaptive.bandwidth <- # this is about adaptive bandwidth
  bw.GWPR.step.selection(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                         SDF = cityLocationSpatialPoint, adaptive = T, p = 2, bigdata = F,
                         upperratio = 0.10, effect = "individual", model = "within", approach = "CV",
                         kernel = "bisquare",doParallel = T, cluster.number = 6, gradientIncrecement = T,
                         GI.step = 1, GI.upper = 100, GI.lower = 2)
GWPR.FEM.Adaptive.bandwidth.step.list <- GWPR.FEM.Adaptive.bandwidth
plot(GWPR.FEM.Adaptive.bandwidth.step.list[,1], GWPR.FEM.Adaptive.bandwidth.step.list[,2])
save(GWPR.FEM.Adaptive.bandwidth.step.list,
     file = "C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/04_Results/GWPR_Adaptive_BW_setp_list.Rdata")
# turning point is 7 and 22. minimum is 2
# we test them both

run <- F
if(run){
  GWPR.FEM.Adaptive.bandwidth.golden <- 
    bw.GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
            SDF = cityLocationSpatialPoint, adaptive = T, p = 2, bigdata = T,
            upperratio = 0.15, effect = "individual", model = "within", approach = "CV",
            kernel = "bisquare",doParallel = T, cluster.number = 6)
}

GWPR.FEM.CV.A.result <- GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                             SDF = cityLocationSpatialPoint, bw = 7, adaptive = T,
                             p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                             model = "within")
run <- F
if(run){
  GWPR.FEM.CV.A.result <- GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                               SDF = cityLocationSpatialPoint, bw = 22, adaptive = T,
                               p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                               model = "within")
}
GWPR.FEM.CV.A.result$SDF@data %>% view()
save(GWPR.FEM.CV.A.result, file = "C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/04_Results/GWPR_FEM_CV_A_result.Rdata")

# let us test pooled regression
GWPR.OLS.CV.A.result <- GWPR(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                             SDF = cityLocationSpatialPoint, bw = 7, adaptive = T,
                             p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                             model = "pooling")
save(GWPR.OLS.CV.A.result, file = "C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/04_Results/GWPR_OLS_CV_A_result.Rdata")

formula.total <- no2_measured_ug.m3 ~ ug_m2_total_no2 + 
  ter_pressure + temp +
  ndvi + precipitation + PBLH +
  Y2016 + Y2017 + Y2018 + Y2019 + Y2020 + Y2021
GWPR.FEM.CV.A.result.total <- GWPR(formula = formula.total, data = usedDataset, index = c("CityCode", "period"),
                             SDF = cityLocationSpatialPoint, bw = 7, adaptive = T,
                             p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                             model = "within")

run <- F
if(run){
  GWPR.plmtest.Fixed.result <-
    GWPR.plmtest(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                 SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
                 p = 2, kernel = "bisquare", longlat = F)
  tm_shape(GWPR.plmtest.Fixed.result$SDF) +
    tm_dots(col = "p.value", breaks = c(0, 0.1, 1))
  ### this indicate that OLS is better than REM in most samples
  
  GWPR.pFtest.Fixed.result <- 
    GWPR.pFtest(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
                p = 2, effect = "individual", kernel = "bisquare", longlat = F)
  tm_shape(GWPR.pFtest.Fixed.result$SDF) +
    tm_dots(col = "p.value", breaks = c(0, 0.1, 1))
  ### this indicate that FEM is better than OLS in most samples
  
  GWPR.phtest.Fixed.result <- 
    GWPR.phtest(formula = formula, data = usedDataset, index = c("CityCode", "period"),
                SDF = cityLocationSpatialPoint, bw = GWPR.FEM.bandwidth, adaptive = F,
                p = 2, effect = "individual", kernel = "bisquare", longlat = F, 
                random.method = "amemiya")
  tm_shape(GWPR.phtest.Fixed.result$SDF) +
    tm_dots(col = "p.value", breaks = c(0, 0.05, 1))
  ### this indicate that FEM is better than REM in most samples
}

# the result are included in the GWPR.result$resid
rawCrossValidationDataset <- usedDataset %>% 
  dplyr::select("CityCode", "period", all.vars(formula))
meanValueOfVariables <- stats::aggregate(rawCrossValidationDataset[,all.vars(formula)],
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

meanValueOfVariables.use <- meanValueOfVariables
coef.CV1 <- GWPR.FEM.CV.A.result$SDF@data
coef.CV1 <- coef.CV1[,1:17]
colnames(coef.CV1) <- paste0(colnames(coef.CV1), "_Coef")
colnames(coef.CV1)[1] <- "CityCode"
colnames(meanValueOfVariables.use) <- paste0(colnames(meanValueOfVariables.use), "_mean")
colnames(meanValueOfVariables.use)[1] <- "CityCode"
meanValueOfVariables.use <- meanValueOfVariables.use %>% dplyr::select(-"period_mean") %>% distinct()

data.predict <- left_join(femTransformationDataset, coef.CV1, by = "CityCode")
data.predict <- left_join(data.predict, meanValueOfVariables.use, by = "CityCode")
data.predict <- data.predict %>%
  mutate(predictNo2 = ug_m2_troposphere_no2_Coef * (ug_m2_troposphere_no2) + 
           ter_pressure_Coef * (ter_pressure) + 
           temp_Coef * temp +
           ndvi_Coef * (ndvi) +
           precipitation_Coef * (precipitation) +
           PBLH_Coef * (PBLH) +
           Y2016_Coef * (Y2016) + Y2017_Coef * (Y2017) +
           Y2018_Coef * (Y2018) + Y2019_Coef * (Y2019) +
           Y2020_Coef * (Y2020) + Y2021_Coef * (Y2021) + no2_measured_ug.m3_mean
  )
data.predict$no2_measured_ug.m3.ori <- data.predict$no2_measured_ug.m3 + data.predict$no2_measured_ug.m3_mean
cor.test(data.predict$predictNo2, data.predict$no2_measured_ug.m3.ori)
save(data.predict, file = "04_Results/dataPrediction.Rdata")
