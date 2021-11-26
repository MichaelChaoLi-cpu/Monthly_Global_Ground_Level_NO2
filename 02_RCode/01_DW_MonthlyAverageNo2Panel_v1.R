# Author: M.L.

# output: PanelNo2Dataset.Rdata

# PanelNo2Dataset.Rdata: "no2" monthly average no2 concentration ppm.

# end

library(tidyverse)
library(dplyr)

setwd("D:/10_Article/01_RawData/11_AirPollutionRawTable/")

#2015.01 - 2021.1123 global AQI: https://aqicn.org/data-platform/covid19/verify/57c94d81-c396-481a-a12c-93a0065c705f
filelist <- list.files()
filelist <- filelist[2:length(filelist)] #the first file always the reminder.

totalAirPollutionDataset <- data.frame(Doubles=double(),
                                       Ints=integer(),
                                       Factors=factor(),
                                       Logicals=logical(),
                                       Characters=character(),
                                       stringsAsFactors=FALSE)
for (singlefile in filelist){
  Global_AQI <- read.csv(file = singlefile, skip = 4, encoding = "UTF-8")
  Global_AQI$year <- str_sub(Global_AQI$Date, 1, 4) %>% as.numeric()
  Global_AQI$month <- str_sub(Global_AQI$Date, 6, 7) %>% as.numeric()
  totalAirPollutionDataset <- rbind(totalAirPollutionDataset, Global_AQI)
  rm(Global_AQI)
}

totalNo2Dataset <- totalAirPollutionDataset %>%
  filter(Specie == "no2") %>% 
  filter(median < 500) %>% 
  dplyr::select(Country, City, median, year, month) 
totalNo2Dataset <- aggregate(totalNo2Dataset$median, 
                             by = list(totalNo2Dataset$Country, totalNo2Dataset$City,
                                       totalNo2Dataset$year, totalNo2Dataset$month), 
                             FUN = "mean", na.rm = T)
colnames(totalNo2Dataset) <- c("Country", "City", "year", "month", "no2")

rm(totalAirPollutionDataset)
rm(filelist)
rm(singlefile)

setwd("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub")
save.image("C:/Users/li.chao.987@s.kyushu-u.ac.jp/OneDrive - Kyushu University/10_Article/08_GitHub/03_Rawdata/PanelNo2Dataset.Rdata")
