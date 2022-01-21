# Author: M.L.

# input: GWPR_FEM_CV_A_result.Rdata
# note: this is the result of GWPR based on FEM. R2 is 0.7943. Adaptive bandwidth is 7.

# output: judgement.score.csv
# judgement.score.csv: "Year" the year.
# judgement.score.csv: "N" the year data size
# judgement.score.csv: "R2" R2, calculated following the article
# judgement.score.csv: "RMSE" RMSE, calculated following the article
# judgement.score.csv: "MAE" MAE, calculated following the article
# judgement.score.csv: "r" correlation coefficient
# judgement.score.csv: "Slope" the slope from the OLS between measured and predicted values
# judgement.score.csv: "Intercept" the intercept from the OLS between measured and predicted values

# end

library(tidyverse)
library(dplyr)

year.judgement.score <- function(residual.dataset){
  r2 <- 1 - sum( (residual.dataset$resid)^2 ) /
    sum((residual.dataset$y - mean(residual.dataset$y))^2)
  rmse <- sqrt(sum(residual.dataset$resid^2)/nrow(residual.dataset))
  mae <- mean(abs(residual.dataset$resid))
  cor.score <- cor.test(residual.dataset$y, residual.dataset$yhat)
  cor.score <- cor.score$estimate %>% as.numeric()
  reg <- lm(yhat ~ y, data = residual.dataset )
  coeff = coefficients(reg)
  N <- nrow(residual.dataset)
  year <- round(residual.dataset$period[1]/100,0)
  line.result <- c(year, N, r2, rmse, mae, cor.score, coeff[2], coeff[1])
  return(line.result)
}

load("04_Results/GWPR_FEM_CV_A_result.Rdata")

### get judgement score
residual.GWPR <- GWPR.FEM.CV.A.result$GWPR.residuals 
rmse <- sqrt(sum(residual.GWPR$resid^2)/nrow(residual.GWPR))
mae <- mean(abs(residual.GWPR$resid))

residual.GWPR.2015 <- residual.GWPR %>%
  filter(period < 201600)
residual.GWPR.2016 <- residual.GWPR %>%
  filter(period > 201600, period < 201700)
residual.GWPR.2017 <- residual.GWPR %>%
  filter(period > 201700, period < 201800)
residual.GWPR.2018 <- residual.GWPR %>%
  filter(period > 201800, period < 201900)
residual.GWPR.2019 <- residual.GWPR %>%
  filter(period > 201900, period < 202000)
residual.GWPR.2020 <- residual.GWPR %>%
  filter(period > 202000, period < 202100)
residual.GWPR.2021 <- residual.GWPR %>%
  filter(period > 202100)

line.2015 <- year.judgement.score(residual.GWPR.2015)
line.2016 <- year.judgement.score(residual.GWPR.2016)
line.2017 <- year.judgement.score(residual.GWPR.2017)
line.2018 <- year.judgement.score(residual.GWPR.2018)
line.2019 <- year.judgement.score(residual.GWPR.2019)
line.2020 <- year.judgement.score(residual.GWPR.2020)
line.2021 <- year.judgement.score(residual.GWPR.2021)
line.total <- year.judgement.score(residual.GWPR)

judgement.score <- rbind(line.2015, line.2016, line.2017, line.2018,
                         line.2019, line.2020, line.2021, line.total)
colnames(judgement.score) <- c("Year", "N", "R2", "RMSE", "MAE", "r", "Slope", "Intercept")
write.csv(judgement.score, file = "08_Tables/judgement.score.csv")
