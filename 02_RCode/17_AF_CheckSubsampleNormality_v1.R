# Author: M.L.

# end

gwpr_data_prepare <- function(formula, data, index, SDF, bw = NULL, adaptive = FALSE, p = 2,
                              effect = "individual", model = c("pooling", "within", "random"), random.method = "swar",
                              kernel = "bisquare", longlat = FALSE)
{
  if(length(index) != 2)
  {
    stop("The \"index\" have included \"ID\" or/and \"time\" index.")
  }
  if(!((index[1] %in% colnames(data)) & (index[2] %in% colnames(data))))
  {
    stop("The data.frame(data) does not have the index columns.")
  }
  if(!(index[1] %in% colnames(SDF@data)))
  {
    stop("The SDF does not have the \"ID\" columns.")
  }
  if(is.null(bw))
  {
    stop("The bw must be set.")
  }
  if(!(model %in% c("pooling", "within", "random")))
  {
    stop("This version GWPR only accept \"pooling\", \"within\", or \"random\"")
  }
  
  # Data preparation
  varibale_name_in_equation <- all.vars(formula)
  data <- dplyr::select(data, index, varibale_name_in_equation)
  data$raw_order_data <- 1:nrow(data)
  raw_id <- index[1]
  colnames(data)[1] <- "id"
  index[1] <- "id"
  
  # Assuming unbalanced panel, get individuals' ID and max record number of individuals
  .N <- 0
  ID <- dplyr::select(data, index[1])
  ID_num <- data.table::setDT(ID)[,list(Count=.N),names(ID)]
  if(model == "within")
  {
    data <- data
    ID <- dplyr::select(data, index[1])
    ID_num <- data.table::setDT(ID)[,list(Count=.N),names(ID)]
  }
  
  # Judge the data size of calculation
  if (nrow(ID_num) > 1000)
  { # for test 40, real number should be 1000
    message("Dear my friend, thanks for your patience!. We pass the bandwidth\n",
            "selection part. Now, regression! This should be faster. Thanks.\n",
            ".............................................................\n")
    huge_data_size <- TRUE
  }
  else
  {
    huge_data_size <- FALSE
  }
  
  # Panel SDF preparation
  SDF@data <- dplyr::select(SDF@data, dplyr::all_of(raw_id))
  colnames(SDF@data)[1] <- "id"
  dp.locat <- sp::coordinates(SDF)
  coord <- cbind(as.data.frame(dp.locat), SDF@data$id)
  colnames(coord) <- c("X", "Y", "id")
  data <- dplyr::left_join(data, coord, by = "id")
  
  lvl1_data <- data
  
  ### transform
  lvl1_data_location <- lvl1_data %>% dplyr::select(X, Y, raw_order_data)
  lvl1_data <- lvl1_data %>% dplyr::select(-X, -Y, -raw_order_data)
  
  mean_data <- lvl1_data %>% dplyr::select(-period)
  mean_data <- mean_data %>% aggregate(by = list(mean_data$id), FUN = mean)
  mean_data <- mean_data %>% dplyr::select(-"Group.1")
  
  id_column <- lvl1_data %>% dplyr::select("id")
  mean_data <- left_join(id_column, mean_data, by = 'id')
  
  transformed_data <- lvl1_data[,3:ncol(lvl1_data)] - mean_data[,2:ncol(mean_data)]
  transformed_data <- cbind(lvl1_data[,1:2], transformed_data)
  transformed_data <- cbind(transformed_data, lvl1_data_location)
  
  return(transformed_data)
}

get_ID_list <- function(formula, data, index, SDF, bw = NULL, adaptive = FALSE, p = 2,
                        effect = "individual", model = c("pooling", "within", "random"), random.method = "swar",
                        kernel = "bisquare", longlat = FALSE)
{
  # Data preparation
  varibale_name_in_equation <- all.vars(formula)
  data <- dplyr::select(data, index, varibale_name_in_equation)
  data$raw_order_data <- 1:nrow(data)
  raw_id <- index[1]
  colnames(data)[1] <- "id"
  index[1] <- "id"
  
  # Assuming unbalanced panel, get individuals' ID and max record number of individuals
  .N <- 0
  ID <- dplyr::select(data, index[1])
  ID_num <- data.table::setDT(ID)[,list(Count=.N),names(ID)]
  if(model == "within")
  {
    data <- data
    ID <- dplyr::select(data, index[1])
    ID_num <- data.table::setDT(ID)[,list(Count=.N),names(ID)]
  }
  return(ID_num)
}

gwpr_A_normality <- function(bw, data, SDF, ID_list, formula, p, longlat, adaptive,
                             model, index, kernel = "bisquare", effect = "individual",
                             random.method = "swar")
{
  GW.arguments <- list(formula = formula, individual.number = nrow(ID_list), bw = bw,
                       kernel = kernel, adaptive = adaptive, p = p, longlat = longlat)
  message("************************ GWPR Begin *************************\n",
          "Formula: ", paste(as.character(formula)[2], " = ", as.character(formula)[3]), " -- Individuals: ", nrow(ID_list), "\n",
          "Bandwidth: ", bw, " ---- ", "Adaptive: ", adaptive, "\n",
          "Model: ", model, " ---- ", "Effect: ", effect, "\n")
  ID_list_single <- as.vector(ID_list[[1]])
  output_result <- data.frame(Doubles = double())
  y_yhat_resid <- data.frame(Doubles = double())
  loop_times <- 1
  wgt = 0
  for (ID_individual in ID_list_single)
  {
    data$aim[data$id == ID_individual] <- 1
    data$aim[data$id != ID_individual] <- 0
    subsample <- data
    subsample <- subsample[order(-subsample$aim),]
    dp_locat_subsample <- dplyr::select(subsample, 'X', 'Y')
    dp_locat_subsample <- as.matrix(dp_locat_subsample)
    dMat <- GWmodel::gw.dist(dp.locat = dp_locat_subsample, rp.locat = dp_locat_subsample,
                             focus = 1, p=p, longlat=longlat)
    subsample$dist <- as.vector(dMat)
    subsample <- subsample[order(subsample$dist),]
    id_subsample <- dplyr::select(subsample, "id")
    id_subsample <- id_subsample[!duplicated(id_subsample$id),]
    id_subsample <- as.data.frame(id_subsample)
    id_subsample <- id_subsample[1:bw,]
    id_subsample <- as.data.frame(id_subsample)
    colnames(id_subsample) <- "id"
    id_subsample <- dplyr::mutate(id_subsample, flag = 1)
    subsample <- dplyr::inner_join(subsample, id_subsample, by = "id")
    bw_of_total <- nrow(subsample)
    weight <- GWmodel::gw.weight(as.numeric(subsample$dist), bw=bw_of_total, kernel=kernel, adaptive=adaptive)
    #subsample$wgt <- as.vector(weight)
    #subsample["no2_measured_ug.m3"] <- subsample["no2_measured_ug.m3"] * subsample['wgt']
    #subsample["ug_m2_troposphere_no2"] <- subsample["ug_m2_troposphere_no2"] * subsample['wgt']
    #subsample["ter_pressure"] <- subsample["ter_pressure"] * subsample['wgt']
    #subsample["temp"] <- subsample["temp"] * subsample['wgt']
    #subsample["ndvi"] <- subsample["ndvi"] * subsample['wgt']
    #subsample["precipitation"] <- subsample["precipitation"] * subsample['wgt']
    #subsample["PBLH"] <- subsample["PBLH"] * subsample['wgt']
    
    subsample["no2_measured_ug.m3"] <- subsample["no2_measured_ug.m3"]
    subsample["ug_m2_troposphere_no2"] <- subsample["ug_m2_troposphere_no2"]
    subsample["ter_pressure"] <- subsample["ter_pressure"]
    subsample["temp"] <- subsample["temp"]
    subsample["ndvi"] <- subsample["ndvi"]
    subsample["precipitation"] <- subsample["precipitation"]
    subsample["PBLH"] <- subsample["PBLH"]
    
    #value_1 <- ifelse(shapiro.test(subsample$no2_measured_ug.m3)$p.value > 0.05, 1, 0)
    #value_2 <- ifelse(shapiro.test(subsample$ug_m2_troposphere_no2)$p.value > 0.05, 1, 0)
    #value_3 <- ifelse(shapiro.test(subsample$ter_pressure)$p.value > 0.05, 1, 0)
    #value_4 <- ifelse(shapiro.test(subsample$temp)$p.value > 0.05, 1, 0)
    #value_5 <- ifelse(shapiro.test(subsample$ndvi)$p.value > 0.05, 1, 0)
    #value_6 <- ifelse(shapiro.test(subsample$precipitation)$p.value > 0.05, 1, 0)
    #value_7 <- ifelse(shapiro.test(subsample$PBLH)$p.value > 0.05, 1, 0)
    
    value_1 <- shapiro.test(subsample$no2_measured_ug.m3)$p.value
    value_2 <- shapiro.test(subsample$ug_m2_troposphere_no2)$p.value
    value_3 <- shapiro.test(subsample$ter_pressure)$p.value
    value_4 <- shapiro.test(subsample$temp)$p.value
    value_5 <- shapiro.test(subsample$ndvi)$p.value
    value_6 <- shapiro.test(subsample$precipitation)$p.value
    value_7 <- shapiro.test(subsample$PBLH)$p.value
    
    line = c(ID_individual, value_1, value_2, value_3, value_4, value_5, value_6, value_7)
    output_result <- rbind(output_result, line)
  }
  colnames(output_result) <- c("id", all.vars(formula))
  #SDF <- sp::merge(SDF, output_result, by = "id")
  return(output_result)
}
