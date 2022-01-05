# Author: M.L.

# note: these two function are from GWPR.light v0.1.1. Since the creator did not take the pooled
#       regression without intercept into account. Therefore, we here revise the code
# output: "GWPR.user", a function used in Cross Validation, R 06_AN

# end

GWPR.user <- function(formula, data, index, SDF, bw = NULL, adaptive = FALSE, p = 2,
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
    data <- drop_ID_with_single_observation(data, ID_num)
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
  
  if(huge_data_size)
  {
    message("Data Prepared! Go!............................................\n")
  }
  
  # GWPRegression
  if (adaptive)
  {
    result <- gwpr_A.user(bw = bw, data = lvl1_data, SDF, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
  }
  else
  {
    result <- gwpr_F.user(bw = bw, data = lvl1_data, SDF, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
  }
  message("The R2 is: ", result$R2,"\n")
  message("Note: in order to avoid mistakes, we forced a rename of the individuals'ID as \"id\". \n")
  return(result)
}

gwpr_F.user <- function(bw = bw, data, SDF, ID_list,
                   formula = formula, p = p, longlat = longlat, adaptive = F,
                   model = model, index = index, kernel = kernel, effect = effect,
                   random.method = random.method, huge_data_size = huge_data_size)
{
  GW.arguments <- list(formula = formula, individual.number = nrow(ID_list), bw = bw,
                       kernel = kernel, adaptive = adaptive, p = p, longlat = longlat)
  message("************************ GWPR Begin *************************\n",
          "Formula: ", paste(as.character(formula)[2], " = ", as.character(formula)[3]), " -- Individuals: ", nrow(ID_list), "\n",
          "Bandwidth: ", bw, " ---- ", "Adaptive: ", adaptive, "\n",
          "Model: ", model, " ---- ", "Effect: ", effect, "\n")
  global_plm <- plm::plm(formula=formula, model=model, data = data,
                         effect = effect, index=index, random.method = random.method)
  ID_list_single <- as.vector(ID_list[[1]])
  output_result <- data.frame(Doubles = double())
  y_yhat_resid <- data.frame(Doubles = double())
  loop_times <- 1
  wgt <- 0
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
    weight <- GWmodel::gw.weight(as.numeric(dMat), bw=bw, kernel=kernel, adaptive=adaptive)
    subsample$wgt <- as.vector(weight)
    subsample <- subsample[(subsample$wgt > 0),]
    Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                   stringsAsFactors = default.stringsAsFactors())
    plm_subsample <- plm::plm(formula=formula, model=model, data=Psubsample,
                              effect = effect, index=index, weights = wgt,
                              random.method = random.method)
    coefMat <- lmtest::coeftest(plm_subsample)
    local_r2 <- plm::r.squared(plm_subsample)
    result_line <- c(ID_individual, coefMat[,1], coefMat[,2], coefMat[,3], local_r2)
    output_result <- rbind(output_result, result_line)
    dataset_add_resid <- cbind(Psubsample, plm_subsample$residuals)
    dataset_add_resid <- as.data.frame(dataset_add_resid)
    varibale_name_in_equation <- all.vars(formula)
    dataset_add_resid <- dplyr::select(dataset_add_resid, dplyr::all_of(index), dplyr::all_of(varibale_name_in_equation)[1],
                                       "plm_subsample$residuals")
    colnames(dataset_add_resid) <- c(index, "y", "resid")
    dataset_add_resid$yhat <- dataset_add_resid$y - dataset_add_resid$resid
    dataset_add_resid <- dplyr::filter(dataset_add_resid, id == ID_individual)
    y_yhat_resid <- rbind(y_yhat_resid, dataset_add_resid)
    if (huge_data_size == T)
    {
      progress_bar(loop_times = loop_times, nrow(ID_list))
      loop_times <- loop_times + 1
    }
  }
  varibale_name_in_equation <- all.vars(formula)
  if (model == "pooling")
  {
    varibale_name_in_equation_out <- varibale_name_in_equation[2:length(varibale_name_in_equation)]
  }
  else
  {
    varibale_name_in_equation_out <- varibale_name_in_equation
    varibale_name_in_equation_out[1] <- "Intercept"
  }
  colnames(output_result) <- c("id", varibale_name_in_equation_out, paste0(varibale_name_in_equation_out,"_SE"),
                               paste0(varibale_name_in_equation_out,"_TVa"), "Local_R2")
  SDF <- sp::merge(SDF, output_result, by = "id")
  y_yhat_resid[,1] <- as.numeric(as.character(y_yhat_resid[,1]))
  y_yhat_resid[,2] <- as.numeric(as.character(y_yhat_resid[,2]))
  r2 <- 1 - sum(y_yhat_resid$resid^2)/(sum((y_yhat_resid$y - mean(y_yhat_resid$y))^2))
  result_list <- list(GW.arguments = GW.arguments, R2 = r2, index = index, plm.result = global_plm,
                      raw.data = data, GWPR.residuals = y_yhat_resid, SDF = SDF)
  return(result_list)
}

gwpr_A.user <- function(bw, data, SDF, ID_list, formula, p, longlat, adaptive,
                   model, index, kernel = "bisquare", effect = "individual",
                   random.method = "swar", huge_data_size = huge_data_size)
{
  GW.arguments <- list(formula = formula, individual.number = nrow(ID_list), bw = bw,
                       kernel = kernel, adaptive = adaptive, p = p, longlat = longlat)
  message("************************ GWPR Begin *************************\n",
          "Formula: ", paste(as.character(formula)[2], " = ", as.character(formula)[3]), " -- Individuals: ", nrow(ID_list), "\n",
          "Bandwidth: ", bw, " ---- ", "Adaptive: ", adaptive, "\n",
          "Model: ", model, " ---- ", "Effect: ", effect, "\n")
  global_plm <- plm::plm(formula=formula, model=model, data=data,
                         effect = effect, index=index, random.method = random.method)
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
    subsample$wgt <- as.vector(weight)
    Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                   stringsAsFactors = default.stringsAsFactors())
    plm_subsample <- plm::plm(formula=formula, model=model, data=Psubsample,
                              effect = effect, index=index, weights = wgt,
                              random.method = random.method)
    coefMat <- lmtest::coeftest(plm_subsample)
    local_r2 <- plm::r.squared(plm_subsample)
    result_line <- c(ID_individual, coefMat[,1], coefMat[,2], coefMat[,3], local_r2)
    output_result <- rbind(output_result, result_line)
    dataset_add_resid <- cbind(Psubsample, plm_subsample$residuals)
    dataset_add_resid <- as.data.frame(dataset_add_resid)
    varibale_name_in_equation <- all.vars(formula)
    dataset_add_resid <- dplyr::select(dataset_add_resid, dplyr::all_of(index), dplyr::all_of(varibale_name_in_equation)[1],
                                       "plm_subsample$residuals")
    colnames(dataset_add_resid) <- c(index, "y", "resid")
    dataset_add_resid$yhat <- dataset_add_resid$y - dataset_add_resid$resid
    dataset_add_resid <- dplyr::filter(dataset_add_resid, id == ID_individual)
    y_yhat_resid <- rbind(y_yhat_resid, dataset_add_resid)
    if (huge_data_size == T)
    {
      progress_bar(loop_times = loop_times, nrow(ID_list))
      loop_times <- loop_times + 1
    }
  }
  varibale_name_in_equation <- all.vars(formula)
  if (model == "pooling")
  {
    varibale_name_in_equation_out <- varibale_name_in_equation[2:length(varibale_name_in_equation)]
  }
  else
  {
    varibale_name_in_equation_out <- varibale_name_in_equation
    varibale_name_in_equation_out[1] <- "Intercept"
  }
  colnames(output_result) <- c("id", varibale_name_in_equation_out, paste0(varibale_name_in_equation_out,"_SE"),
                               paste0(varibale_name_in_equation_out,"_TVa"), "Local_R2")
  SDF <- sp::merge(SDF, output_result, by = "id")
  y_yhat_resid[,1] <- as.numeric(as.character(y_yhat_resid[,1]))
  y_yhat_resid[,2] <- as.numeric(as.character(y_yhat_resid[,2]))
  r2 <- 1 - sum(y_yhat_resid$resid^2)/(sum((y_yhat_resid$y - mean(y_yhat_resid$y))^2))
  result_list <- list(GW.arguments = GW.arguments, R2 = r2, index = index, plm.result = global_plm,
                      raw.data = data, GWPR.residuals = y_yhat_resid, SDF = SDF)
  return(result_list)
}
