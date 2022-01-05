# Author: M.L.

# note: these two function are from GWPR.light v0.1.1. Since the creator did not take the pooled
#       regression without intercept into account. Therefore, we here revise the code
# output: "GWPR.user", a function used in Cross Validation, R 03_AN

# end

bw.GWPR.step.selection <- function(formula, data, index, SDF, adaptive = FALSE, p = 2,bigdata = FALSE, upperratio = 0.25,
                    effect = "individual", model = c("pooling", "within", "random"), random.method = "swar",
                    approach = c("CV","AIC"), kernel = "bisquare", longlat = FALSE, doParallel = FALSE,
                    cluster.number = 2, human.set.range = FALSE, h.upper = NULL, h.lower = NULL,
                    gradientIncrecement = FALSE, GI.step = NULL, GI.upper = NULL, GI.lower = NULL)
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
  if(!(model %in% c("pooling", "within", "random")))
  {
    stop("This version GWPR only accept \"pooling\", \"within\", or \"random\"")
  }
  if(!(approach %in% c("CV", "AIC")))
  {
    stop("This version GWPR only accept \"CV\"and \"AIC\" approach ")
  }
  
  # Data preparation
  varibale_name_in_equation <- all.vars(formula)
  data <- dplyr::select(data, index, varibale_name_in_equation)
  data$raw_order_data <- 1:nrow(data)
  raw_id <- index[1]
  colnames(data)[1] <- "id"
  index[1] <- "id"
  
  # Assuming unbalanced panel, get individuals' ID and max record number of individuals
  ID <- dplyr::select(data, "id")
  .N <- 0
  ID_num <- data.table::setDT(ID)[,list(Count = .N),names(ID)]
  if(model == "within" )
  {
    data <- drop_ID_with_single_observation(data, ID_num)
    ID <- dplyr::select(data, "id")
    ID_num <- data.table::setDT(ID)[,list(Count = .N),names(ID)]
  }
  
  # Judge the datasize of calculation
  if ((nrow(ID_num) > 1000) & !doParallel)
  {
    message("Dear my friend, more than 1,000 individuals in your dataset.\n",
            "It would be time-consuming. We use \"-\" and \"|\" to inform you where\n",
            "we are. A \"-\" is equal to 2.5% in once score calculation. A \"|\" is\n",
            "25%. Do not feel nervous or boring! Your research is on the way!\n",
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
  SDF <- 0 # drop the SDF
  colnames(coord) <- c("X", "Y", "id")
  data <- dplyr::left_join(data, coord, by = "id")
  lvl1_data <- data # data put into calculation
  
  if(human.set.range)
  {
    message("...............................................................................................\n",
            "Now, the range of bandwidth selection is set by the user\n",
            "We assume that the user is familiar with bandwidth selection, and uses this setting to reduce calculation time.\n",
            "If not, please stop the current calculation, and set \"human.set.range\" as FALSE\n",
            "...............................................................................................\n")
    if(is.null(h.lower))
    {
      stop("Need to set the lower boundary (h.lower)!")
    }
    if(is.null(h.upper))
    {
      stop("Need to set the upper boundary (h.upper)!")
    }
    if(h.lower > h.upper)
    {
      stop("h.lower should be smaller than h.upper")
    }
    lower <- h.lower
    upper <- h.upper
  }
  else
  {
    # decide upper and lower boundary
    lower.freedom <- protect_model_with_enough_freedom(formula = formula, data = lvl1_data, ID_list = ID_num,
                                                       index = index, p = p, longlat = longlat)
    message("To make sure every subsample have enough freedom, the minimum number of individuals is ",lower.freedom, "\n")
    if(adaptive)
    {
      upper <- nrow(ID_num)
      lower <- lower.freedom + 1
      if((model == "random")&(random.method == "swar"))
      {
        lower <- length(varibale_name_in_equation) + 1
      }
    }
    else
    {
      b.box <- sp::bbox(dp.locat)
      upper <- sqrt((b.box[1,2]-b.box[1,1])^2+(b.box[2,2]-b.box[2,1])^2)
      lower <- upper/5000
      if ((model == "random")&(random.method == "swar"))
      {
        lower <- protect_model_with_least_individuals(lvl1_data, ID_num, index, kernel, p, longlat,
                                                      bw_panel = (length(varibale_name_in_equation) + 1))
      }
      else
      {
        lower <- protect_model_with_least_individuals(lvl1_data, ID_num, index,
                                                      kernel, p, longlat, bw_panel = lower.freedom)
      }
    }
    
    if(bigdata)
    {
      upper <- upper * upperratio
      lower <- lower
    }
    
    if(bigdata)
    {
      message("You set the \"bigdata\" is: ", bigdata, ". The ratio is: ", upperratio, "\n",
              "Now the lower boundary of the bandwidth selection is ", lower, ", and upper boundary is ", upper,".\n",
              "Note: if the optimal bandwidth is close to the upper boundary, you need to increase the ratio.\n",
              "However, you should also know that the larger ratio requires more memory. Please, balance them and enjoy your research.\n")
    }
    if(huge_data_size)
    {
      message("Data Prepared! Go!............................................\n")
    }
  }
  message("The upper boundary is ", upper,", and the lower boundary is ", lower,"\n")
  if(doParallel)
  {
    message("..................................................................................\n")
    message("You use parallel process, so be careful about your memory usage. Cluster number: ", cluster.number,"\n")
    if(cluster.number > parallel::detectCores())
    {
      stop("The cluster number exceeds the cores you have")
    }
    else
    {
      if(cluster.number > (parallel::detectCores()-2))
      {
        warning("You might use too many cluster, only one left for other task.")
      }
    }
    #0.1.2
    if (gradientIncrecement)
    {
      if (is.null(GI.upper) | is.null(GI.lower) | is.null(GI.step))
      {
        stop("Please input upper, lower boundaries (GI.upper and GI.lower) and step length (GI.step) of GI")
      }
      if (adaptive)
      {
        message("Since GI method is used, so the GI.upper is the real upper boundary: ",
                GI.upper, " lower boundary: ", GI.lower," step length: ", GI.step)
        BandwidthVector <- c()
        ScoreVector <- c()
        if(approach == "CV")
        {
          bw.now <- GI.lower
          while (bw.now < GI.upper)
          {
            BandwidthVector <- append(BandwidthVector, bw.now)
            Score <- CV_A_para.step(bw = bw.now, data = lvl1_data, ID_list = ID_num,
                                    formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                                    model = model, index = index, kernel = kernel, effect = effect,
                                    random.method = random.method,  cluster.number = cluster.number)
            ScoreVector <- append(ScoreVector, Score)
            bw.now = bw.now + GI.step
          }
          BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
          return(BandwidthSocreTable)
        }
        else
        {
          message("AIC is not coming")
        }
      } 
      else
      {
        message("Since GI method is used, so the GI.upper is the real upper boundary: ",
                GI.upper, " lower boundary: ", GI.lower," step length: ", GI.step)
        BandwidthVector <- c()
        ScoreVector <- c()
        if(approach == "CV")
        {
          bw.now <- GI.lower
          while (bw.now < GI.upper)
          {
            BandwidthVector <- append(BandwidthVector, bw.now)
            Score <- CV_F_para.step(bw = bw.now, data = lvl1_data, ID_list = ID_num,
                                    formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                                    model = model, index = index, kernel = kernel, effect = effect,
                                    random.method = random.method,  cluster.number = cluster.number)
            ScoreVector <- append(ScoreVector, Score)
            bw.now = bw.now + GI.step
          }
          BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
          return(BandwidthSocreTable)
        }
        else
        {
          message("AIC is not coming")
        }
      }
    }
    #0.1.2 /|\
    else
    {
      if (adaptive)
      {
        if(approach == "CV")
        {
          bw <- gold(CV_A_para, xL = lower, xU = upper, adapt.bw = adaptive, data = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method,  cluster.number = cluster.number)
        }
        else
        {
          bw <- gold(AIC_A_para, xL = lower, xU = upper, adapt.bw = adaptive, data_input = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, cluster.number = cluster.number)
        }
      }
      else
      {
        if(approach == "CV")
        {
          bw <- gold(CV_F_para, xL = lower, xU = upper, adapt.bw = adaptive, data = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method,  cluster.number = cluster.number)
        }
        else
        {
          bw <- gold(AIC_F_para, xL = lower, xU = upper, adapt.bw = adaptive, data_input = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method,  cluster.number = cluster.number)
        }
      }
      return(bw)
    }
  }
  else
  {
    #0.1.2
    if (gradientIncrecement)
    {
      if (adaptive)
      {
        stop("Grandient Increcement Selection only accept Fixed bandwidth!")
      }
      if (is.null(GI.upper) | is.null(GI.lower) | is.null(GI.step))
      {
        stop("Please input upper, lower boundaries (GI.upper and GI.lower) and step length (GI.step) of GI")
      }
      message("Since GI method is used, so the GI.upper is the real upper boundary: ",
              GI.upper, " lower boundary: ", GI.lower," step length: ", GI.step)
      BandwidthVector <- c()
      ScoreVector <- c()
      if(approach == "CV")
      {
        bw.now <- GI.lower
        while (bw.now < GI.upper)
        {
          BandwidthVector <- append(BandwidthVector, bw.now)
          Score <- CV_F.step(bw = bw.now, data = lvl1_data, ID_list = ID_num,
                              formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                              model = model, index = index, kernel = kernel, effect = effect,
                              random.method = random.method, huge_data_size = huge_data_size)
          ScoreVector <- append(ScoreVector, Score)
          bw.now = bw.now + GI.step
        }
        BandwidthSocreTable <- cbind(BandwidthVector, ScoreVector)
        return(BandwidthSocreTable)
      }
      else
      {
        message("AIC is not coming")
      }
    }
    #0.1.2 /|\
    else
    {
      if (adaptive)
      {
        if(approach == "CV")
        {
          bw <- gold(CV_A, xL = lower, xU = upper, adapt.bw = adaptive, data = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
        }
        else
        {
          bw <- gold(AIC_A, xL = lower, xU = upper, adapt.bw = adaptive, data_input = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
        }
      }
      else
      {
        if(approach == "CV")
        {
          bw <- gold(CV_F, xL = lower, xU = upper, adapt.bw = adaptive, data = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
        }
        else
        {
          bw <- gold(AIC_F, xL = lower, xU = upper, adapt.bw = adaptive, data_input = lvl1_data, ID_list = ID_num,
                     formula = formula, p = p, longlat = longlat, adaptive = adaptive,
                     model = model, index = index, kernel = kernel, effect = effect,
                     random.method = random.method, huge_data_size = huge_data_size)
        }
      }
    }
  }
  #v0.1.1
  #return(bw)
}

CV_F.step <- function(bw, data, ID_list, formula, p, longlat, adaptive, kernel,
                      model = model, index = index, effect = effect,
                      random.method = random.method, huge_data_size)
{
  #  v0.1.1 the loss function is based on local r2
  #  CVscore_vector <- c()
  #v0.1.2
  residualsVector <- c()
  ID_list_single <- as.vector(ID_list[[1]])
  loop_times <- 1
  wgt <- 0
  varibale_name_in_equation <- all.vars(formula)
  for (ID_individual in ID_list_single)
  {
    data$aim[data$id == ID_individual] <- 1
    data$aim[data$id != ID_individual] <- 0
    subsample <- data
    #v0.1.2
    numberOfAim <- nrow(subsample[subsample$aim == 1,])
    subsample <- subsample[order(-subsample$aim),]
    dp_locat_subsample <- dplyr::select(subsample, 'X', 'Y')
    dp_locat_subsample <- as.matrix(dp_locat_subsample)
    dMat <- GWmodel::gw.dist(dp.locat = dp_locat_subsample, rp.locat = dp_locat_subsample,
                             focus = 1, p=p, longlat=longlat)
    weight <- GWmodel::gw.weight(as.numeric(dMat), bw=bw, kernel=kernel, adaptive=adaptive)
    subsample$wgt <- as.vector(weight)
    subsample <- subsample[(subsample$wgt > 0.01),]
    Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                   stringsAsFactors = default.stringsAsFactors())
    plm_subsample <- try(plm::plm(formula=formula, model = model, data = Psubsample,
                                  effect = effect, index = index, weights = wgt,
                                  random.method = random.method), silent = TRUE)
    #0.1.2
    if(!inherits(plm_subsample, "try-error"))
    {
      residualsLocalAim <-  plm_subsample$residuals[1:numberOfAim]
    }
    else
    {
      residualsLocalAim <- Inf
    }
    residualsVector <- append(residualsVector, residualsLocalAim)
    #    v0.1.1 the loss function is based on local r2
    #    if(!inherits(plm_subsample, "try-error"))
    #    {
    #      CVscore <- nrow(subsample) * sum(plm_subsample$residuals^2) /
    #        (nrow(subsample) - length(varibale_name_in_equation) + 1)^2
    #    }
    #    else
    #    {
    #      CVscore <- Inf
    #    }
    #    CVscore_vector <- append(CVscore_vector, CVscore)
    
    if (huge_data_size == T)
    {
      progress_bar(loop_times = loop_times, nrow(ID_list))
      loop_times <- loop_times + 1
    }
  }
  #  v0.1.1 the loss function is based on local r2
  #  mean_CVscore <- mean(CVscore_vector)
  #  cat("Fixed Bandwidth:", bw, "CV score:", mean_CVscore, "\n")
  #  return(mean_CVscore)
  #v0.1.2
  CVscore <- nrow(data) * sum(residualsVector^2) /
    (nrow(data) - length(varibale_name_in_equation) + 1)^2
  cat("Fixed Bandwidth:", bw, "CV score:", CVscore, "\n")
  return(CVscore)
}

CV_F_para.step <- function(bw, data, ID_list, formula, p, longlat, adaptive, kernel,
                            model = model, index = index, effect = effect,
                            random.method = random.method, cluster.number = cluster.number)
{
  ID_list_single <- as.vector(ID_list[[1]])
  wgt <- 0
  ID_individual <- 0
  varibale_name_in_equation <- all.vars(formula)
  cl <- parallel::makeCluster(cluster.number)
  doParallel::registerDoParallel(cl)
  #  v0.1.1 the loss function is based on local r2
  # CVscore_vector <- foreach(ID_individual = ID_list_single, .combine = c) %dopar%
  # v0.1.2
  residualsVector <- foreach(ID_individual = ID_list_single, .combine = c) %dopar%
    {
      data$aim[data$id == ID_individual] <- 1
      data$aim[data$id != ID_individual] <- 0
      subsample <- data
      #v0.1.2
      numberOfAim <- nrow(subsample[subsample$aim == 1,])
      subsample <- subsample[order(-subsample$aim),]
      dp_locat_subsample <- dplyr::select(subsample, 'X', 'Y')
      dp_locat_subsample <- as.matrix(dp_locat_subsample)
      dMat <- GWmodel::gw.dist(dp.locat = dp_locat_subsample, rp.locat = dp_locat_subsample,
                               focus = 1, p=p, longlat=longlat)
      weight <- GWmodel::gw.weight(as.numeric(dMat), bw=bw, kernel=kernel, adaptive=adaptive)
      subsample$wgt <- as.vector(weight)
      subsample <- subsample[(subsample$wgt > 0.01),]
      Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                     stringsAsFactors = default.stringsAsFactors())
      plm_subsample <- try(plm::plm(formula=formula, model=model, data=Psubsample,
                                    effect = effect, index=index, weights = wgt,
                                    random.method = random.method), silent = TRUE)
      # v0.1.1
      #    if(!inherits(plm_subsample, "try-error"))
      #    {
      #      CVscore <- nrow(subsample) * sum(plm_subsample$residuals^2) /
      #        (nrow(subsample) - length(varibale_name_in_equation) + 1)^2
      #    }
      #    else
      #    {
      #      CVscore <- Inf
      #    }
      #0.1.2
      if(!inherits(plm_subsample, "try-error"))
      {
        residualsLocalAim <-  plm_subsample$residuals[1:numberOfAim]
      }
      else
      {
        residualsLocalAim <- Inf
      }
    }
  parallel::stopCluster(cl)
  #  v0.1.1 the loss function is based on local r2
  #  mean_CVscore <- mean(CVscore_vector)
  #  cat("Fixed Bandwidth:", bw, "CV score:", mean_CVscore, "\n")
  #  return(mean_CVscore)
  #v0.1.2
  CVscore <- nrow(data) * sum(residualsVector^2) /
    (nrow(data) - length(varibale_name_in_equation) + 1)^2
  cat("Fixed Bandwidth:", bw, "CV score:", CVscore, "\n")
  return(CVscore)
}

drop_ID_with_single_observation <- function(data, ID_num)
{
  data <- dplyr::left_join(data, ID_num, by = "id")
  data <- data[(data$Count != 1),]
  data <- dplyr::select(data, -"Count")
  return(data)
}

protect_model_with_enough_freedom <- function(formula, data, ID_list, index,
                                              p, longlat)
{
  ID_list_single <- as.vector(ID_list[[1]])
  step_increase_lower <- 1
  lower <- 1
  go_out <- T
  required_freedom <- length(all.vars(formula))
  while((step_increase_lower < 1001)&(go_out == T))
  {
    step_increase_lower <- step_increase_lower + 1 #required at least two individuals
    lower <- lower + 1
    for (ID_individual in ID_list_single)
    {
      go_out <- F
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
      id_subsample <- as.data.frame(id_subsample) #TestCode
      id_subsample <- id_subsample[1:lower,]
      id_subsample <- as.data.frame(id_subsample)
      colnames(id_subsample) <- "id"
      id_subsample <- dplyr::mutate(id_subsample, flag = 1)
      subsample <- dplyr::inner_join(subsample, id_subsample, by = "id")
      if(nrow(subsample) < required_freedom)
      {
        go_out <- T
        break
      }
    }
  }
  return(lower)
}

protect_model_with_least_individuals <- function(data, ID_list, index,
                                                 kernel, p, longlat, bw_panel)
{
  ID_list_single <- as.vector(ID_list[[1]])
  max_dist <- c()
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
    id_subsample <- id_subsample[1:bw_panel,]
    id_subsample <- as.data.frame(id_subsample)
    colnames(id_subsample) <- "id"
    id_subsample <- dplyr::mutate(id_subsample, flag = 1)
    subsample <- dplyr::inner_join(subsample, id_subsample, by = "id")
    max_dist <- append(max_dist ,max(subsample$dist))
  }
  lower <- max(max_dist) * 1.011 # because individuals with the weight lower than 0.01 would be ignored,
  # to guarantee all the individuals used in local panel model, we use 1.011 here.
  return(lower)
}

CV_A_para.step <- function(bw, data, ID_list, formula, p, longlat, adaptive, kernel,
                           model = model, index = index, effect = effect,
                           random.method = random.method, cluster.number = cluster.number)
{
  ID_list_single <- as.vector(ID_list[[1]])
  wgt <- 0
  ID_individual <- 0
  varibale_name_in_equation <- all.vars(formula)
  cl <- parallel::makeCluster(cluster.number)
  doParallel::registerDoParallel(cl)
  #  v0.1.1 the loss function is based on local r2
  # CVscore_vector <- foreach(ID_individual = ID_list_single, .combine = c) %dopar%
  # v0.1.2
  residualsVector <- foreach(ID_individual = ID_list_single, .combine = c) %dopar%
    {
      data$aim[data$id == ID_individual] <- 1
      data$aim[data$id != ID_individual] <- 0
      subsample <- data
      #v0.1.2
      numberOfAim <- nrow(subsample[subsample$aim == 1,])
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
      id_subsample$flag <- 1
      subsample <- dplyr::inner_join(subsample, id_subsample, by = "id")
      bw_to_total <- nrow(subsample)
      weight <- GWmodel::gw.weight(as.numeric(subsample$dist), bw=bw_to_total, kernel=kernel, adaptive=adaptive)
      subsample$wgt <- as.vector(weight)
      Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                     stringsAsFactors = default.stringsAsFactors())
      plm_subsample <- plm::plm(formula=formula, model=model, data=Psubsample,
                                effect = effect, index=index, weights = wgt,
                                random.method = random.method)
      # v0.1.1
      #    if(!inherits(plm_subsample, "try-error"))
      #    {
      #      CVscore <- nrow(subsample) * sum(plm_subsample$residuals^2) /
      #        (nrow(subsample) - length(varibale_name_in_equation) + 1)^2
      #    }
      #    else
      #    {
      #      CVscore <- Inf
      #    }
      #0.1.2
      if(!inherits(plm_subsample, "try-error"))
      {
        residualsLocalAim <-  plm_subsample$residuals[1:numberOfAim]
      }
      else
      {
        residualsLocalAim <- Inf
      }
    }
  parallel::stopCluster(cl)
  #  v0.1.1 the loss function is based on local r2
  #  mean_CVscore <- mean(CVscore_vector)
  #  cat("Fixed Bandwidth:", bw, "CV score:", mean_CVscore, "\n")
  #  return(mean_CVscore)
  #v0.1.2
  CVscore <- nrow(data) * sum(residualsVector^2) /
    (nrow(data) - length(varibale_name_in_equation) + 1)^2
  cat("Fixed Bandwidth:", bw, "CV score:", CVscore, "\n")
  return(CVscore)
}