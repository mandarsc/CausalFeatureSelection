# Causal Inference Packages
library(pcalg)
library(graph)

# Local Causal Discovery Algorithms
library(bnlearn)

# Feature Selection
library(FSelector)
library(glmnet)

# Data processing
library(pracma)
library(binhf)
library(XML)
library(tools)
library(permute)

# Multicollinearity
library(pls)

# Classifier libraries
library(C50)
library(randomForest)

localCausalFS <- function(X_train, alg, response_idx){
    col_names <- colnames(X_train)
    blacklist <- matrix(0, nrow=(ncol(X_train)-1), ncol=2)
    blacklist[, 1] <- col_names[response_idx]
    blacklist[, 2] <- col_names[-response_idx]

    # cat('Local Feature Selection: ', alg, '\n')
    if(alg=='MMMB'){
      obj <- mmpc(X_train, blacklist=blacklist, alpha=0.05, optimized=TRUE, test='zf')
      selectedFeatures <- mmpc.features <- getFeatures(obj, response_idx)
    }else if(alg=='HITON-MB'){
      obj <- si.hiton.pc(X_train, blacklist=blacklist, alpha=0.05, optimized=TRUE, test='zf')
      selectedFeatures <- hiton.features <- getFeatures(obj, response_idx)
    }else if(alg=='ARACNE'){
        obj <- aracne(X_train, blacklist=blacklist, mi='mi-g')
        selectedFeatures <- aracne.features <- getFeatures(obj, response_idx)
    }else if(alg=='MMHC'){
        obj <- mmhc(X_train, blacklist=blacklist)
        selectedFeatures <- mmhc.features <- getFeatures(obj, response_idx)
    }
    return(selectedFeatures)
}

univariateFS <- function(X_train, y_idx, n_features){
    y_var <- colnames(X_train)[y_idx]
    x_var <- colnames(X_train)[-y_idx]
    formula <- as.formula(paste(y_var, paste(x_var, collapse='+'), sep='~'))
    spearman_corr <- rank.correlation(formula, data=X_train)
    pearson_corr <- linear.correlation(formula, data=X_train)
    info_gain <- information.gain(formula, data=X_train)

    spearman_features <- cutoff.k(spearman_corr, n_features)
    pearson_features <- cutoff.k(pearson_corr, n_features)
    info_gain_features <- cutoff.k(info_gain, n_features)
    return(list(spearman_features, pearson_features, info_gain_features))
}

regressionFS <- function(X_train, y_idx){
    y_var <- colnames(X_train)[y_idx]
    x_var <- colnames(X_train)[-y_idx]
    formula <- as.formula(paste(y_var, paste(x_var, collapse='+'), sep='~'))

    lasso_model <- glmnet(as.matrix(X_train[, -ncol(X_train)]), X_train[, y_idx], alpha=1) #alpha=1 is the lasso penalty
    lasso_features <- (coef(lasso_model, s=0.01) != 0)@i               #coeffs for min lambda for Lasso --> 0.01

    formula_null <- as.formula(paste(y_var, 1, sep='~'))
    null = lm(formula_null, data=X_train)
    full = lm(formula, data=X_train)
    
    step_fit <- step(null, scope=list(lower=null, upper=full), direction="both", trace=0, showWarnings = FALSE)
    stepr_features <- names((step_fit$coefficients)[2:length(step_fit$coefficients)])
    
    return(list(lasso_features, stepr_features))
}

oneRFS <- function(X_train, y_idx, n_features){
    y_var <- colnames(X_train)[y_idx]
    x_var <- colnames(X_train)[-y_idx]
    formula <- as.formula(paste(y_var, paste(x_var, collapse='+'), sep='~'))

    oneR_fit <- oneR(formula, data=X_train)
    oneR_features <- cutoff.k(oneR_fit, n_features)
    return(list(oneR_features))
}

mmhcObj <- function(X_train, y_idx){
  col_names <- colnames(X_train)
  blacklist <- matrix(0, nrow=(ncol(X_train)-1), ncol=2)
  blacklist[, 1] <- col_names[y_idx]
  blacklist[, 2] <- col_names[-y_idx]
  obj <- mmhc(X_train, blacklist=blacklist)
  return(as.graphNEL(obj))     
}

calCI <- function(y_pred, y_true=NULL, n_boot=1000){
  if(sum(is.na(y_pred)) > 0){
    return(NA)
  }
  n <- length(y_pred)
  if(is.null(y_true)){
    boot_sample <- sample(y_pred, n*n_boot, replace=TRUE)
    boot_sample_mat <- matrix(boot_sample, nrow=n, ncol=n_boot)  
    xbar <- mean(y_pred)
    boot_means <- colMeans(boot_sample_mat)
    delta <- boot_means - xbar
    d <- quantile(delta, c(0.05, 0.95))
    ci <- xbar - c(d[2], d[1])    
  }else{
    boot_sample <- sample(1:n, n*n_boot, replace=TRUE)
    boot_sample_mat <- matrix(boot_sample, nrow=n, ncol=n_boot)    
    x_rmse <- sqrt(mean((y_true-y_pred)^2))
    boot_rmse <- apply(boot_sample_mat, 2, function(x) sqrt(mean((y_true[x] - y_pred[x])^2)))
    delta <- boot_rmse - x_rmse
    d <- quantile(delta, c(0.05, 0.95))
    ci <- x_rmse - c(d[2], d[1])    
  }
  return(ci)
}

getFeatures <- function(obj, response_idx){
  features <- union(obj$nodes[[response_idx]]$mb, obj$nodes[[response_idx]]$nbr)
  return(features)
}

preprocessData <- function(i, orgClimData, numYears, rainfall_cols){
  idx <- numYears[-i]
  sampleClimData <- as.matrix(orgClimData[-i, -which(names(orgClimData) %in% rainfall_cols)])
  testClimData <- as.matrix(orgClimData[i, -which(names(orgClimData) %in% rainfall_cols)])
  rainfallSeasonData <- as.matrix(orgClimData[-i, rainfall_cols])

  for(k in 1:ncol(sampleClimData)){
    missingVal <- which(sampleClimData[, k] <= -99)
    if(length(missingVal) > 0){
      # print(k)
      sampleClimData[missingVal, k] <- mean(sampleClimData[-missingVal, k])
    }
    if(length(which(testClimData[, k] <= -99)) > 0)
      testClimData[, k] <- mean(sampleClimData[, k])
  }

  detrend_sampleClim <- detrendData(i, idx, sampleClimData, testClimData)
  sampleClimData <- detrend_sampleClim[[1]]
  testClimData <- detrend_sampleClim[[2]]

  detrend_rainfall <- detrendData(i, idx, rainfallSeasonData)
  rainfallSeasonData <- detrend_rainfall[[1]]

  idx.stats <- list()
  idx.stats[[1]] <- apply(sampleClimData, 2, mean)
  idx.stats[[2]] <- apply(sampleClimData, 2, sd)
    
  sampleClimData <- normalizeData(sampleClimData, 'TRAIN')  # Normalize data after detrending
  testClimData <- normalizeData(testClimData, 'TEST', idx.stats)  # Get the year left out to build test data and normalize it 

  rainfallData <- matrix(0, nrow=nrow(rainfallSeasonData), ncol=1)  
  colnames(rainfallData) <- "q_RIBFLV"
  
  for(rain_idx in 1:ncol(rainfallSeasonData)){  # Aggregate the detrended values of rainfall season into one column
    rainfallData <- rainfallData + rainfallSeasonData[, rain_idx]
  }

  rainfallNormData <- normalizeData(rainfallData, 'TRAIN')    # Normalized the aggregated rainfall values

  return(list(sampleClimData, testClimData, rainfallNormData))
}

discVariable <- function(data){
  
  discVector <- matrix(0, nrow=nrow(data), ncol=1)
  # Find the HIGH anomalous years for Var1_Month[i]
  sortedVector <- sort(data, decreasing=TRUE,index.return=TRUE)
  # Assign 1 for HIGH anomaly   
  highIdx <- which(sortedVector$x > quantile(sortedVector$x,0.67))
  discVector[sortedVector$ix[highIdx]] <- 1
  # Assign 3 for LOW anomaly   
  lowIdx <- which(sortedVector$x < quantile(sortedVector$x,0.33))
  discVector[sortedVector$ix[lowIdx]] <- 3
  # Assign 2 for NORMAL phase 
  normalIdx <-setdiff(1:length(sortedVector$x),c(highIdx, lowIdx))
  discVector[sortedVector$ix[normalIdx]] <- 2
  return(discVector)
}

extract <- function(string){
  a=unlist(strsplit(string," "))
}

tetradToR <- function(filename){
    data = xmlParse(filename)
    xml_data = xmlToList(data)
    gr = new("graphNEL", nodes=unlist(xml_data$variable), edgemode = "directed")

    import_edges <- matrix(unlist(lapply(xml_data$edges, FUN=extract)),ncol=3,byrow =TRUE)

    m_edge_processed = matrix(ncol=2)
    dir.edges.idx <- which(import_edges[,2]=='-->', arr.ind=TRUE)
    bidir.wgts <- rep(1, length(dir.edges.idx))
    gr <- addEdge(import_edges[dir.edges.idx, 1] , import_edges[dir.edges.idx, 3], gr, bidir.wgts)

    undir.edges.idx <- which(import_edges[,2]=='---')
    undir.wgts <- rep(1, length(undir.edges.idx))
    gr <- addEdge(import_edges[undir.edges.idx, 1] , import_edges[undir.edges.idx, 3], gr, undir.wgts)
    gr <- addEdge(import_edges[undir.edges.idx, 3] , import_edges[undir.edges.idx, 1], gr, undir.wgts)

    bidir.edges.idx <- which(import_edges[,2]=='<->')
    bidir.wgts <- rep(2, length(bidir.edges.idx))
    gr <- addEdge(import_edges[bidir.edges.idx, 1] , import_edges[bidir.edges.idx, 3], gr, bidir.wgts)
    gr <- addEdge(import_edges[bidir.edges.idx, 3] , import_edges[bidir.edges.idx, 1], gr, bidir.wgts)

    return(gr)
}

detrendData <- function(rowNo, idx, trainData, testData=NULL){
  for(k in 1:ncol(trainData)){
      tempData <- cbind(trainData[, k], idx)
      regModel <- lm(V1~., data.frame(tempData))
      coeff <- coefficients(regModel)
      for(l in 1:length(idx))
          trainData[l, k] <- (trainData[l, k] - (1*coeff[1] + idx[l]*coeff[2]))
    if(!is.null(testData)) 
      testData[, k] <- (testData[, k] - (1*coeff[1] + rowNo*coeff[2]))
  }
  return(list(trainData, testData))
}

normalizeData <- function(data, flag, stats=0){
  if(flag=='TRAIN' || flag=='ALL'){
    for(i in 1:ncol(data)){
      data[, i] <- (data[, i] - mean(data[, i]))/sd(data[, i])
    }
  }else{
    for(i in 1:ncol(data)){
      data[, i] <- (data[, i] - stats[[1]][i])/(stats[[2]][i])  
    }
  }
  return(data)
}

getDiscretizedResponse <- function(data, cols){
    y_data <- data[, cols]
    # Detrend and aggregate
    y_data <- apply(y_data, 2, detrend)
    y_data <- as.matrix(apply(y_data, 1, sum))
  
    # Normalize and discretize rainfall
    y_norm_data <- normalizeData(y_data, 'TRAIN')
    y_disc_data <- discVariable(y_norm_data)
    return(y_disc_data)
}

buildModel <- function(X_train, X_test, label, response, y_var){
  x_var <- colnames(X_train)
  formula <- as.formula(paste(y_var, paste(x_var, collapse='+'), sep='~'))
  label <- as.factor(label)
  build_c50 <- C5.0(x=X_train, y=label, method="class", trials=1)
  predict_c50 <- predict(build_c50, data.frame(X_test, check.names=FALSE))

  build_rf <- randomForest(x=X_train, y=label)
  predict_rf <- predict(build_rf, data.frame(X_test))

  build_lm <- lm(formula, data=cbind(X_train, response))
  predict_lm <- predict(build_lm, data.frame(X_test))

  # rpart_fit <- rpart(formula, data=cbind(X_train, response))
  # predict_rpart <- predict(rpart_fit, data.frame(X_test))
  return(list(build_c50, predict_c50, build_lm, predict_lm))
}

computePValue2 <- function(idx_name, beta, data, pred_idx, response_idx, flag){
    N <- nrow(data)
    response_vector <- data[, response_idx]
    # Initialize count to 1 for the original estimated causal effect
    count <- 1
    set.seed(1234201)
    alpha <- 0.05
    rand_perm <- 1000
    beta_hat_vec <- vector(length = rand_perm)
    randomized_row_set <- as.matrix(shuffleSet(n=nrow(data), nset=rand_perm))
    y_var <- colnames(data)[response_idx]
    x_var <- colnames(data)[pred_idx]
    formula <- as.formula(paste(y_var, paste(x_var, collapse='+'), sep='~'))
    beta_hat_vec <- unlist(lapply(1:rand_perm, function(x){ 
                      rand_idx <- randomized_row_set[x, ]
                      data[, response_idx] <- data[rand_idx, response_idx]
                      reg_model <- pcr(formula, data=data[, c(pred_idx, response_idx)])
                      max_var_comp <-  which.max(explvar(reg_model))
                      beta_hat <- data.frame(reg_model$coefficients)[idx_name, max_var_comp]
                      data[, response_idx] <- response_vector
                      beta_hat
                      }
                    ))

    for(i in 1:rand_perm){
      rand_idx <- randomized_row_set[i, ]
      data[, response_idx] <- data[rand_idx, response_idx]
      if(flag==1)
        reg_model <- pcr(formula, data=data[, c(pred_idx, response_idx)])
      else
        reg_model <- plsr(formula., data=data[, c(pred_idx, response_idx)], method="oscorespls")
      max_var_comp <-  which.max(explvar(reg_model))
      beta_hat <- data.frame(reg_model$coefficients)[idx_name, max_var_comp]
      beta_hat_vec[i] = beta_hat
      data[, response_idx] <- response_vector
    }
    # if(beta > 0){
    #   p_value <- (sum(beta_hat_vec > beta) + 1)/(rand_perm + 1)
    # }else if(beta < 0){
    #   p_value <- (sum(beta_hat_vec < beta) + 1)/(rand_perm + 1)
    # }else{
    #   p_value <- 1
    # }
    p_value <- (sum(abs(beta_hat_vec) > abs(beta)) + 1)/(rand_perm + 1)

    print(cat('Index:', idx_name, ' p_value: ', p_value))
    result1 <- FALSE
    if(p_value <= alpha){
      return(c(beta, p_value))
    }
    return(NaN)
}

computePValue <- function(idx_name, beta, data, pred_idx, response_idx, flag){
    N <- nrow(data)
    response_vector <- data[, response_idx]
    # Initialize count to 1 for the original estimated causal effect
    count <- 1
    set.seed(1234201)
    alpha <- 0.05
    rand_perm <- 5000
    randomized_row_set <- as.matrix(shuffleSet(n=nrow(data), nset=rand_perm))
    y_var <- colnames(data)[response_idx]
    x_var <- colnames(data)[pred_idx]
    formula <- as.formula(paste(y_var, paste(x_var, collapse='+'), sep='~'))
    beta_hat_vec <- unlist(lapply(1:rand_perm, function(x){ 
                      rand_idx <- randomized_row_set[x, ]
                      data[, response_idx] <- data[rand_idx, response_idx]
                      reg_model <- pcr(formula, data=data[, c(pred_idx, response_idx)])
                      max_var_comp <-  which.max(explvar(reg_model))
                      beta_hat <- data.frame(reg_model$coefficients)[idx_name, max_var_comp]
                      data[, response_idx] <- response_vector
                      beta_hat
                      }
                    ))
    p_value <- (sum(abs(beta_hat_vec) > abs(beta)) + 1)/(rand_perm + 1)
    if(p_value <= alpha){
      return(c(beta, p_value))
    }
    return(c(beta, p_value))
}

cfs <- function(train_data, test_data, pc_graph, dir_edge_rank, response_pos, option){
    n <- nrow(train_data)
    p <- ncol(train_data)
    col_names <- colnames(train_data)
    direct_idx <- array(0, dim=n)
    bool_vect <- array(FALSE, dim=p)
    causal_effects <- matrix(0, nrow=p, ncol=2, dimnames=list(col_names, c('Causal_Effects', 'p-value')))
    count <- 0
    remove_idx <- na_idx <- c()
    names_edges <- names_direct_edges <- c()  
    remove_idx <- na_idx <- direct_edge_train <- direct_edge_idx <- direct_edge_test <- c()
    all_effects <- data.frame()

    for(j in 1:nrow(dir_edge_rank)){

        from <- dir_edge_rank[j, 1]
        to <- dir_edge_rank[j, 2]
        
        if(from==response_pos || from %in% na_idx || to %in% na_idx){
          # print(sprintf('%d has response as cause', i))
          count <- count + 1
          remove_idx <- union(remove_idx, j)
          dir_edge_rank[j, 3] <- -1  
          next
        }

        # Identify directed edges which have arrow head at rainfall
        if(to == response_pos){
          if(!bool_vect[from]){#check if the causal effect is already calculated or not.
            bool_vect[from] <- TRUE
            from_effects <- computeCausalEffects(from, train_data, pc_graph, response_pos)  # Use all years to get the causal effect
            # from_effects <- getCausalEffects.plsr.2(from, train_data, pc_graph, response_pos)
            if(length(from_effects)==1 && is.na(from_effects)){
              na_idx <- union(na_idx, from)
              remove_idx <- union(remove_idx, j)
              dir_edge_rank[j, 3] <- -1
              next
            }
            # causal_effects[from, ] <- from_effects[1, ]
            all_effects <- rbind(all_effects, cbind(from_effects, rep(colnames(train_data)[from])))
          }
          next
        }

        # For each directed edge A->B, get the causal effect of indices in a directed causal edge on rainfall
        if(!bool_vect[from]){
            # print(cat('From: ', colnames(train_data)[from]))
            bool_vect[from] <- TRUE
            from_effects <- computeCausalEffects(from, train_data, pc_graph, response_pos)
            # from_effects <- getCausalEffects.plsr.2(from, train_data, pc_graph, response_pos)

            if(length(from_effects)==1 && is.na(from_effects)){
              causal_effects[from, ] <- 0
              na_idx <- union(na_idx, from)
              remove_idx <- union(remove_idx, j)
            }else{
              # causal_effects[from, ] <- from_effects[1,]
              all_effects <- rbind(all_effects, cbind(from_effects, rep(colnames(train_data)[from])))              
            }
        }

        if(!bool_vect[to]){
            # print(cat('To: ',colnames(train_data)[to]))
            bool_vect[to] <- TRUE
            to_effects <- computeCausalEffects(to, train_data, pc_graph, response_pos)
            # to_effects <- getCausalEffects.plsr.2(to, train_data, pc_graph, response_pos)

            if(length(to_effects)==1 && is.na(to_effects)){
              causal_effects[to, ] <- 0
              na_idx <- union(na_idx, to)
              remove_idx <- union(remove_idx, j)
            }else{
              # causal_effects[to, ] <- to_effects[1,]
              all_effects <- rbind(all_effects, cbind(to_effects, rep(colnames(train_data)[to])))
            }
        }
        # Sum the absolute values of causal effects to rank the edges
        # if(sum(abs(causal_effects[from, 1])) > 0 && sum(abs(causal_effects[to, 1])) > 0){
        #   dir_edge_rank[j, 3] <- sum(abs(causal_effects[from, 1]), abs(causal_effects[to, 1]))
        # }else{
        #   dir_edge_rank[j, 3] <- -1 
        # }
    }
    # t1 <- causal_effects[causal_effects[, 2]!=0, ]
    p_val_adj <- p.adjust(as.numeric(as.character(all_effects[,2])), method="BH")
    causal_features <- unique(as.character(all_effects[which(p_val_adj < 0.05), 3]))
    # cluster_features <- clusterFeatures(all_edges, train_data, test_data, causal_effects)
    # return(list(cluster_features[[1]], cluster_features[[2]])) 
    return(list(train_data[, causal_features], test_data[, causal_features]))
}

clusterFeatures <- function(all_edges, train_data, test_data, causal_effects){
    row.names(causal_effects) <- colnames(train_data)
    non_zero_features <- which(causal_effects[,1] != 0)
    causal_effects <- causal_effects[non_zero_features,]

    ratios <- array(0, dim=(length(causal_effects[,1])-2))
    clusters <- list()
    clusters[[1]] <- 'temp'
    temp <- as.matrix(causal_effects[, 1])
    for(i in 2:(length(causal_effects[, 1])-1)){
        # print(cat('row num ', length(causal_effects[, 1]), ' i ', i))
        clusters[[i]] <- kmeans(x = temp, centers = i)
        ratios[i] <- 1 - clusters[[i]]$tot.withinss / clusters[[i]]$totss
    }
    opt_k <- find.maximum.distance.point(ratios[-1])+1
    feature_names <- paste('Feature', 1:opt_k, sep='_')
    combined.newFeature.train <- matrix(0, nrow=nrow(train_data), ncol=opt_k, dimnames=list(c(), feature_names))
    combined.newFeature.test <- matrix(0, nrow=nrow(test_data), ncol=opt_k, dimnames=list(c(), feature_names))

    myclust <- data.frame(clusters[[opt_k]]$cluster)
    selected_features <- c()
    for(i in 1:opt_k){
        features <- which(myclust == i)
        lowest_pvalue <- min(causal_effects[features, 'p-value'])
        all_stat_sig_features <- as.matrix(which(causal_effects[features, 'p-value']==lowest_pvalue))
        row.names(all_stat_sig_features) <- row.names(causal_effects)[features[all_stat_sig_features]]
        selected_features <- c(selected_features, row.names(all_stat_sig_features))
        train_data_subset <- as.matrix(train_data[, row.names(all_stat_sig_features)])
        test_data_subset <- as.matrix(test_data[, row.names(all_stat_sig_features)])
        if(i==1){
          combined.newFeature.train <- train_data_subset
          combined.newFeature.test  <- test_data_subset
        }else{
          combined.newFeature.train <- cbind(combined.newFeature.train, train_data_subset)
          combined.newFeature.test  <- cbind(combined.newFeature.test, test_data_subset)
        }
    }
    colnames(combined.newFeature.train) <- colnames(combined.newFeature.test) <- selected_features
    return(list(combined.newFeature.train, combined.newFeature.test))
}

find.maximum.distance.point <- function(
  y, 
  x=1:length(y)){   
  allCoord = rbind(y, x)
  
  firstPoint = allCoord[,1]
  lineVec = allCoord[,length(y)] - firstPoint
  lineVecN = lineVec / sqrt(sum(lineVec^2))
  
  vecFromFirst = allCoord - firstPoint
  scalarProduct = lineVecN %*% vecFromFirst
  
  vecFromFirstParallel = t(scalarProduct) %*% lineVecN
  vecToLine = t(vecFromFirst) - vecFromFirstParallel
  distToLine = sqrt(rowSums(vecToLine^2,2))
  which.max(distToLine)
}

kfold_cv <- function(data, k){
  num_years <- nrow(data)
  num_classes <- 3
  subsample_size <- round(num_years/k)
  org_data <- data
  sample_size <- nrow(data)
  season <- 10:12
  col_names <- paste("q_RIBFLV", season, sep="_")

  # Comment this for loop to run LOOCV
  # for(j in 1:10){
    
    # dir_name <- 'G:\\NC\\summer\\Mandar\\fwddeliverables\\climateData\\sahel'
    # dir_name <- paste(dir_name, '/CV_Run_', j, sep='')
    # dir.create(dir_name)

    randomized_idx <- sample(1:num_years, num_years, replace=FALSE)
    randomized_data <- data[randomized_idx, ]
    subsample_list <- list()

    for(i in 1:k){
      if(i==k){
        subsample_list[[i]] <- subsample_data <- randomized_data
      }else{
        subsample_list[[i]] <- subsample_data <- randomized_data[1:subsample_size, ]
        randomized_data <- randomized_data[-(1:subsample_size), ]
      }
    }

    for(i in 1:length(subsample_list)){
      # Leave the ith fold out and build train and test data
      train_data <- do.call(rbind.data.frame, subsample_list[-i])
      test_data <- subsample_list[[i]]

      all_data <- rbind(train_data, test_data)
      all_data <- all_data[order(row.names(all_data)), ]

      all_years <- matrix(1:nrow(all_data), nrow=nrow(all_data), ncol=1, dimnames=list(c(row.names(all_data))))
      test_years <- all_years[row.names(test_data), 1]

      # Data pre-processing for bio data as mentioned in IDA paper, standardize the data so that all the covariates have unit variance
      train_data_sd <- apply(train_data, 2, sd)
      process_train_data <- matrix(0, nrow=nrow(train_data), ncol=ncol(train_data), dimnames=list(row.names(train_data), colnames(train_data)))
      process_test_data <- matrix(0, nrow=nrow(test_data), ncol=ncol(test_data), dimnames=list(row.names(test_data), colnames(test_data)))

      # Preprocess climate indices: detrend and normalize.
      process_data <- preprocessData(test_years, all_data, all_years, col_names)
      process_train_data <- process_data[[1]]
      process_test_data <- process_data[[2]]
      process_rainfall_data <- process_data[[3]]
      colnames(process_rainfall_data) <- 'q_RIBFLV'
      norm_train_data <- cbind(process_train_data, process_rainfall_data)

      data_dir_name <-''
      write.table(norm_train_data, paste('ClimData_CV_', i, '.txt', sep=""), row.names=FALSE, sep='\t')
      write.table(norm_train_data, paste('TrainData_CV_', i, '.txt', sep=""), row.names=TRUE, sep='\t')
      write.table(process_test_data, paste('TestData_CV_', i, '.txt', sep=""), row.names=TRUE, sep='\t')

    }
}

computeCausalEffects <- function(idx, sampleClimData, pcGraph, response_idx){
  causalEffectStats <- ida.pcr.pvalue(idx, response_idx, sampleClimData, pcGraph, verbose=FALSE, y.notparent=TRUE, method="local")
  causalEffects <- causalEffectStats[[1]]
  causalEffectPValue <- causalEffectStats[[2]]
  newStats <- cbind(causalEffects, causalEffectPValue)

  no_na_causalEffects <- which(!is.na(newStats[, 1]))
  if(length(no_na_causalEffects)>0){
    causalEffects <- matrix(0, nrow=length(no_na_causalEffects), ncol=ncol(newStats))
    causalEffects[, 1:2] <- newStats[no_na_causalEffects, ]
  }
  else{
    causalEffects <- NULL
  }

  if(length(causalEffects)==0)
    return(NaN)
  else{
    idx <- sort.int(abs(causalEffects[, 1]), index.return=TRUE)$ix[1]
    return(causalEffects)
  }
}


check.new.coll <- function(amat,amatSkel,x,pa1,pa2.t,pa2.f) {
  ## Check if undirected edges that are pointed to x create a new v-structure
  ## Additionally check, if edges that are pointed away from x create
  ## new v-structure; i.e. x -> pa <- papa would be problematic
  ## pa1 are definit parents of x
  ## pa2 are undirected "parents" of x
  ## pa2.t are the nodes in pa2 that are directed towards pa2
  ## pa2.f are the nodes in pa2 that are directed away from pa2
  ## Value is TRUE, if new collider is introduced
  res <- FALSE
  if ((length(pa2.t)>0) & any(!is.na(pa2.t))) {
    ## check whether all pa1 and all pa2.t are connected;
    ## if no, there is a new collider
    if ((length(pa1)>0) & any(!is.na(pa1))) {
      res <- (min(amatSkel[pa1,pa2.t])==0) ## TRUE if new collider
    }
    ## in addition, all pa2.t have to be connected
    if ((length(pa2.t)>1) & (!res)) {
      tmp <- amatSkel[pa2.t,pa2.t]
      diag(tmp) <- 1
      res2 <- (min(tmp)==0) ## TRUE if new collider
      res <- (res|res2)
    }
  }
  if (!res & ((length(pa2.f)>0) & any(!is.na(pa2.f)))) {
    ## consider here only the DIRECTED Parents of pa2.f
    ## remove undirected edges
    amatTmp <- amat
    amatTmp <- amatTmp-t(amatTmp)
    amatTmp[amatTmp<0] <- 0
    tmp <- amatTmp[pa2.f,,drop=FALSE]
    ## find parents of pa2.f
    papa <- setdiff(which(apply(tmp,2,sum)!=0),x)
    ## if any node in papa is not directly connected to x, there is a new
    ## collider
    if (length(papa)==0) {
      res3 <- FALSE
    } else {
      res3 <- (min(amatSkel[x,papa])==0) ## TRUE if new collider
    }
    res <- (res|res3)
  }
  res
}
