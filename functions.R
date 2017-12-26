# Causal Inference Package
library(pcalg)
library(graph)

# Local Causal Discovery Algorithm
library(bnlearn)

# Data processing
library(pracma)
library(binhf)
library(XML)
library(tools)
library(permute)

# Multicollinearity
library(pls)

# Classifier libraries
library(rpart)
library(C50)
library(randomForest)
library(e1071)
library(class)
library(nnet)
library(gbm)
library(kernlab)

localCausalFS <- function(topKFeatures, alg, response_idx){
  
    col_names <- colnames(topKFeatures)
  blacklist <- matrix(0, nrow=(ncol(topKFeatures)-1), ncol=2)
    blacklist[, 1] <- col_names[response_idx]
    blacklist[, 2] <- col_names[-response_idx]

    # cat('Local Feature Selection: ', alg, '\n')
    if(alg=='MMMB'){
      obj <- mmpc(topKFeatures, blacklist=blacklist, alpha=0.05, optimized=TRUE, test='zf')
      selectedFeatures <- mmpc.features <- getFeatures(obj, response_idx)
    }else if(alg=='HITON-MB'){
      obj <- si.hiton.pc(topKFeatures, blacklist=blacklist, alpha=0.05, optimized=TRUE, test='zf')
      selectedFeatures <- hiton.features <- getFeatures(obj, response_idx)
    }else if(alg=='ARACNE'){
        obj <- aracne(topKFeatures, blacklist=blacklist, mi='mi-g')
        selectedFeatures <- aracne.features <- getFeatures(obj, response_idx)
    }else if(alg=='MMHC'){
        obj <- mmhc(topKFeatures, blacklist=blacklist, alpha=0.05)
        selectedFeatures <- mmhc.features <- getFeatures(obj, response_idx)
    }
    return(selectedFeatures)
}

getFeatures <- function(obj, response_idx){
  features <- union(obj$nodes[[response_idx]]$mb, obj$nodes[[response_idx]]$nbr)
  # features <- union(c(), obj$nodes[[response_idx]]$nbr)
  return(features)
}

calcAccuracy <- function (vec1, vec2){
  vec <- vec1 == vec2
  acc <- length(which(vec == TRUE)) / length(vec1)
  return(acc * 100)
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
  colnames(rainfallData) <- "Rainfall"
  
  for(rain_idx in 1:ncol(rainfallSeasonData)){  # Aggregate the detrended values of rainfall season into one column
    rainfallData <- rainfallData + rainfallSeasonData[, rain_idx]
  }

  rainfallNormData <- normalizeData(rainfallData, 'TRAIN')    # Normalized the aggregated rainfall values

  return(list(sampleClimData, testClimData, rainfallNormData))
}

discVariable <- function(rainfallData){
  
  discVector <- matrix(0, nrow=nrow(rainfallData), ncol=1)

  # Find the HIGH anomalous years for Var1_Month[i]
  sortedVector <- sort(rainfallData, decreasing=TRUE,index.return=TRUE)
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

normalizeData <- function(climData, flag, stats=0){
  if(flag=='TRAIN' || flag=='ALL'){
    for(i in 1:ncol(climData)){
      climData[, i] <- (climData[, i] - mean(climData[, i]))/sd(climData[, i])
    }
  }else{
    for(i in 1:ncol(climData)){
      climData[, i] <- (climData[, i] - stats[[1]][i])/(stats[[2]][i])  
    }
  }
  return(climData)
}

getDiscretizedRainfall <- function(orgClimData, rainfall_cols){
  orgRainfallSeasonData <- orgClimData[, rainfall_cols]
  detrendOrgRainfallData <- matrix(0, nrow=nrow(orgRainfallSeasonData), ncol=ncol(orgRainfallSeasonData))
  orgRainfallData <- matrix(0, nrow=nrow(orgRainfallSeasonData), ncol=1)
  # Detrend and aggregate
  detrendOrgRainfallData <- apply(orgRainfallSeasonData, 2, detrend)

  for(i in 1:ncol(detrendOrgRainfallData)){
    orgRainfallData <- orgRainfallData + detrendOrgRainfallData[, i]
  }
  
  # Normalize and discretize rainfall
  orgNormRainfallData <- normalizeData(orgRainfallData, 'TRAIN')
  discOrgRainfallData <- discVariable(orgNormRainfallData)
  return(discOrgRainfallData)
}

buildModel <- function(topKFeatures, oneYearDataKFeatures, features=NULL, response, featureSelection=FALSE){
  # print(featureSelection)
  if(featureSelection){
    print("Doing feature selection")
    trainData <- topKFeatures[, c(features, 'Rainfall')]
    testData <- oneYearDataKFeatures[, features, drop=FALSE]
    colnames(testData) <- features
    print(dim(trainData))
  }else{
    # print('No feature selection')
    trainData <- topKFeatures
    testData <- oneYearDataKFeatures
  }
  # set.seed(123)

  response_idx <- which(colnames(trainData)=='Rainfall')

  predictor <- as.matrix(trainData[, -response_idx])
  colnames(predictor) <- colnames(trainData)[-response_idx]
  trainData$Rainfall <- as.factor(response)
    
  # build_c50 <- C5.0(Rainfall~., trainData, method="class", trials=1)
  build_c50 <- C5.0(x=predictor, y=trainData$Rainfall, method="class", trials=1)
  predict_c50 <- predict(build_c50, data.frame(testData, check.names=FALSE))

  # build_ruleset <- C5.0(Rainfall~., trainData, method="class", trials=1, rules=TRUE)
  # build_ruleset <- C5.0(x=predictor, y=trainData$Rainfall, method="class", trials=1, rules=TRUE)
  # predict_ruleset <- predict(build_ruleset, data.frame(testData, check.names=FALSE))

  # build_nb <- naiveBayes(Rainfall~., trainData)
  # build_nb <- naiveBayes(x=predictor, y=trainData$Rainfall)
  # predict_nb <- predict(build_nb, data.frame(testData, check.names=FALSE))

  # build_svm <- svm(x=predictor, y=trainData$Rainfall)
  # predict_svm <- predict(build_svm, data.frame(testData, check.names=FALSE))

  # build_ksvm_1 <- ksvm(x=predictor, y=trainData$Rainfall, kernel='vanilladot')
  # predict_ksvm_1 <- predict(build_ksvm_1, data.frame(testData, check.names=FALSE))

  # build_ksvm_2 <- ksvm(x=predictor, y=trainData$Rainfall, kernel='rbfdot')
  # predict_ksvm_2 <- predict(build_ksvm_2, data.frame(testData, check.names=FALSE))    

  # build_knn <- knn1(as.matrix(trainData[, -response_idx]), testData, response)

  # build_boosted_tree <- C5.0(x=predictor, y=trainData$Rainfall, method="class", trials=10)
  # predict_boosted_tree <- predict(build_boosted_tree, data.frame(testData, check.names=FALSE))

  # gbm_fit <- gbm.fit(x=predictor, y=trainData$Rainfall, distribution="multinomial", n.trees=50)
  # gbm_predictions <- predict(gbm_fit, data.frame(testData, check.names=FALSE), n.trees=gbm_fit$n.trees, type="response")
  # gbm_predictions <- apply(gbm_predictions, 1, which.max)

  return(list(build_c50, predict_c50, build_svm=NULL, predict_svm=NULL, build_ksvm_1=NULL, predict_ksvm_1=NULL, build_ksvm_2=NULL, predict_ksvm_2=NULL, build_knn=NULL))
}

# Do not use this function
testSignificance <- function(idx_name, beta, data, pred_idx, response_idx, flag){
  alpha <- 0.05
  result1 <- computePValue(idx_name, beta, data, pred_idx, response_idx, flag)
  result2 <- performZTest(data, beta, alpha, response_idx)
  if(!is.na(result1) && !is.na(result2)){
    return (beta)
  }
}

performZTest <- function(data, beta, alpha, response_idx){
  response_vector <- data
  mu_data <- mean(response_vector)
  sd_data <- sd(response_vector)
  prob <- pnorm(q = beta, mean = mu_data, sd = sd_data, lower.tail = FALSE)
  if(prob < alpha){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

computePValue <- function(idx_name, beta, data, pred_idx, response_idx, flag){

    N <- nrow(data)
    response_vector <- data[, response_idx]
    # Initialize count to 1 for the original estimated causal effect
    count <- 1
    set.seed(1234201)
    alpha <- 0.05
    rand_perm <- 100
    beta_hat_vec <- vector(length = rand_perm)
    # print("Beta")
    # print(beta)
    randomized_row_set <- as.matrix(shuffleSet(n=nrow(data), nset=rand_perm))
    # Do multi-core using foreach instead of for
    for(i in 1:rand_perm){
      rand_idx <- randomized_row_set[i, ]
      data[, response_idx] <- data[rand_idx, response_idx]
      if(flag==1)
        reg_model <- pcr(Rainfall~., data=data[, c(pred_idx, response_idx)])
      else
        reg_model <- plsr(Rainfall~., data=data[, c(pred_idx, response_idx)], method="oscorespls")
      max_var_comp <-  which.max(explvar(reg_model))
      beta_hat <- data.frame(reg_model$coefficients)[idx_name, max_var_comp]
      beta_hat_vec[i] = beta_hat
      # print(beta_hat)
      # Hypothesis: New causal effect beta_hat, is atleast the estimated causal effect beta
      if(abs(beta_hat) >= abs(beta)){
        count <- count + 1 
      }
      #if((count/(rand_perm+1)) > alpha){
      #    return(NA) 
      #}
            
      data[, response_idx] <- response_vector
    }
    p_value <- count/(rand_perm+1)
    # print("Printing count")
    print(cat('Index', idx_name, ' p_value: ', p_value))
    result1 <- FALSE
    result2 <- FALSE
    if(p_value <= alpha){
      #print('significant causality coefficient -> pvalue test')
      result1 <- TRUE
    }
#     result2 <- performZTest(beta_hat_vec, beta, alpha, response_idx)
    result2 <- TRUE
    print(paste('result1 ',result1,'result2 ',result2, sep = '-'))
    if(result1 && result2){
      return(beta)
    }else{
      return(NA)
    }
}

cdfd <- function(train_data, test_data, pc_graph, dir_edge_rank, response_pos, option){
  n <- nrow(train_data)
  p <- ncol(train_data)
  col_names <- colnames(train_data)
  direct_idx <- array(0, dim=n)
  bool_vect <- array(FALSE, dim=p)
  causal_effects <- matrix(0, nrow=p, ncol=1)

  count <- 0
  remove_idx <- na_idx <- c()
  names_edges <- names_direct_edges <- c()  
  remove_idx <- na_idx <- direct_edge_train <- direct_edge_idx <- direct_edge_test <- c()

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
      # print('Direct edge with rainfall', colnames(train_data)[from])
      if(!bool_vect[from]){
        bool_vect[from] <- TRUE
        from_effects <- getCausalEffects.pcr.2(from, train_data, pc_graph, response_pos)  # Use all years to get the causal effect
        # from_effects <- getCausalEffects.plsr.2(from, train_data, pc_graph, response_pos)
        if(length(from_effects)==1 && is.na(from_effects)){
          na_idx <- union(na_idx, from)
          remove_idx <- union(remove_idx, j)
          dir_edge_rank[j, 3] <- -1
          next
        }
        causal_effects[from, ] <- from_effects
      }
      # Construct a feature CE(A) x A for A -> Rainfall and CE(A) = Causal Effect of A on rainfall
      direct_edge_idx <- c(direct_edge_idx, j)
      direct_idx <- causal_effects[from, ]*train_data[, from]
      # direct_idx <- train_data[, from]
      direct_edge_train <- cbind(direct_edge_train, direct_idx)

      direct_edge_test <- cbind(direct_edge_test, as.matrix(causal_effects[from, ]*test_data[, from]))
      # direct_edge_test <- cbind(direct_edge_test, test_data[, from])      
      
      names_direct_edges <- c(names_direct_edges, paste(colnames(train_data)[from], 'NEW', sep='_'))
      # names_direct_edges <- c(names_direct_edges, colnames(train_data)[from])
      next
    }

    # For each directed edge A->B, get the causal effect of indices in a directed causal edge on rainfall
      if(!bool_vect[from]){
        # print(cat('From: ', colnames(train_data)[from]))
        bool_vect[from] <- TRUE
        from_effects <- getCausalEffects.pcr.2(from, train_data, pc_graph, response_pos)
        # from_effects <- getCausalEffects.plsr.2(from, train_data, pc_graph, response_pos)

        if(length(from_effects)==1 && is.na(from_effects)){
          causal_effects[from, ] <- 0
          na_idx <- union(na_idx, from)
          remove_idx <- union(remove_idx, j)
        }else
          causal_effects[from, ] <- from_effects
      }
      
      if(!bool_vect[to]){
        # print(cat('To: ',colnames(train_data)[to]))
        bool_vect[to] <- TRUE
        to_effects <- getCausalEffects.pcr.2(to, train_data, pc_graph, response_pos)
        # to_effects <- getCausalEffects.plsr.2(to, train_data, pc_graph, response_pos)
        
        if(length(to_effects)==1 && is.na(to_effects)){
          causal_effects[to, ] <- 0
          na_idx <- union(na_idx, to)
          remove_idx <- union(remove_idx, j)
        }else
          causal_effects[to, ] <- to_effects
      }
      # Sum the absolute values of causal effects to rank the edges
      if(sum(abs(causal_effects[from,])) > 0 && sum(abs(causal_effects[to,])) > 0){
          dir_edge_rank[j, 3] <- sum(abs(causal_effects[from, ]), abs(causal_effects[to, ]))
      }else{
          dir_edge_rank[j, 3] <- -1 
      }
  }

  if(length(which(dir_edge_rank[,1] == 0)))
    print('Not all directed edges are in dir_edge_rank')

  delete_edges <- union(direct_edge_idx, remove_idx)
  if(length(delete_edges) > 0){
    temp <- nrow(dir_edge_rank)-length(delete_edges)
    new_dir_edge_rank <- matrix(0, nrow=temp, ncol=3)
    new_dir_edge_rank[1:temp, ] <- dir_edge_rank[-delete_edges, ]
  }else{
    new_dir_edge_rank <- dir_edge_rank
  }
          
  new_dir_edge_rank[1:temp,] <- new_dir_edge_rank[order(new_dir_edge_rank[, 3], decreasing=TRUE), ]
    
  if(as.numeric(option)==1 || as.numeric(option)==3){
    K <- nrow(new_dir_edge_rank)  # Take all edges as features 
    all_edges <- matrix(0, nrow=K, ncol=ncol(new_dir_edge_rank))
    all_edges <- new_dir_edge_rank
  }else{
    if(nrow(new_dir_edge_rank)>1)
      temp_indices <- which(new_dir_edge_rank[,3] > median(new_dir_edge_rank[,3]))  # Take edges with causal effects greater than the median of all causal effects
    else 
      temp_indices<-1
    K <- length(temp_indices)
    all_edges <- matrix(0, nrow=K, ncol=ncol(new_dir_edge_rank))
    all_edges[1:K, ] <- new_dir_edge_rank[temp_indices, ]
  }

  if(nrow(all_edges) > 0){
    selected_features <- union(all_edges[, 1], all_edges[, 2])
    names_edges <- paste(pc_graph@nodes[all_edges[, 1]], '->', pc_graph@nodes[all_edges[, 2]], sep="")
    # names_edges <- col_names[selected_features]
    train_edge_features <- constructFeatures.2(all_edges, train_data, causal_effects, response_pos)  # Construct features CE(A) x z-score(A) + CE(B) x z-score(B) 
    test_edge_features <- constructFeatures.test.2(all_edges, test_data, causal_effects) # Create test data for one year
    train_edge_features <- cbind(train_edge_features, direct_edge_train)   # Aggregate data (edges + direct relations) and column names for training and testing
    test_edge_features <- cbind(test_edge_features, direct_edge_test)
  }else if(nrow(all_edges) == 0 && is.null(direct_edge_train)){
      print(sprintf("No new features in network %d",i))
      return(NULL)
  }else{
      train_edge_features <- direct_edge_train
      test_edge_features <- direct_edge_test
  } 
    
  if(as.numeric(option)==3){
    # train_data[, response_pos] <- tr
    # Aggregate the constructed features in training data with existing ones and their column names
    train_edge_features <- cbind(train_data, train_edge_features)    
    colnames(train_edge_features) <- c(col_names, names_edges, names_direct_edges)

    # Combine the new constructed features (edges + direct relations + climate indices) and normalize them
    test_edge_features <- cbind(test_data, test_edge_features)
    test_edge_features <- data.frame(test_edge_features, check.names=FALSE)

    colnames(test_edge_features) <- c(col_names[-response_pos], names_edges, names_direct_edges)
  }else{
    train_edge_features <- cbind(train_edge_features, train_data[, response_pos])
    colnames(train_edge_features) <- c(names_edges, names_direct_edges, col_names[response_pos])

    colnames(test_edge_features) <- c(names_edges, names_direct_edges)
     
    train_edge_features <- data.frame(train_edge_features, check.names=FALSE)
    test_edge_features <- data.frame(test_edge_features, check.names=FALSE)
  }
  return(list(train_edge_features, test_edge_features)) 
}

cdfd_clust <- function(train_data, test_data, pc_graph, dir_edge_rank, response_pos, option){
  n <- nrow(train_data)
  p <- ncol(train_data)
  col_names <- colnames(train_data)
  direct_idx <- array(0, dim=n)
  bool_vect <- array(FALSE, dim=p)
  causal_effects <- matrix(0, nrow=p, ncol=1)

  count <- 0
  remove_idx <- na_idx <- c()
  names_edges <- names_direct_edges <- c()  
  remove_idx <- na_idx <- direct_edge_train <- direct_edge_idx <- direct_edge_test <- c()

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
      # print('Direct edge with rainfall', colnames(train_data)[from])
      if(!bool_vect[from]){#check if the causal effect is already calculated or not.
        bool_vect[from] <- TRUE
        from_effects <- getCausalEffects.pcr.2(from, train_data, pc_graph, response_pos)  # Use all years to get the causal effect
        # from_effects <- getCausalEffects.plsr.2(from, train_data, pc_graph, response_pos)
        if(length(from_effects)==1 && is.na(from_effects)){
          na_idx <- union(na_idx, from)
          remove_idx <- union(remove_idx, j)
          dir_edge_rank[j, 3] <- -1
          next
        }
        causal_effects[from, ] <- from_effects
      }
      # Construct a feature CE(A) x A for A -> Rainfall and CE(A) = Causal Effect of A on rainfall
      direct_edge_idx <- c(direct_edge_idx, j)
      direct_idx <- causal_effects[from, ]*train_data[, from]
      # direct_idx <- train_data[, from]
      direct_edge_train <- cbind(direct_edge_train, direct_idx)

      direct_edge_test <- cbind(direct_edge_test, as.matrix(causal_effects[from, ]*test_data[, from]))
      # direct_edge_test <- cbind(direct_edge_test, test_data[, from])      
      
      names_direct_edges <- c(names_direct_edges, paste(colnames(train_data)[from], 'NEW', sep='_'))
      # names_direct_edges <- c(names_direct_edges, colnames(train_data)[from])
      next
    }

    # For each directed edge A->B, get the causal effect of indices in a directed causal edge on rainfall
      if(!bool_vect[from]){
        # print(cat('From: ', colnames(train_data)[from]))
        bool_vect[from] <- TRUE
        from_effects <- getCausalEffects.pcr.2(from, train_data, pc_graph, response_pos)
        # from_effects <- getCausalEffects.plsr.2(from, train_data, pc_graph, response_pos)

        if(length(from_effects)==1 && is.na(from_effects)){
          causal_effects[from, ] <- 0
          na_idx <- union(na_idx, from)
          remove_idx <- union(remove_idx, j)
        }else
          causal_effects[from, ] <- from_effects
      }
      
      if(!bool_vect[to]){
        # print(cat('To: ',colnames(train_data)[to]))
        bool_vect[to] <- TRUE
        to_effects <- getCausalEffects.pcr.2(to, train_data, pc_graph, response_pos)
        # to_effects <- getCausalEffects.plsr.2(to, train_data, pc_graph, response_pos)
        
        if(length(to_effects)==1 && is.na(to_effects)){
          causal_effects[to, ] <- 0
          na_idx <- union(na_idx, to)
          remove_idx <- union(remove_idx, j)
        }else
          causal_effects[to, ] <- to_effects
      }
      # Sum the absolute values of causal effects to rank the edges
      if(sum(abs(causal_effects[from,])) > 0 && sum(abs(causal_effects[to,])) > 0){
          dir_edge_rank[j, 3] <- sum(abs(causal_effects[from, ]), abs(causal_effects[to, ]))
      }else{
          dir_edge_rank[j, 3] <- -1 
      }
  }

  if(length(which(dir_edge_rank[,1] == 0)))
    print('Not all directed edges are in dir_edge_rank')

  delete_edges <- union(direct_edge_idx, remove_idx)
  if(length(delete_edges) > 0){
    temp <- nrow(dir_edge_rank)-length(delete_edges)
    new_dir_edge_rank <- matrix(0, nrow=temp, ncol=3)
    new_dir_edge_rank[1:temp, ] <- dir_edge_rank[-delete_edges, ]
  }else{
    new_dir_edge_rank <- dir_edge_rank
  }
          
  new_dir_edge_rank[1:temp,] <- new_dir_edge_rank[order(new_dir_edge_rank[, 3], decreasing=TRUE), ]
    
  if(as.numeric(option)==1 || as.numeric(option)==3){
    K <- nrow(new_dir_edge_rank)  # Take all edges as features 
    all_edges <- matrix(0, nrow=K, ncol=ncol(new_dir_edge_rank))
    all_edges <- new_dir_edge_rank
  }else{
    if(nrow(new_dir_edge_rank)>1)
      temp_indices <- which(new_dir_edge_rank[,3] > median(new_dir_edge_rank[,3]))  # Take edges with causal effects greater than the median of all causal effects
    else 
      temp_indices<-1
    K <- length(temp_indices)
    all_edges <- matrix(0, nrow=K, ncol=ncol(new_dir_edge_rank))
    all_edges[1:K, ] <- new_dir_edge_rank[temp_indices, ]
  }

  if(nrow(all_edges) > 0){
#     selected_features <- union(all_edges[, 1], all_edges[, 2])
#     names_edges <- paste(pc_graph@nodes[all_edges[, 1]], '->', pc_graph@nodes[all_edges[, 2]], sep="")
#     # names_edges <- col_names[selected_features]
#     train_edge_features <- constructFeatures.2(all_edges, train_data, causal_effects, response_pos)  # Construct features CE(A) x z-score(A) + CE(B) x z-score(B) 
#     test_edge_features <- constructFeatures.test.2(all_edges, test_data, causal_effects) # Create test data for one year
#     train_edge_features <- cbind(train_edge_features, direct_edge_train)   # Aggregate data (edges + direct relations) and column names for training and testing
#     test_edge_features <- cbind(test_edge_features, direct_edge_test)
    
      new_data <- constructFeatures.3(all_edges, train_data, test_data, causal_effects)  # Construct features CE(A) x z-score(A) + CE(B) x z-score(B) 
      train_edge_features <- new_data[[1]]
      test_edge_features <- new_data[[2]]
#       print('Got new features')
#       print(cat('New features %d',dim(train_edge_features)))
#       print(head(train_edge_features))
#       print(colnames(train_edge_features))
      new_features <- colnames(train_edge_features)
#       test_edge_features <- constructFeatures.test.3(all_edges, test_data, causal_effects) # Create test data for one year
#       train_edge_features <- cbind(train_edge_features, direct_edge_train)   # Aggregate data (edges + direct relations) and column names for training and testing
#       test_edge_features <- cbind(test_edge_features, direct_edge_test)
  }else if(nrow(all_edges) == 0 && is.null(direct_edge_train)){
      print(sprintf("No new features in network %d",i))
      return(NULL)
  }else{
      train_edge_features <- direct_edge_train
      test_edge_features <- direct_edge_test
  } 
    
  if(as.numeric(option)==3){
    # train_data[, response_pos] <- tr
    # Aggregate the constructed features in training data with existing ones and their column names
    train_edge_features <- cbind(train_data, train_edge_features)    
    colnames(train_edge_features) <- c(col_names, names_edges, names_direct_edges)

    # Combine the new constructed features (edges + direct relations + climate indices) and normalize them
    test_edge_features <- cbind(test_data, test_edge_features)
    test_edge_features <- data.frame(test_edge_features, check.names=FALSE)

    colnames(test_edge_features) <- c(col_names[-response_pos], names_edges, names_direct_edges)
  }else{
    # print('error')
    train_edge_features <- cbind(train_edge_features, train_data[, response_pos])
#     print(dim(train_edge_features))
#     print(length(new_features))
#     colnames(train_edge_features) <- c(names_edges, names_direct_edges, col_names[response_pos])
# colnames(test_edge_features) <- c(names_edges, names_direct_edges)

    colnames(train_edge_features) <- c(new_features, col_names[response_pos])
    colnames(test_edge_features) <-  new_features
     
    train_edge_features <- data.frame(train_edge_features, check.names=FALSE)
    test_edge_features <- data.frame(test_edge_features, check.names=FALSE)
  }
  return(list(train_edge_features, test_edge_features)) 
}

constructFeatures.3 <- function(all_edges, train_data, test_data, causal_effects){
  row.names(causal_effects) <- colnames(train_data)
  non_zero_features <- which(causal_effects[,1] != 0)
  causal_effects <- causal_effects[non_zero_features,]
    
#   ratios <- array(0, dim=(length(causal_effects)-2))
#   clusters <- list()
#   clusters[[1]] <- 'temp'
#   for(i in 2:(length(causal_effects)-2)){
# #     print(cat('row num ', length(causal_effects), ' i ', i))
#     temp <- causal_effects
#     clusters[[i]] <- kmeans(x = temp, centers = i)
#     ratios[i] <- 1 - clusters[[i]]$tot.withinss / clusters[[i]]$totss
#   }
#   opt_k <- find.maximum.distance.point(ratios[-1])+1
#   print('opt_k########################')
#   print(opt_k)
  # plot(1:length(ratios), x = ratios, type = 'l')
  # clusters <- kmeans(x = causal_effects, centers = opt_k)
  res <- NbClust(causal_effects, distance = "euclidean", min.nc=2, max.nc=(2/3*nrow(causal_effects)), method = "kmeans", index = "ch")
  opt_k <- res$Best.nc

  feature_names <- paste('Feature', 1:opt_k, sep='_')
  combined.newFeature.train <- matrix(0, nrow=nrow(train_data), ncol=opt_k, dimnames=list(c(), feature_names))
  combined.newFeature.test <- matrix(0, nrow=nrow(test_data), ncol=opt_k, dimnames=list(c(), feature_names))
  
  myclust <- data.frame(clusters[[opt_k]]$cluster)
  for(i in 1:opt_k){
    #print(which(myclust == i))
    
    features <- which(res$Best.partition == i)
    train_data_subset <- as.matrix(train_data[, row.names(myclust)[features]])
    test_data_subset <- as.matrix(test_data[, row.names(myclust)[features]])
    #print(data_subset[1,])
#     combined.newFeature.train[, i] <- linearCombination(train_data_subset, features, causal_effects)
#     combined.newFeature.test[, i] <- linearCombination(test_data_subset, features, causal_effects)

    if(length(features) > 1){
#       combined.newFeature.train[, i] <- linearCombination(train_data_subset, features, causal_effects)
#       combined.newFeature.test[, i] <- linearCombination(test_data_subset, features, causal_effects)
      combined.newFeature.train[,i] = rowMeans(train_data_subset)
      combined.newFeature.test[,i] = rowMeans(test_data_subset)
    }else{
      combined.newFeature.train[, i] <- train_data_subset
      combined.newFeature.test[, i] <- test_data_subset
    }
    #print(combined.newFeature[1,i])
  }
#   plot(x = causal_effects, col = clusters$cluster)
#   points(clusters$centers,col = 1:4, pch = 8)
#   print(clusters$centers)
# print('returning new features')
  return(list(combined.newFeature.train, combined.newFeature.test))
}

linearCombination <- function(train_data_subset, features, causal_effects){
  temp <- 0
  for(i in 1:length(features)){
    temp <- temp + train_data_subset[, i]*causal_effects[features[i]]
  }
  return(temp)
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
  season <- 7:9
  col_names <- paste("Rainfall", season, sep="_")

  for(j in 1:10){
    
    dir_name <- 'G:\\NC\\summer\\Mandar\\fwddeliverables\\climateData\\sahel'
    dir_name <- paste(dir_name, '/CV_Run_', j, sep='')
    dir.create(dir_name)

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
      colnames(process_rainfall_data) <- 'Rainfall'
      norm_train_data <- cbind(process_train_data, process_rainfall_data)

      data_dir_name <- paste(dir_name, '/CV_', i, sep='')
      dir.create(data_dir_name)
      write.table(norm_train_data, paste(data_dir_name, '/ClimData_CV_', i, '.txt', sep=""), row.names=FALSE, sep='\t')
      write.table(norm_train_data, paste(data_dir_name, '/TrainData_CV_', i, '.txt', sep=""), row.names=TRUE, sep='\t')
      write.table(process_test_data, paste(data_dir_name, '/TestData_CV_', i, '.txt', sep=""), row.names=TRUE, sep='\t')
    }
  }
}

getCausalEffects.pcr.2 <- function(idx, sampleClimData, pcGraph, response_idx){
  #print(response_idx)
  causalEffect <- ida.pcr.pvalue(idx, response_idx, sampleClimData, pcGraph, y.notparent=TRUE, method="local")
  causalEffect <- causalEffect[which(!is.na(causalEffect))]

  if(length(causalEffect)==0)
    return(NaN)
  else{
    # print(causalEffect)
    idx <- sort.int(abs(causalEffect), index.return=TRUE)$ix[1]
    return(causalEffect[idx])
  }
}

constructFeatures.2 <- function(topKEdges, climData, causalEffects, response_idx){
  
  combined.newFeature <- matrix(0, nrow=nrow(climData), ncol=nrow(topKEdges))
  from.newFeature <- array(0, dim=nrow(climData))
  to.newFeature <- array(0, dim=nrow(climData))
  
  for(i in 1:nrow(topKEdges)){
    from <- topKEdges[i, 1]
    to <- topKEdges[i, 2]
    
    from.newFeature <- causalEffects[from, ]*climData[, from]
    to.newFeature <- causalEffects[to, ]*climData[, to]

    combined.newFeature[, i] <- from.newFeature + to.newFeature
  }
  return(combined.newFeature)
}

constructFeatures.test.2 <- function(topKEdges, testData, causalEffects){
  combined.newFeature <- matrix(0, nrow=nrow(testData), ncol=nrow(topKEdges))
  
  for(j in 1:nrow(topKEdges)){
    from <- topKEdges[j, 1]
    to <- topKEdges[j, 2]
    from.newFeature <- causalEffects[from, ]*testData[, from]
    to.newFeature <- causalEffects[to, ]*testData[, to]
    combined.newFeature[, j] <- as.matrix(from.newFeature + to.newFeature)
  }
  return(combined.newFeature)  
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
