rm(list=ls())
setwd('G:\\Mandar\\NCSU\\Independent Study\\Climate Project\\Code\\CausalFeatureSelection')
source("functions.R")
source("ida.pcr.pvalue.R")
home <- 'G:\\Mandar\\NCSU\\Independent Study\\Climate Project\\Code'

region_option <- 1
# Choose only one option to select the dataset below
# Option 1: Sahel
dataset <- "Sahel"
col_names <- paste("Rainfall", 7:9, sep="_")
file_path <- file.path(paste(home, dataset, sep="\\"))
file_name <- file.path(paste(file_path, 'Climate Raw Data Sahel.csv', sep="\\"))

# Option 2: EA
dataset <- "EA"
col_names <- paste("Rainfall", 10:12, sep="_")
file_path <- file.path(paste(home, dataset, sep="\\"))
file_name <- file.path(paste(file_path, 'Climate Raw Data EA.csv', sep="\\"))

# Option 3: Riboflavin
dataset <- "Riboflavin"
col_names <- 'q_RIBFLV'
file_path <- file.path(paste(home, dataset, sep="\\"))
file_name <- file.path(paste(file_path, 'riboflavinv100.csv', sep="\\"))

dataset <- "ADNI"
col_names <- "MMSE"
gender <- "MALE"
# gender <- "FEMALE"
file <- paste('ADNI_', gender, '_MMSE.csv', sep='')
file_path <- file.path(paste(home, dataset, sep="\\"))
file_name <- file.path(paste(file_path, file, sep="\\"))


print(cat('Reading files from...',file_path))
if(dataset=="Sahel" || dataset=="EA" || dataset=="Riboflavin"){
    raw_data <- read.table(file=file_name, header=T, row.names=1, check.names=0, sep=",")
}else{
    raw_data <- read.table(file=file_name, header=T, row.names=1, check.names=1, sep=",")
}

if(dataset=="Sahel" || dataset=="EA"){
    n_samples <- nrow(raw_data)
    row_ids <- matrix(1:n_samples, nrow=n_samples, ncol=1, dimnames=list(1951:2007))
    y_var <- 'Rainfall'
    graph_prefix <- 'Clim'    
}

if(dataset == "Riboflavin"){
    raw_data <- t(raw_data)
    n_samples <- nrow(raw_data)
    row_ids <- matrix(1:n_samples, nrow=n_samples, ncol=1, dimnames=list(row.names(raw_data)))    
    y_var <- 'q_RIBFLV'
    graph_prefix <- 'Bio'
}

if(dataset == "ADNI"){
    n_samples <- nrow(raw_data)
    row_ids <- matrix(1:n_samples, nrow=n_samples, ncol=1, dimnames=list(row.names(raw_data)))    
    y_var <- 'MMSE'
    graph_prefix <- paste('ADNI_', gender, sep='')
}

build_c50   <- build_lm <- list()
predict_c50 <- predict_lm   <- list()

# Generate discretized response variable
disc_response <-matrix(0, nrow=n_samples, ncol=1)
if(dataset == "ADNI"){
    norm_response <- normalizeData(raw_data[, y_var, drop=FALSE], 'TRAIN')
    disc_response[, 1] <- raw_data[, 'DX']
}else if(dataset=="Riboflavin"){
    norm_response <- normalizeData(raw_data[, col_names, drop=FALSE], 'TRAIN')
    disc_response[, 1] <- discVariable(norm_response)
}else{
    y_data <- apply(raw_data[, col_names, drop=FALSE], 2, detrend)
    y_data <- as.matrix(apply(y_data, 1, sum))
  
    norm_response <- normalizeData(y_data, 'TRAIN')
    disc_response[, 1] <- getDiscretizedResponse(raw_data, col_names)
}

# Initialize vairables to record predictions and accuracy
num_folds <- nrow(raw_data)

cfs_accuracy <- mmhc_accuracy <- hiton_accuracy <- pc_accuracy <- sc_accuracy <- ig_accuracy <- lasso_accuracy <- stepr_accuracy <- oner_accuracy <- array(0, num_folds)
cfs_pred <- rf_pred <- mmhc_pred <- hiton_pred <- pc_pred <- sc_pred <- ig_pred <- lasso_pred <- stepr_pred <- oner_pred <- array(NA, num_folds)
y_true <- array(0, num_folds)
avg_num_features <- 0

cfs_features <- list()
hiton_features <- mmhc_features <- sc_features <- pc_features <- ig_features <- oner_features <- stepr_features <- lasso_features <- list()
feature_selection <- TRUE

for(i in 1:num_folds){
    start_time <- Sys.time()
    print(cat('Fold:', i))
    print(start_time)


    if(dataset=="ADNI"){
        X_train_file <- paste(dataset, '_', gender, '_TRAIN_', i, '.txt', sep='')
        X_train <- read.table(paste(file_path, X_train_file, sep='\\'), header=TRUE, sep='\t')        
        X_train_ <- raw_data[-i, colnames(X_train)]
        idx.stats <- list()
        idx.stats[[1]] <- apply(X_train_, 2, mean)
        idx.stats[[2]] <- apply(X_train_, 2, sd)
        X_test_raw <- raw_data[i, colnames(X_train), drop=FALSE]
        X_test <- normalizeData(X_test_raw, 'TEST', idx.stats)

        test_id <- row_ids[row.names(X_test), 1]
        y_true[i] <- norm_response[test_id, 1]
        true_class <- disc_response[i, 1]
        Y_train <- as.matrix(disc_response[-i, 1])
        # Read Causal Graph from training data
        cg_file_path <- file.path(paste(file_path, paste('Causal_Graphs_', 'v10_', gender, sep=''), sep='\\'))
        pc_graph_file <- paste(graph_prefix, '_TRAIN_', i, '_0.05.xml', sep='')

    }else{
        X_train_file <- paste('TrainData_CV_', i, '.txt', sep='')
        X_test_file <- paste('TestData_CV_', i, '.txt', sep='')            
        X_train <- read.table(paste(file_path, X_train_file, sep='\\'), header=TRUE, sep='\t')
        X_test <- read.table(paste(file_path, X_test_file, sep='\\'), header=TRUE, row.names=1, sep='\t')    

        test_id <- row_ids[row.names(X_test), 1]
        y_true[i] <- norm_response[test_id, 1]
        true_class <- disc_response[test_id, 1]
        Y_train <- discVariable(as.matrix(X_train[, y_var]))   # Discretize the normalized rainfall 
        # Read Causal Graph from training data
        cg_file_path <- file.path(paste(file_path, 'Causal_Graphs_v10', sep='\\'))
        pc_graph_file <- paste(graph_prefix, 'Data_CV_', i, '_0.05.xml', sep='')
    }

    # print(X_train_file)
    # colnames(Y_train) <- y_var
    # print(pc_graph_file)
    # print(cat('Reading files from...', file.path(cg_file_path)))

    pc_graph <- tetradToR(paste(cg_file_path, pc_graph_file, sep="\\"))    # Get the pc_graph object from the xml file
    # pc_graph <- mmhcObj(X_train, which(colnames(X_train)==y_var))
    pc_wm <- wgtMatrix(pc_graph)
    pc_wmD <- t(pc_wm - t(pc_wm))

    dir_pc <- which(pc_wmD==1, arr.ind=TRUE)  # Get list of directed edges from the graph
    dir_edge_rank <- matrix(0, nrow=nrow(dir_pc), ncol=3) 
    dir_edge_rank[, 1:2] <- dir_pc
    causal_effects <- matrix(0, nrow=length(pc_graph@nodes), ncol=1)   # Store the causal effect of each index on rainfall in anomalous phases

    y_idx <- which(pc_graph@nodes==y_var)

    constructed_data <- cfs(X_train, X_test, pc_graph, dir_edge_rank, y_idx, option=1)
    if(is.null(constructed_data)){
        print(cat('Graph ', pc_graph_file[i], ' does not find any directed edges as features'))
        next
    }
    X_train_ <- constructed_data[[1]]
    X_test_ <- constructed_data[[2]]
    cfs_features[[i]] <- colnames(X_train_)

    model <- buildModel(X_train_, X_test_, label=Y_train, response=X_train[, y_var, drop=FALSE], y_var)
    build_c50[[i]]  <- model[[1]];      predict_c50[[i]]    <- model[[2]]
    build_lm[[i]]   <- model[[3]];      predict_lm[[i]]     <- model[[4]]
    avg_num_features <- avg_num_features + length(cfs_features[[i]])

    cfs_accuracy[i] <- predict_c50[[i]] == true_class
    cfs_pred[i] <- predict_lm[[i]]
}

print(cat('CFS (Accuracy): ', mean(cfs_accuracy), " CI:", calCI(cfs_accuracy)))
print(cat('CFS (RMSE): ', sqrt(mean((y_true - cfs_pred)^2)), " CI:", calCI(cfs_pred, y_true)))
print(cat('Avg num of features', round(avg_num_features/num_folds), '\n'))

# Compute avg number of features used by CFS and use them for univariate feature selection methods
avg_num_features <- round(avg_num_features/num_folds)

# Perform feature selection
if(feature_selection){
    for(i in 1:num_folds){
        print(cat('Fold: ', i))
        if(dataset=="ADNI"){
            X_train_file <- paste(dataset, '_', gender, '_TRAIN_', i, '.txt', sep='')
            X_train <- read.table(paste(file_path, X_train_file, sep='\\'), header=TRUE, sep='\t')        
            X_train_ <- raw_data[-i, colnames(X_train)]
            idx.stats <- list()
            idx.stats[[1]] <- apply(X_train_, 2, mean)
            idx.stats[[2]] <- apply(X_train_, 2, sd)
            X_test_raw <- raw_data[i, colnames(X_train), drop=FALSE]
            X_test <- normalizeData(X_test_raw, 'TEST', idx.stats)

            test_id <- row_ids[row.names(X_test), 1]
            y_true[i] <- norm_response[test_id, 1]
            true_class <- disc_response[i, 1]
            Y_train <- as.matrix(disc_response[-i, 1])
        }else{
            X_train_file <- paste('TrainData_CV_', i, '.txt', sep='')
            X_test_file <- paste('TestData_CV_', i, '.txt', sep='')            
            X_train <- read.table(paste(file_path, X_train_file, sep='\\'), header=TRUE, sep='\t')
            X_test <- read.table(paste(file_path, X_test_file, sep='\\'), header=TRUE, row.names=1, sep='\t')    
            test_id <- row_ids[row.names(X_test), 1]
            y_true[i] <- norm_response[test_id, 1]
            true_class <- disc_response[test_id, 1]

            Y_train <- discVariable(as.matrix(X_train[, y_var]))   # Discretize the normalized rainfall 
        }
        y_idx <- which(colnames(X_train)==y_var) 
        colnames(Y_train) <- y_var

        # print('Local causal feature selection')
        # hiton_features[[i]] <- localCausalFS(X_train, 'HITON-MB', ncol(X_train))
        mmhc_features[[i]] <- localCausalFS(X_train, 'MMHC', ncol(X_train))
        # if(length(hiton_features[[i]]) > 0){
        #     model_fit1 <- buildModel(X_train[, hiton_features[[i]], drop=FALSE], X_test[, hiton_features[[i]], drop=FALSE], label=Y_train, response=X_train[, y_var, drop=FALSE], y_var)
        #     hiton_accuracy[i] = model_fit1[[2]]==true_class
        #     hiton_pred[i]  <- model_fit1[[4]]                   
        # }
        if(length(mmhc_features[[i]]) > 0){
            model_fit2 <- buildModel(X_train[, mmhc_features[[i]], drop=FALSE], X_test[, mmhc_features[[i]], drop=FALSE], label=Y_train, response=X_train[, y_var, drop=FALSE], y_var)
            mmhc_accuracy[i] <- model_fit2[[2]]==true_class
            mmhc_pred[i] <- model_fit2[[4]]
        }

        # print('Univariate feature selection')
        # univariate_features <- univariateFS(X_train, y_idx, avg_num_features)
        # sc_features[[i]] <- univariate_features[[1]]
        # pc_features[[i]] <- univariate_features[[2]]
        # ig_features[[i]] <- univariate_features[[3]]
        # model_fit3 <- buildModel(X_train[, sc_features[[i]], drop=FALSE], X_test[, sc_features[[i]], drop=FALSE], label=Y_train, response=X_train[, y_var, drop=FALSE], y_var)
        # sc_accuracy[i] <- model_fit3[[2]]==true_class        
        # sc_pred[i] <- model_fit3[[4]]
        # model_fit4 <- buildModel(X_train[, pc_features[[i]], drop=FALSE], X_test[, pc_features[[i]], drop=FALSE], label=Y_train, response=X_train[, y_var, drop=FALSE], y_var)
        # pc_accuracy[i] <- model_fit4[[2]]==true_class
        # pc_pred[i] <- model_fit4[[4]]        
        # model_fit5 <- buildModel(X_train[, ig_features[[i]], drop=FALSE], X_test[, ig_features[[i]], drop=FALSE], label=Y_train, response=X_train[, y_var, drop=FALSE], y_var)
        # ig_accuracy[i] <- model_fit5[[2]]==true_class        
        # ig_pred[i] <- model_fit5[[4]]

        # print('Regression-based feature selection')
        # regression_features <- regressionFS(X_train, y_idx)
        # lasso_features[[i]] <- regression_features[[1]][-1]
        # stepr_features[[i]] <- regression_features[[2]]
        # if(length(lasso_features[[i]]) > 0){
        #     model_fit6 <- buildModel(X_train[, lasso_features[[i]], drop=FALSE], X_test[, lasso_features[[i]], drop=FALSE], label=Y_train, response=X_train[, y_var, drop=FALSE], y_var)
        #     lasso_accuracy[i] <- model_fit6[[2]]==true_class                            
        #     lasso_pred[i] <- model_fit6[[4]]
        # }
        # if(length(stepr_features[[i]]) > 0){
        #     model_fit7 <- buildModel(X_train[, stepr_features[[i]], drop=FALSE], X_test[, stepr_features[[i]], drop=FALSE], label=Y_train, response=X_train[, y_var, drop=FALSE], y_var)
        #     stepr_accuracy[i] <- model_fit7[[2]]==true_class                    
        #     stepr_pred[i] <- model_fit7[[4]]
        # }

        # print('oneR feature selection')
        # oner_features[[i]] <- oneRFS(X_train, y_idx, avg_num_features)[[1]]
        # if(length(oner_features[[i]]) > 0){
        #     model_fit8 <- buildModel(X_train[, oner_features[[i]], drop=FALSE], X_test[, oner_features[[i]], drop=FALSE], label=Y_train, response=X_train[, y_var, drop=FALSE], y_var)
        #     oner_accuracy[i] <- model_fit8[[2]]==true_class
        #     oner_pred[i] <- model_fit8[[4]]                                    
        # }
    }
    print(cat('Pearson:', mean(pc_accuracy), " (RMSE):", sqrt(mean((y_true - pc_pred)^2))))
    print(cat('CI (Accuracy):', calCI(pc_accuracy), " CI (RMSE):", calCI(pc_pred, y_true)))
    print(cat('Spearman:', mean(sc_accuracy), " (RMSE):", sqrt(mean((y_true - sc_pred)^2))))
    print(cat('CI (Accuracy):', calCI(sc_accuracy), " CI (RMSE):", calCI(sc_pred, y_true)))    
    print(cat('OneR:', mean(oner_accuracy), " (RMSE):", sqrt(mean((y_true - oner_pred)^2))))
    print(cat('CI (Accuracy):', calCI(oner_accuracy), " CI (RMSE):", calCI(oner_pred, y_true)))
    print(cat('StepR:', mean(stepr_accuracy), " (RMSE):", sqrt(mean((y_true - stepr_pred)^2))))
    print(cat('CI (Accuracy):', calCI(stepr_accuracy), " CI (RMSE):", calCI(stepr_pred, y_true)))
    print(cat('Lasso:', mean(lasso_accuracy), " (RMSE):", sqrt(mean((y_true - lasso_pred)^2))))
    print(cat('CI (Accuracy):', calCI(lasso_accuracy), " CI (RMSE):", calCI(lasso_pred, y_true)))
    print(cat('HITON (Accuracy):', mean(hiton_accuracy), " (RMSE):", sqrt(mean((y_true - hiton_pred)^2))))
    print(cat('CI (Accuracy):', calCI(hiton_accuracy), " CI (RMSE):", calCI(hiton_pred, y_true)))
    print(cat('Info Gain:', mean(ig_accuracy), " (RMSE):", sqrt(mean((y_true - ig_pred)^2))))
    print(cat('CI (Accuracy):', calCI(ig_accuracy), " CI (RMSE):", calCI(ig_pred, y_true)))
    print(cat('MMHC:', mean(mmhc_accuracy), " (RMSE):", sqrt(mean((y_true - mmhc_pred)^2))))
    print(cat('CI (Accuracy):', calCI(mmhc_accuracy), " CI (RMSE):", calCI(mmhc_pred, y_true)))
}

end_time <- Sys.time()
print(end_time)
