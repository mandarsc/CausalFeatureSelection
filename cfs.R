rm(list=ls())
#home <- '/home/mschaudh/CausalInference_Experiments/ClimateData'
home <- 'G:\\Mandar\\NCSU\\Independent Study\\Climate Project\\Dhara\\R Code'
setwd(home)
source("functions.R")
source("ida.pcr.pvalue_2.R")
home <- 'G:\\Mandar\\NCSU\\Independent Study\\Climate Project\\Dhara\\8Fold\\Sahel'
setwd(home)

library(FSelector)
library(class)
library(C50)
library(e1071)
library(enpls)
library(pracma)
# library(caret)
library(glmnet)
col_names <- paste("Rainfall", 7:9, sep="_")
  
file_path <- file.path(paste(home, sep="/"))
print(cat('Reading files from...',file_path))
org_clim_data <- read.table(file=file.path(paste(file_path, 'Climate Raw Data Sahel.csv', sep="/")), header=T, row.names=1, check.names=0,
  sep=",")

num_years <- nrow(org_clim_data)
years <- matrix(1:57, nrow=num_years, ncol=1, dimnames=list(1951:2007))

HIGH <- 1
NORMAL <- 2
LOW <- 3
option <- 1
# Store classification models and results from discClimData and sample_clim_data
actual_class <- array(0, 3)
build_c50   <- build_boosted_tree      <- build_gbm  <- build_nb   <- build_knn  <- build_ksvm_1 <- build_ksvm_2 <- list()
predict_c50 <- predict_gbm<- predict_boosted_tree    <- predict_nb <-  predict_ksvm_1 <- predict_ksvm_2 <- list()
result_c50  <- result_gbm <- result_boosted_tree     <- result_nb  <- result_knn <- result_ksvm_1 <- result_ksvm_2 <- list()
c50_class_count <- ruleset_class_count  <- svm_class_count  <- nb_class_count <- knn_class_count  <- array(0, 3)

disc_org_rainfall_data <-matrix(0, nrow=num_years, ncol=1)
disc_org_rainfall_data[, 1] <- getDiscretizedRainfall(org_clim_data, col_names)

cdfd_sc_accuracy <- matrix(0, nrow=10, ncol=40)
cdfd_pc_accuracy <- matrix(0, nrow=10, ncol=40)
cdfd_cs_accuracy <- matrix(0, nrow=10, ncol=40)
cdfd_ig_accuracy <- matrix(0, nrow=10, ncol=40)
cdfd_gr_accuracy <- matrix(0, nrow=10, ncol=40)
cdfd_su_accuracy <- matrix(0, nrow=10, ncol=40)
cdfd_rf_accuracy <- matrix(0, nrow=10, ncol=40)
cdfd_lasso_accuracy <- matrix(0, nrow=10, ncol=8)
cdfd_accuracy <- matrix(0, nrow=10, ncol=8)

for(cv_run in 1:10){
    # cv_run <- as.numeric(readline(prompt="Enter CV run no: "))
    cv_run_dir <- paste('CV_Run', cv_run, sep='_')

    #setwd(paste(home, region, 'Cross_Validation', cv_run_dir, sep='/'))
    setwd(paste(home, cv_run_dir, sep='/'))

    cv_dirs <- grep('CV_[1-8]', list.dirs(path='.', recursive=FALSE), value=TRUE)
    num_folds <- length(cv_dirs)

    dt_accuracy <- boosted_tree_accuracy <-  gbm_accuracy <- ksvm_1_accuracy <- ksvm_2_accuracy <- array(0, num_folds)

    dt_cdfd_sc_accuracy <- matrix(0, nrow=num_folds, ncol=40)
    dt_cdfd_pc_accuracy <- matrix(0, nrow=num_folds, ncol=40)
    dt_cdfd_cs_accuracy <- matrix(0, nrow=num_folds, ncol=40)
    dt_cdfd_ig_accuracy <- matrix(0, nrow=num_folds, ncol=40)
    dt_cdfd_gr_accuracy <- matrix(0, nrow=num_folds, ncol=40)
    dt_cdfd_su_accuracy <- matrix(0, nrow=num_folds, ncol=40)
    dt_cdfd_rf_accuracy <- matrix(0, nrow=num_folds, ncol=40)
    dt_cdfd_lasso_accuracy <- matrix(0, nrow=num_folds, ncol=1)
    # dt_cdfd_hiton_accuracy <- matrix(0, nrow=num_folds, ncol=1)
    # dt_cdfd_mmmb_accuracy <- matrix(0, nrow=num_folds, ncol=1)
    # dt_cdfd_aracne_accuracy <- matrix(0, nrow=num_folds, ncol=1)
    # dt_cdfd_sr_accuracy <- matrix(0, nrow=num_folds, ncol=1)
    # dt_cdfd_log_reg_accuracy <- matrix(0, nrow=num_folds, ncol=1)

    dt_sc_accuracy <- matrix(0, nrow=num_folds, ncol=40)
    avg_num_features <- 0

    for(i in 1:num_folds){
    #   i=1
        start_time <- Sys.time()
        print(start_time)

        train_data_file <- list.files(path=cv_dirs[i], pattern='^TrainData*.*.txt')
        test_data_file <- list.files(path=cv_dirs[i], pattern='TestData*.*.txt')

        print(train_data_file)
        
        train_data <- read.table(paste(cv_dirs[i], train_data_file, sep='/'), header=TRUE, sep='\t')
        test_data <- read.table(paste(cv_dirs[i], test_data_file, sep='/'), header=TRUE, row.names=1, sep='\t')

        cg_file_path <- file.path(paste(cv_dirs[i], 'Causal_Graphs_v10', sep='/'))

        # Read Causal Graph from training data
        pc_graph_file <- list.files(path=cg_file_path, pattern='^ClimData_*.*.xml')
        cg_filename <- paste(file_path_sans_ext(pc_graph_file), 'xml', sep='.')
        print(cg_filename)

        print(cat('Reading files from...', file.path(cg_file_path)))

        test_years <- array(0, dim=nrow(test_data))
        test_years <- years[row.names(test_data), 1]
        true_class <- array(0, dim=length(test_years))

        cv_no <- na.omit(as.numeric(unlist(strsplit(pc_graph_file, "[^0-9]+"))))[1]
        pc_graph <- tetradToR(paste(cg_file_path, cg_filename, sep="/"))    # Get the pc_graph object from the xml file

        pc_wm <- wgtMatrix(pc_graph)
        pc_wmD <- t(pc_wm - t(pc_wm))

        dir_pc <- which(pc_wmD==1, arr.ind=TRUE)  # Get list of directed edges from the graph
        dir_edge_rank <- matrix(0, nrow=nrow(dir_pc), ncol=3) # Store the edges and their causal effect on rainfall
        dir_edge_rank[, 1:2] <- dir_pc
        causal_effects <- matrix(0, nrow=length(pc_graph@nodes), ncol=1)   # Store the causal effect of each index on rainfall in anomalous phases

        rainfall_idx <- which(pc_graph@nodes=='Rainfall')

        disc_rainfall_data <- discVariable(as.matrix(train_data[, 'Rainfall']))   # Discretize the normalized rainfall 
        colnames(disc_rainfall_data) <- 'Rainfall'
        true_class <- disc_org_rainfall_data[test_years, 1]

        constructed_data <- cdfd_clust(train_data, test_data, pc_graph, dir_edge_rank, rainfall_idx, option=1)
        if(is.null(constructed_data)){
            print(cat('Graph ', pc_graph_file[i], ' does not find any directed edges as features'))
            next
        }
        train_edge_features <- constructed_data[[1]]
        test_edge_features <- constructed_data[[2]]
        write.csv(train_edge_features, paste(cv_dirs[i], 'CDFD_ST_TrainData.csv', sep='/'))
        write.csv(train_edge_features, paste(cv_dirs[i], 'CDFD_ST_TestData.csv', sep='/'))
        # Lasso
        lasso_model <- glmnet(as.matrix(train_edge_features[, -ncol(train_edge_features)]), train_data$Rainfall, alpha=1) #alpha=1 is the lasso penalty
        coefficient_set <- (coef(lasso_model, s=0.01) != 0)@i               #coeffs for min lambda for Lasso --> 0.01
        new_training_data <- train_edge_features[, coefficient_set]
        new_testing_data <-  test_edge_features[,coefficient_set]
        new_training_data$Rainfall <- disc_rainfall_data
        model_lasso <- buildModel(new_training_data, new_testing_data, NULL, response=disc_rainfall_data)

        #Local Feature selection
        # response_idx <- which(colnames(train_edge_features)=='Rainfall')
        # features_hiton  <- localCausalFS(train_edge_features, 'HITON-MB', response_idx)
        # features_mmmb   <- localCausalFS(train_edge_features, 'MMMB', response_idx)
        # features_aracne <- localCausalFS(train_edge_features, 'ARACNE', response_idx)
        # dt_cdfd_hiton_accuracy[cv_no, 1] <- 0
        # dt_cdfd_mmmb_accuracy[cv_no, 1] <- 0
        # dt_cdfd_aracne_accuracy[cv_no, 1] <- 0

        # if(length(features_hiton)>0){
        #     new_training_data <- as.matrix(train_edge_features[, c(features_hiton, 'Rainfall')])
        #     new_testing_data <-  as.matrix(test_edge_features[, features_hiton])
        #     colnames(new_testing_data) <- features_hiton
        #     new_training_data[, 'Rainfall'] <- disc_rainfall_data
        #     model_hiton <- buildModel(data.frame(new_training_data), data.frame(new_testing_data), NULL, disc_rainfall_data)
        #     build_tree_hiton    <- model_hiton[[1]];    predict_cdfd_hiton <- model_hiton[[2]]
        #     result_cdfd_hiton <- predict_cdfd_hiton == true_class
        #     dt_cdfd_hiton_accuracy[cv_no, 1] <- length(which(result_cdfd_hiton==TRUE))/length(test_years)*100
        # }
        # if(length(features_mmmb)>0){
        #     new_training_data <- as.matrix(train_edge_features[, c(features_mmmb, "Rainfall")])
        #     new_testing_data <-  as.matrix(test_edge_features[,features_mmmb])
        #     colnames(new_testing_data) <- features_mmmb
        #     new_training_data[, 'Rainfall'] <- disc_rainfall_data
        #     model_mmmb <- buildModel(data.frame(new_training_data), data.frame(new_testing_data), NULL, disc_rainfall_data)
        #     build_tree_mmmb    <- model_mmmb[[1]];    predict_cdfd_mmmb <- model_mmmb[[2]]
        #     result_cdfd_mmmb <- predict_cdfd_mmmb == true_class
        #     dt_cdfd_mmmb_accuracy[cv_no, 1] <- length(which(result_cdfd_mmmb==TRUE))/length(test_years)*100
        # }
        # if(length(features_aracne)>0){
        #     new_training_data <- as.matrix(train_edge_features[, c(features_aracne, 'Rainfall')])
        #     new_testing_data <-  as.matrix(test_edge_features[, features_aracne])
        #     colnames(new_testing_data) <- features_aracne
        #     new_training_data[, 'Rainfall'] <- disc_rainfall_data
        #     model_aracne <- buildModel(data.frame(new_training_data), data.frame(new_testing_data), NULL, disc_rainfall_data)
        #     build_tree_aracne    <- model_aracne[[1]];    predict_cdfd_aracne <- model_aracne[[2]]
        #     result_cdfd_aracne <- predict_cdfd_aracne == true_class
        #     dt_cdfd_aracne_accuracy[cv_no, 1] <- length(which(result_cdfd_aracne==TRUE))/length(test_years)*100
        # }
    #     
        build_tree_lasso    <- model_lasso[[1]];    predict_cdfd_lasso <- model_lasso[[2]]
      
    #    build_tree_sr    <- model_sr[[1]];    predict_cdfd_sr <- model_sr[[2]]
    #    build_tree_log_reg  <- model_log_reg[[1]];  predict_cdfd_log_reg <- model_log_reg[[2]]

        result_cdfd_lasso <- predict_cdfd_lasso == true_class
        
    #     result_cdfd_sr <- predict_cdfd_sr == true_class
    #    result_cdfd_log_reg <- predict_cdfd_log_reg== true_class

        dt_cdfd_lasso_accuracy[cv_no, 1] <- length(which(result_cdfd_lasso==TRUE))/length(test_years)*100
        
    #     dt_cdfd_sr_accuracy[cv_no, 1] <- length(which(result_cdfd_sr==TRUE))/length(test_years)*100
    #     dt_cdfd_log_reg_accuracy[cv_no, 1] <- length(which(result_cdfd_log_reg==TRUE))/length(test_years)*100

        topK <- 1:(ncol(train_edge_features)-1)
        rf_train_data <- train_edge_features
        rf_test_data <- test_edge_features
        colnames(rf_train_data)[topK] <- paste("Feature", topK, sep='_')
        colnames(rf_test_data)[topK] <- paste("Feature", topK, sep='_')

        CDFD_SPEARMAN = rank.correlation(Rainfall~., train_edge_features)
        CDFD_PEARSON = rank.correlation(Rainfall~.,train_edge_features)
        CDFD_CHI_SQUARE = chi.squared(Rainfall~.,train_edge_features)
        CDFD_INFO_GAIN = information.gain(Rainfall~.,train_edge_features)
        CDFD_GAIN_RATIO = gain.ratio(Rainfall~.,train_edge_features)
        CDFD_SYMM_UNCERTAINTY = symmetrical.uncertainty(Rainfall~.,train_edge_features)
        CDFD_RANDOM_FOREST = random.forest.importance(Rainfall~., rf_train_data)

        # for(idx in 1:length(topK)){
        #     build_model_k_1 <- buildModel(train_edge_features[, c(cutoff.k(CDFD_SPEARMAN, topK[idx]), 'Rainfall')], test_edge_features, response=disc_rainfall_data)
        #     build_model_k_2 <- buildModel(train_edge_features[, c(cutoff.k(CDFD_PEARSON, topK[idx]), 'Rainfall')], test_edge_features, response=disc_rainfall_data)
        #     build_model_k_3 <- buildModel(train_edge_features[, c(cutoff.k(CDFD_CHI_SQUARE, topK[idx]), 'Rainfall')], test_edge_features, response=disc_rainfall_data)
        #     build_model_k_4 <- buildModel(train_edge_features[, c(cutoff.k(CDFD_INFO_GAIN, topK[idx]), 'Rainfall')], test_edge_features, response=disc_rainfall_data)
        #     build_model_k_5 <- buildModel(train_edge_features[, c(cutoff.k(CDFD_GAIN_RATIO, topK[idx]), 'Rainfall')], test_edge_features, response=disc_rainfall_data)
        #     build_model_k_6 <- buildModel(train_edge_features[, c(cutoff.k(CDFD_SYMM_UNCERTAINTY, topK[idx]), 'Rainfall')], test_edge_features, response=disc_rainfall_data)
        #     build_model_k_7 <- buildModel(rf_train_data[, c(cutoff.k(CDFD_RANDOM_FOREST, topK[idx]), 'Rainfall')], rf_test_data, response=disc_rainfall_data)

        #     build_tree_1 <- build_model_k_1[[1]]; predict_tree_1 <- build_model_k_1[[2]]
        #     build_tree_2 <- build_model_k_2[[1]]; predict_tree_2 <- build_model_k_2[[2]]
        #     build_tree_3 <- build_model_k_3[[1]]; predict_tree_3 <- build_model_k_3[[2]]
        #     build_tree_4 <- build_model_k_4[[1]]; predict_tree_4 <- build_model_k_4[[2]]
        #     build_tree_5 <- build_model_k_5[[1]]; predict_tree_5 <- build_model_k_5[[2]]
        #     build_tree_6 <- build_model_k_6[[1]]; predict_tree_6 <- build_model_k_6[[2]]
        #     build_tree_7 <- build_model_k_7[[1]]; predict_tree_7 <- build_model_k_7[[2]]
            
        #     result_cdfd_sc <-  predict_tree_1 == true_class
        #     result_cdfd_pc <-  predict_tree_2 == true_class
        #     result_cdfd_cs <-  predict_tree_3 == true_class
        #     result_cdfd_ig <-  predict_tree_4 == true_class
        #     result_cdfd_gr <-  predict_tree_5 == true_class
        #     result_cdfd_su <-  predict_tree_6 == true_class
        #     result_cdfd_rf <-  predict_tree_7 == true_class

        #     dt_cdfd_sc_accuracy[cv_no, idx] <- (length(which(result_cdfd_sc==TRUE))/length(test_years))*100
        #     dt_cdfd_pc_accuracy[cv_no, idx] <- length(which(result_cdfd_pc==TRUE))/length(test_years)*100
        #     dt_cdfd_cs_accuracy[cv_no, idx] <- length(which(result_cdfd_cs==TRUE))/length(test_years)*100
        #     dt_cdfd_ig_accuracy[cv_no, idx] <- length(which(result_cdfd_ig==TRUE))/length(test_years)*100
        #     dt_cdfd_gr_accuracy[cv_no, idx] <- length(which(result_cdfd_gr==TRUE))/length(test_years)*100
        #     dt_cdfd_su_accuracy[cv_no, idx] <- length(which(result_cdfd_su==TRUE))/length(test_years)*100
        #     dt_cdfd_rf_accuracy[cv_no, idx] <- length(which(result_cdfd_rf==TRUE))/length(test_years)*100 
        # }

        # dt_cdfd_sc_accuracy[cv_no, idx:40] <- dt_cdfd_sc_accuracy[cv_no, idx]
        # dt_cdfd_pc_accuracy[cv_no, idx:40] <- dt_cdfd_pc_accuracy[cv_no, idx]
        # dt_cdfd_cs_accuracy[cv_no, idx:40] <- dt_cdfd_cs_accuracy[cv_no, idx]
        # dt_cdfd_ig_accuracy[cv_no, idx:40] <- dt_cdfd_ig_accuracy[cv_no, idx]
        # dt_cdfd_gr_accuracy[cv_no, idx:40] <- dt_cdfd_gr_accuracy[cv_no, idx]
        # dt_cdfd_su_accuracy[cv_no, idx:40] <- dt_cdfd_su_accuracy[cv_no, idx]
        # dt_cdfd_rf_accuracy[cv_no, idx:40] <- dt_cdfd_rf_accuracy[cv_no, idx] 

        train_edge_features$Rainfall <- as.factor(disc_rainfall_data)
          
        model <- buildModel(train_edge_features, test_edge_features, response=disc_rainfall_data)
        build_c50[[cv_no]]     <- model[[1]];      predict_c50[[cv_no]]     <- model[[2]]
        avg_num_features <- avg_num_features + nrow(as.matrix(C5imp(build_c50[[cv_no]])))

        # model <- buildModel(train_data, test_data, response=disc_rainfall_data)
        # build_c50[[cv_no]]     <- model[[1]];      predict_c50[[cv_no]]     <- model[[2]]
        # avg_num_features <- avg_num_features + nrow(as.matrix(C5imp(build_c50[[cv_no]])))
        # build_gbm[[cv_no]] <- model[[3]];      predict_gbm[[cv_no]] <- model[[4]]
        # build_ksvm_1[[cv_no]]     <- model[[5]];      predict_ksvm_1[[cv_no]]     <- model[[6]]
        # build_ksvm_2[[cv_no]]      <- model[[7]];      predict_ksvm_2[[cv_no]]      <- model[[8]]
        # build_knn[[cv_no]]     <- model[[9]];

        result_c50[[cv_no]]    <- predict_c50[[cv_no]] == true_class
        # result_gbm[[cv_no]]<- predict_gbm[[cv_no]] == true_class
        # result_ksvm_1[[cv_no]]    <- predict_ksvm_1[[cv_no]] == true_class
        # result_ksvm_2[[cv_no]]     <- predict_ksvm_2[[cv_no]] == true_class
        # result_knn[[cv_no]]    <- build_knn[[cv_no]] == true_class

        dt_accuracy[cv_no] <- length(which(result_c50[[cv_no]]==TRUE))/length(test_years)
        # boosted_tree_accuracy[cv_no] <- length(which(result_boosted_tree[[cv_no]]==TRUE))/length(test_years)
        # gbm_accuracy[cv_no] <- length(which(result_gbm[[cv_no]]==TRUE))/length(test_years)
        # ksvm_1_accuracy[cv_no] <- length(which(result_ksvm_1[[cv_no]]==TRUE))/length(test_years)
        # ksvm_2_accuracy[cv_no] <- length(which(result_ksvm_2[[cv_no]]==TRUE))/length(test_years)
    }

    end_time <- Sys.time()
    print(end_time)

    print('CDFD DT Accuracy')
    print(mean(dt_accuracy))

    dir.create('CDFD_Cluster_MT')
    write.csv(cdfd_accuracy[cv_run, ] <- mean(dt_accuracy)*100, 'CDFD_Cluster_MT\\CDFD_Accuracy.csv')
    # write.csv(cdfd_sc_accuracy[cv_run, ] <- apply(dt_cdfd_sc_accuracy, 2, mean), 'CDFD_Cluster_MT\\CDFD_SC_Accuracy.csv', row.names=FALSE)
    # write.csv(cdfd_pc_accuracy[cv_run, ] <- apply(dt_cdfd_pc_accuracy, 2, mean), 'CDFD_Cluster_MT\\CDFD_PC_Accuracy.csv', row.names=FALSE)
    # write.csv(cdfd_cs_accuracy[cv_run, ] <- apply(dt_cdfd_cs_accuracy, 2, mean), 'CDFD_Cluster_MT\\CDFD_CS_Accuracy.csv', row.names=FALSE)
    # write.csv(cdfd_ig_accuracy[cv_run, ] <- apply(dt_cdfd_ig_accuracy, 2, mean), 'CDFD_Cluster_MT\\CDFD_IG_Accuracy.csv', row.names=FALSE)
    # write.csv(cdfd_gr_accuracy[cv_run, ] <- apply(dt_cdfd_gr_accuracy, 2, mean), 'CDFD_Cluster_MT\\CDFD_GR_Accuracy.csv', row.names=FALSE)
    # write.csv(cdfd_su_accuracy[cv_run, ] <- apply(dt_cdfd_su_accuracy, 2, mean), 'CDFD_Cluster_MT\\CDFD_SU_Accuracy.csv', row.names=FALSE)
    # write.csv(cdfd_rf_accuracy[cv_run, ] <- apply(dt_cdfd_rf_accuracy, 2, mean), 'CDFD_Cluster_MT\\CDFD_RF_Accuracy.csv', row.names=FALSE)

    # write.csv(cdfd_lasso_accuracy[cv_run, ] <- dt_cdfd_lasso_accuracy, 'CDFD_Cluster_MT\\CDFD_Lasso_Accuracy.csv', row.names=FALSE)

    # write.csv(dt_cdfd_hiton_accuracy, 'CDFD_HITON_Accuracy.csv', row.names=FALSE)
    # write.csv(dt_cdfd_mmmb_accuracy, 'CDFD_MMMB_Accuracy.csv', row.names=FALSE)
    # write.csv(dt_cdfd_aracne_accuracy, 'CDFD_ARACNE_Accuracy.csv', row.names=FALSE)
    #write.csv(dt_cdfd_log_reg_accuracy, 'CDFD_Log_Reg_Accuracy.csv', row.names=FALSE)

    # print('CDFD and SC Accuracy')
    # print(apply(dt_cdfd_sc_accuracy, 2, mean))
    # print('SC Accuracy')
    # print(apply(dt_sc_accuracy, 2, mean))
    print(cat('Avg num of features', avg_num_features/8, '\n'))
    # cdfd_sc_ea_1 <- dt_cdfd_sc_accuracy
    # print('Boosted DT Accuracy')
    # print(mean(boosted_tree_accuracy))
    # print('GBM Accuracy')
    # print(mean(gbm_accuracy))
    save.image(paste('CDFD_Clust_MT', cv_run, '_.RData', sep=''))
    #setwd(paste(home, region, 'Cross_Validation', sep='/'))
    setwd(home)    
}

write.csv(cdfd_accuracy, "CDFD_MoM_Accuracy.csv", row.names=FALSE)
# write.csv(cdfd_sc_accuracy, "CDFD_SC_MoM_Accuracy.csv", row.names=FALSE)
# write.csv(cdfd_pc_accuracy, "CDFD_PC_MoM_Accuracy.csv", row.names=FALSE)
# write.csv(cdfd_cs_accuracy, "CDFD_CS_MoM_Accuracy.csv", row.names=FALSE)
# write.csv(cdfd_ig_accuracy, "CDFD_IG_MoM_Accuracy.csv", row.names=FALSE)
# write.csv(cdfd_gr_accuracy, "CDFD_GR_MoM_Accuracy.csv", row.names=FALSE)
# write.csv(cdfd_su_accuracy, "CDFD_SU_MoM_Accuracy.csv", row.names=FALSE)
# write.csv(cdfd_rf_accuracy, "CDFD_RF_MoM_Accuracy.csv", row.names=FALSE)
# write.csv(cdfd_lasso_accuracy, "CDFD_LASSO_MoM_Accuracy.csv", row.names=FALSE)