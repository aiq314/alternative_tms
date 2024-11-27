# Function to calculate sigmoid
sigmoid <- function(x){
  return(1.0/(1.0 + exp(-x)))
}

# Function to calculate variable B
calculate_B1 <- function(alpha1, alpha2, beta1, beta2, feature1, feature2) {
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  return(sigmoid(z1) - 0.5 * sigmoid(z2))
}

calculate_B2 <- function(alpha1, alpha2, beta1, beta2, feature1, feature2) {
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  return(1.0 + sigmoid(z1) - 2.0 * sigmoid(z2))
}

# Function to calculate sigmoid
sigmoid <- function(x){
  return(1.0/(1.0 + exp(-x)))
}

# Function to calculate variable B
calculate_B1 <- function(alpha1, alpha2, beta1, beta2, feature1, feature2) {
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  return(sigmoid(z1) - 0.5 * sigmoid(z2))
}

calculate_B2 <- function(alpha1, alpha2, beta1, beta2, feature1, feature2) {
  z1 <- alpha1 - beta1 * feature1 - beta2 * feature2
  z2 <- alpha2 - beta1 * feature1 - beta2 * feature2
  return(1.0 + sigmoid(z1) - 2.0 * sigmoid(z2))
}

run_all_testing <- function(results_folder = 'results',
                             training = 'trainingdf',
                             testing = 'testingdf',
                             features_vector = 'features',
                             units_vector = 'units',
                             is_normalized = FALSE,
                             max_attempts = 1000,
                             num_tests = 10000){
  # Determine if the analysis involves a single feature or multiple features
  is_single <- length(features_vector) == 1
  
  # Construct a log file name based on the number of features
  log_filename <- if (!is_single) {
    paste(features_vector, collapse="-")
  } else {
    features_vector
  }
  logfile <- file(paste(results_folder,"/",paste(log_filename,"training","log.txt", sep="_"), sep = ""), open="wt")
  
  # Dynamically select columns for the analysis based on features
  # Include common columns like "label", "drug_name", "Sample_ID"
  common_cols <- c("label", "drug_name", "Sample_ID")
  training_cols <- as.character(c(features_vector, common_cols))
  testing_cols <- as.character(c(features_vector, common_cols))
  
  # Subset training and testing dataframes based on selected columns
  training <- training[, training_cols]
  testing <- testing[, testing_cols]
  
  # Construct dataset for training and testing
  levels <- c("low", "intermediate", "high")
  values <- c(1, 2, 3)
  training$label <- as.integer(factor(training$label, levels = levels, labels = values))
  testing$label <- as.integer(factor(testing$label, levels = levels, labels = values))
  training$label <- as.factor(training$label)
  testing$label <- as.factor(testing$label)
  
  # Remove errors
  training <- na.omit(training)
  testing <- na.omit(testing)
  
  dimension <- length(features_vector)
  if (is_normalized) {
    # Calculate mean and SD for the first m columns of ordinaldf
    means <- sapply(training[1:dimension], mean)
    sds <- sapply(training[1:dimension], sd)
    tempdf <- training
    tempdf[1:dimension] <- mapply(function(x, mean, sd) (x - mean) / sd, training[1:dimension], means, sds)
    training <- tempdf
    tempdf <- testing
    tempdf[1:dimension] <- mapply(function(x, mean, sd) (x - mean) / sd, testing[1:dimension], means, sds)
    testing <- tempdf
  }
  
  # =======================================================
  # Testing performance stochastically depending on "num_tests"
  # =======================================================
  # Fit the ordinal logistic regression model
  # Prepare the formula for logistic regression based on the number of features
  formula_string <- paste("label ~", paste(features_vector[1:dimension], collapse = " + "))
  formula <- as.formula(formula_string)
  
  # Specify the number of attempts
  fit_results <- fitandgetTMS(training_data = training,
                              testing_data = testing,
                              formula = formula,
                              dimension = dimension,
                              max_attempts = max_attempts,
                              logfile = logfile)
  converged <- fit_results$converged
  if (converged) {
    training <- fit_results$training_data
    testing <- fit_results$testing_data
    mod <- fit_results$mod
    alphas <- as.numeric(mod$zeta)  # Alpha values
    betas <- as.numeric(mod$coefficients)  # Beta coefficients for all features
    training <- cbind(training,predict(mod, newdata = training, type = "probs"))
    training["pred"] <- predict(mod, newdata = training, type = "class")
    testing <- cbind(testing,predict(mod, newdata = testing, type = "probs"))
    testing["pred"] <- predict(mod, newdata = testing, type = "class")
    all_or_loocv <- FALSE
    if (dimension == 2) {
      scatterplotfun(data = training,
                     mod = mod,
                     feature_1 = as.character(features_vector[1]),
                     feature_2 = as.character(features_vector[2]),
                     unit_1 = as.character(units_vector[1]),
                     unit_2 = as.character(units_vector[2]),
                     is_converged = converged,
                     is_training = TRUE,
                     is_legend = FALSE,
                     results_folder = results_folder,
                     is_normalized = is_normalized,
                     all_or_loocv = all_or_loocv)
      scatterplotfun(data = testing,
                     mod = mod,
                     feature_1 = as.character(features_vector[1]),
                     feature_2 = as.character(features_vector[2]),
                     unit_1 = as.character(units_vector[1]),
                     unit_2 = as.character(units_vector[2]),
                     is_converged = converged,
                     is_training = FALSE,
                     is_legend = FALSE,
                     results_folder = results_folder,
                     is_normalized = is_normalized,
                     all_or_loocv = all_or_loocv)
    }
    # Constants of the TMS
    # if (!is_single) {
    th1 <- (alphas[2] - alphas[1] ) / 2.0 + log(1.0 - 2.0 * exp(alphas[1] - alphas[2]))
    th2 <- -th1
    
    # Plot TMS for training dataset
    label_colors <- c("low" = "green", "intermediate" = "blue", "high" = "red")
    training$risk <- ifelse(training$label == 1, "low",
                            ifelse(training$label == 2, "intermediate", "high"))
    
    # Plot TMS for testing dataset
    testing$risk <- ifelse(testing$label == 1, "low",
                            ifelse(testing$label == 2, "intermediate", "high"))
    
    # if (!is_single) {
    tms_name <- "TMS"
    filename_training <- paste(results_folder,"/",paste(log_filename, "training_development_tms.jpg", sep="_"), sep = "")
    filename_testing <- paste(results_folder,"/",paste(log_filename, "testing_development_tms.jpg", sep="_"), sep = "")
    
    # Plot TMS for training dataset
    tmsplotfun(data = training,
               th1 = th1,
               th2 = th2,
               label_colors = label_colors,
               title = "Training dataset",
               file_name = filename_training,
               tms_name = "TMS")
    tmsplotfun(data = testing,
               th1 = th1,
               th2 = th2,
               label_colors = label_colors,
               title = "Testing dataset",
               file_name = filename_testing,
               tms_name = "TMS")
  } else {
    if (dimension == 2) {
      scatterplotfun(data = training,
                     mod = NA,
                     feature_1 = as.character(features_vector[1]),
                     feature_2 = as.character(features_vector[2]),
                     unit_1 = as.character(units_vector[1]),
                     unit_2 = as.character(units_vector[2]),
                     is_converged = converged,
                     is_training = TRUE,
                     is_legend = FALSE,
                     results_folder = results_folder,
                     is_normalized = is_normalized,
                     all_or_loocv = all_or_loocv)
      scatterplotfun(data = testing,
                     mod = NA,
                     feature_1 = as.character(features_vector[1]),
                     feature_2 = as.character(features_vector[2]),
                     unit_1 = as.character(units_vector[1]),
                     unit_2 = as.character(units_vector[2]),
                     is_converged = converged,
                     is_training = FALSE,
                     is_legend = FALSE,
                     results_folder = results_folder,
                     is_normalized = is_normalized,
                     all_or_loocv = all_or_loocv)
    }
    temp_return <- returnNA(log_filename,dimension)
    temp_return["Rank_score_training"] <- NA
    temp_return["Accuracy_cipa_1_training"] <- NA
    temp_return["Accuracy_cipa_2_training"] <- NA
    temp_return["AUC_cipa_1_training"] <- NA
    temp_return["AUC_cipa_2_training"] <- NA
    temp_return["Sensitivity_cipa_1_training"] <- NA
    temp_return["Sensitivity_cipa_2_training"] <- NA
    temp_return["Specificity_cipa_1_training"] <- NA
    temp_return["Specificity_cipa_2_training"] <- NA
    temp_return["LR_positive_cipa_1_training"] <- NA
    temp_return["LR_positive_cipa_2_training"] <- NA
    temp_return["LR_negative_cipa_1_training"] <- NA
    temp_return["LR_negative_cipa_2_training"] <- NA
    temp_return["F1score_cipa_1_training"] <- NA
    temp_return["F1score_cipa_2_training"] <- NA
    temp_return["Classification_error_training"] <- NA
    
    return(temp_return)
  }
  
  # Preallocate the pmeasures dataframe for testing drugs
  pmeasures_testing <- vector("list", num_tests)
  
  # Combine all data
  all_data <- rbind(training,testing)
  
  # Calculate the classification error from the whole dataset
  predicted_labels <- predict(mod, newdata = testing, type = "class")
  pred_err_testing <- (abs(as.integer(predicted_labels) - as.integer(testing$label)))
  
  # Calculate the pmeasures for the whole dataset
  training_drugs <- unique(training$drug_name)
  testing_drugs <- unique(testing$drug_name)
  training$is_training <- 1
  testing$is_training <- 0
  # Iterate through the num_tests tests
  set.seed(1)
  for (i in 1:num_tests) {
    # Sample rows for testing data
    testing_data <- do.call(rbind, lapply(unique(testing$drug_name), function(drug) {
      drug_data <- testing[testing$drug_name == drug, ]
      drug_data[sample(nrow(drug_data), 1), ]
    }))
    
    # Sample rows for training data
    training_data <- do.call(rbind, lapply(unique(training$drug_name), function(drug) {
      drug_data <- training[training$drug_name == drug, ]
      drug_data[sample(nrow(drug_data), 1), ]
    }))
    
    temp_data <- rbind(training_data, testing_data)

    # Calculate all pmeasures
    pmeasures_testing[[i]] <- pmeasuresfun_loocv(data = temp_data,
                                       label_values = values,
                                       all_or_loocv = FALSE,
                                       pairwise_test = TRUE)
  } # i-th test
  
  # Combine all pmeasures results at once
  pmeasures_testing <- do.call(rbind, pmeasures_testing)
  
  writepmeasuresfun(pmeasures = pmeasures_testing,
                    pred_error = pred_err_testing,
                    logfile = logfile,
                    all_or_loocv = all_or_loocv,
                    is_loocv = is_loocv,
                    th1 = th1,
                    th2 = th2)
  
  # Compute pmeasures for whole training data
  pmeasures_training <- pmeasuresfun_loocv(data = training,
                                           label_values = values,
                                           all_or_loocv = TRUE,
                                           pairwise_test = FALSE)
  
  # Calculate the classification error from the whole training dataset
  predicted_labels_training <- predict(mod, newdata = training, type = "class")
  pred_err_training <- (abs(as.integer(predicted_labels_training) - as.integer(training$label)))
  
  writepmeasuresfun(pmeasures = pmeasures_training,
                    pred_error = pred_err_training,
                    logfile = logfile,
                    all_or_loocv = TRUE,
                    is_loocv = is_loocv,
                    th1 = th1,
                    th2 = th2)
  
  # Close the connection to the text file
  close(logfile)
  
  # Rank score for the OLR model
  rank_score_testing <- rankscorefun(pmeasures = pmeasures_testing,
                             pred_error = pred_err_testing,
                             is_normalized = TRUE)
  rank_score_training <- rankscorefun_training(pmeasures = pmeasures_training,
                                      pred_error = pred_err_training,
                                      is_normalized = TRUE)
  
  # Store summary of pmeasures to summarydf
  feature_pair_name <- log_filename
  
  summarydf <- data.frame(
    Feature_Pair = feature_pair_name,
    Accuracy_cipa_1 = quantile(pmeasures_testing$Accuracy_cipa_1, 0.025),
    Accuracy_cipa_2 = quantile(pmeasures_testing$Accuracy_cipa_2, 0.025),
    AUC_cipa_1 = quantile(pmeasures_testing$AUC_cipa_1, 0.025),
    AUC_cipa_2 = quantile(pmeasures_testing$AUC_cipa_2, 0.025),
    Sensitivity_cipa_1 = quantile(pmeasures_testing$Sensitivity_cipa_1, 0.025),
    Sensitivity_cipa_2 = quantile(pmeasures_testing$Sensitivity_cipa_2, 0.025),
    Specificity_cipa_1 = quantile(pmeasures_testing$Specificity_cipa_1, 0.025),
    Specificity_cipa_2 = quantile(pmeasures_testing$Specificity_cipa_2, 0.025),
    LR_positive_cipa_1 = quantile(pmeasures_testing$LR_positive_cipa_1, 0.025),
    LR_positive_cipa_2 = quantile(pmeasures_testing$LR_positive_cipa_2, 0.025),
    LR_negative_cipa_1 = quantile(pmeasures_testing$LR_negative_cipa_1, 0.975),
    LR_negative_cipa_2 = quantile(pmeasures_testing$LR_negative_cipa_2, 0.975),
    F1score_cipa_1 = quantile(pmeasures_testing$F1score_cipa_1, 0.025),
    F1score_cipa_2 = quantile(pmeasures_testing$F1score_cipa_2, 0.025),
    Classification_error = mean(pred_err_testing) + 1.96 * sd(pred_err_testing) / sqrt(length(pred_err_testing)),
    Pairwise_classification_accuracy = quantile(pmeasures_testing$Pairwise, 0.025),
    Rank_score = rank_score_testing
  )
  # Add normalized log likelihood value
  summarydf['logLik'] <- logLik(mod)/nrow(training)
  
  # Add Alphas
  for (i in 1:length(alphas)) {
    summarydf[[paste0("Alpha_", i)]] <- alphas[i]
  }
  
  # Add Betas
  for (i in 1:length(betas)) {
    summarydf[[paste0("Beta_", i)]] <- betas[i]
  }
  
  # Add training scores
  summarydf["Rank_score_training"] <- rank_score_training
  summarydf["Accuracy_cipa_1_training"] <- quantile(pmeasures_training$Accuracy_cipa_1, 0.025)
  summarydf["Accuracy_cipa_2_training"] <- quantile(pmeasures_training$Accuracy_cipa_2, 0.025)
  summarydf["AUC_cipa_1_training"] <- quantile(pmeasures_training$AUC_cipa_1, 0.025)
  summarydf["AUC_cipa_2_training"] <- quantile(pmeasures_training$AUC_cipa_2, 0.025)
  summarydf["Sensitivity_cipa_1_training"] <- quantile(pmeasures_training$Sensitivity_cipa_1, 0.025)
  summarydf["Sensitivity_cipa_2_training"] <- quantile(pmeasures_training$Sensitivity_cipa_2, 0.025)
  summarydf["Specificity_cipa_1_training"] <- quantile(pmeasures_training$Specificity_cipa_1, 0.025)
  summarydf["Specificity_cipa_2_training"] <- quantile(pmeasures_training$Specificity_cipa_2, 0.025)
  summarydf["LR_positive_cipa_1_training"] <- quantile(pmeasures_training$LR_positive_cipa_1, 0.025)
  summarydf["LR_positive_cipa_2_training"] <- quantile(pmeasures_training$LR_positive_cipa_2, 0.025)
  summarydf["LR_negative_cipa_1_training"] <- quantile(pmeasures_training$LR_negative_cipa_1, 0.975)
  summarydf["LR_negative_cipa_2_training"] <- quantile(pmeasures_training$LR_negative_cipa_2, 0.975)
  summarydf["F1score_cipa_1_training"] <- quantile(pmeasures_training$F1score_cipa_1, 0.025)
  summarydf["F1score_cipa_2_training"] <- quantile(pmeasures_training$F1score_cipa_2, 0.025)
  summarydf["Classification_error_training"] <- mean(pred_err_training) + 1.96 * sd(pred_err_training) / sqrt(length(pred_err_training))
  
  return(summarydf)
}

run_all_training <- function(results_folder = 'results',
                    training = 'trainingdf',
                    testing = 'testingdf',
                    features_vector = 'features',
                    units_vector = 'units',
                    is_normalized = FALSE,
                    max_attempts = 1000){
  # Determine if the analysis involves a single feature or multiple features
  is_single <- length(features_vector) == 1
  
  # Construct a log file name based on the number of features
  log_filename <- if (!is_single) {
    paste(features_vector, collapse="-")
  } else {
    features_vector
  }
  logfile <- file(paste(results_folder,"/",paste(log_filename,"training","log.txt", sep="_"), sep = ""), open="wt")
  
  # Dynamically select columns for the analysis based on features
  # Include common columns like "label", "drug_name", "Sample_ID"
  common_cols <- c("label", "drug_name", "Sample_ID")
  training_cols <- as.character(c(features_vector, common_cols))
  testing_cols <- as.character(c(features_vector, common_cols))
  
  # Subset training and testing dataframes based on selected columns
  training <- training[, training_cols]
  testing <- testing[, testing_cols]
  
  # Construct dataset for training and testing
  levels <- c("low", "intermediate", "high")
  values <- c(1, 2, 3)
  training$label <- as.integer(factor(training$label, levels = levels, labels = values))
  testing$label <- as.integer(factor(testing$label, levels = levels, labels = values))
  training$label <- as.factor(training$label)
  testing$label <- as.factor(testing$label)
  
  # Remove errors
  training <- na.omit(training)
  testing <- na.omit(testing)
  
  # Combine all data
  all_data <- rbind(training,testing)
  
  # Normalize the data
  dimension <- length(features_vector)
  if (is_normalized) {
    # Calculate mean and SD for the first 'dimension' columns of all_data
    means <- sapply(all_data[1:dimension], mean)
    sds <- sapply(all_data[1:dimension], sd)
    tempdf <- all_data
    tempdf[1:dimension] <- mapply(function(x, mean, sd) (x - mean) / sd, all_data[1:dimension], means, sds)
    all_data <- tempdf
  }
  
  # ================================
  # Training performance on all data
  # ================================
  # Fit the ordinal logistic regression model
  # Prepare the formula for logistic regression based on the number of features
  formula_string <- paste("label ~", paste(features_vector[1:dimension], collapse = " + "))
  formula <- as.formula(formula_string)
  
  # Specify the number of attempts
  # max_attempts <- 10000
  fit_results <- fitandgetTMS(training_data = all_data,
                              testing_data = all_data,
                              formula = formula,
                              dimension = dimension,
                              max_attempts = max_attempts,
                              logfile = logfile)
  converged <- fit_results$converged
  if (converged) {
    all_data <- fit_results$training_data
    mod <- fit_results$mod
    alphas <- as.numeric(mod$zeta)  # Alpha values
    betas <- as.numeric(mod$coefficients)  # Beta coefficients for all features
    all_data <- cbind(all_data,predict(mod, newdata = all_data, type = "probs"))
    all_data["pred"] <- predict(mod, newdata = all_data, type = "class")
    all_or_loocv <- TRUE
    if (dimension == 2) {
      scatterplotfun(data = all_data,
                    mod = mod,
                    feature_1 = as.character(features_vector[1]),
                    feature_2 = as.character(features_vector[2]),
                    unit_1 = as.character(units_vector[1]),
                    unit_2 = as.character(units_vector[2]),
                    is_converged = converged,
                    is_training = TRUE,
                    is_legend = FALSE,
                    results_folder = results_folder,
                    is_normalized = is_normalized,
                    all_or_loocv = all_or_loocv)
    }
    # Constants of the TMS
    # if (!is_single) {
    th1 <- (alphas[2] - alphas[1] ) / 2.0 + log(1.0 - 2.0 * exp(alphas[1] - alphas[2]))
    th2 <- -th1
    
    # Plot TMS for training dataset
    label_colors <- c("low" = "green", "intermediate" = "blue", "high" = "red")
    all_data$risk <- ifelse(all_data$label == 1, "low",
                             ifelse(all_data$label == 2, "intermediate", "high"))
    
    # if (!is_single) {
    tms_name <- "TMS"
    filename_training <- paste(results_folder,"/",paste(log_filename, "dataset_tms.jpg", sep="_"), sep = "")
    
    # Plot TMS for training dataset
    tmsplotfun(data = all_data,
               th1 = th1,
               th2 = th2,
               label_colors = label_colors,
               title = "Training all dataset",
               file_name = filename_training,
               tms_name = "TMS")
  } else {
    if (dimension == 2) {
      scatterplotfun(data = all_data,
                    mod = NA,
                    feature_1 = as.character(features_vector[1]),
                    feature_2 = as.character(features_vector[2]),
                    unit_1 = as.character(units_vector[1]),
                    unit_2 = as.character(units_vector[2]),
                    is_converged = converged,
                    is_training = TRUE,
                    is_legend = FALSE,
                    results_folder = results_folder,
                    is_normalized = is_normalized,
                    all_or_loocv = all_or_loocv)
    }
    return(returnNA(log_filename,dimension))
  }
  
  # Preallocate the pmeasures dataframe
  pmeasures <- data.frame()
  
  # Calculate the classification error from the whole training dataset
  predicted_labels <- predict(mod, newdata = all_data, type = "class")
  pred_err <- (abs(as.integer(predicted_labels) - as.integer(all_data$label)))
  
  # Calculate the pmeasures for the whole dataset
  training_drugs <- unique(training$drug_name)
  testing_drugs <- unique(testing$drug_name)
  all_data$is_training <- ifelse(all_data$drug_name %in% training_drugs, 1, 0)
  for (sample_id in training$Sample_ID[training$drug_name == training_drugs[1]]) {
    temp_data <- subset(all_data, all_data$Sample_ID == sample_id)
    temppmeasure <- pmeasuresfun_loocv(data = temp_data,
                                       label_values = values,
                                       all_or_loocv = TRUE,
                                       pairwise_test = TRUE)
    pmeasures <- rbind(pmeasures,temppmeasure)
  }

  writepmeasuresfun(pmeasures = pmeasures,
                    pred_error = pred_err,
                    logfile = logfile,
                    all_or_loocv = all_or_loocv,
                    is_loocv = is_loocv,
                    th1 = th1,
                    th2 = th2)
  
  # Close the connection to the text file
  close(logfile)
  
  # Rank score for the OLR model
  rank_score <- rankscorefun(pmeasures = pmeasures,
                             pred_error = pred_err,
                             is_normalized = TRUE)
  
  # Store summary of pmeasures to summarydf
  feature_pair_name <- log_filename
  
  summarydf <- data.frame(
    Feature_Pair = feature_pair_name,
    Accuracy_cipa_1 = quantile(pmeasures$Accuracy_cipa_1, 0.025),
    Accuracy_cipa_2 = quantile(pmeasures$Accuracy_cipa_2, 0.025),
    AUC_cipa_1 = quantile(pmeasures$AUC_cipa_1, 0.025),
    AUC_cipa_2 = quantile(pmeasures$AUC_cipa_2, 0.025),
    Sensitivity_cipa_1 = quantile(pmeasures$Sensitivity_cipa_1, 0.025),
    Sensitivity_cipa_2 = quantile(pmeasures$Sensitivity_cipa_2, 0.025),
    Specificity_cipa_1 = quantile(pmeasures$Specificity_cipa_1, 0.025),
    Specificity_cipa_2 = quantile(pmeasures$Specificity_cipa_2, 0.025),
    LR_positive_cipa_1 = quantile(pmeasures$LR_positive_cipa_1, 0.025),
    LR_positive_cipa_2 = quantile(pmeasures$LR_positive_cipa_2, 0.025),
    LR_negative_cipa_1 = quantile(pmeasures$LR_negative_cipa_1, 0.975),
    LR_negative_cipa_2 = quantile(pmeasures$LR_negative_cipa_2, 0.975),
    F1score_cipa_1 = quantile(pmeasures$F1score_cipa_1, 0.025),
    F1score_cipa_2 = quantile(pmeasures$F1score_cipa_2, 0.025),
    Classification_error = mean(pred_err) + 1.96 * sd(pred_err) / sqrt(length(pred_err)),
    Pairwise_classification_accuracy = quantile(pmeasures$Pairwise, 0.025),
    Rank_score = rank_score
  )
  # Add normalized log likelihood value
  summarydf['logLik'] <- logLik(mod)/nrow(all_data)
  
  # Add Alphas
  for (i in 1:length(alphas)) {
    summarydf[[paste0("Alpha_", i)]] <- alphas[i]
  }
  
  # Add Betas
  for (i in 1:length(betas)) {
    summarydf[[paste0("Beta_", i)]] <- betas[i]
  }
  
  
  return(summarydf)
}

run_all_loocv <- function(results_folder = 'results',
                          training = 'trainingdf',
                          testing = 'testingdf',
                          features_vector = 'features',
                          units_vector = 'units',
                          is_normalized = FALSE,
                          max_attempts = 1000) {
  # Determine if the analysis involves a single feature or multiple features
  is_single <- length(features_vector) == 1
  
  # Construct a log file name based on the number of features
  log_filename <- if (!is_single) {
    paste(features_vector, collapse="-")
  } else {
    features_vector
  }
  logfile <- file(paste(results_folder,"/",paste(log_filename,"loocv","log.txt", sep="_"), sep = ""), open="wt")
  
  # Dynamically select columns for the analysis based on features
  # Include common columns like "label", "drug_name", "Sample_ID"
  common_cols <- c("label", "drug_name", "Sample_ID")
  training_cols <- as.character(c(features_vector, common_cols))
  testing_cols <- as.character(c(features_vector, common_cols))
  
  # Subset training and testing dataframes based on selected columns
  training <- training[, training_cols]
  testing <- testing[, testing_cols]
  
  # Construct dataset for training and testing
  levels <- c("low", "intermediate", "high")
  values <- c(1, 2, 3)
  training$label <- as.integer(factor(training$label, levels = levels, labels = values))
  testing$label <- as.integer(factor(testing$label, levels = levels, labels = values))
  training$label <- as.factor(training$label)
  testing$label <- as.factor(testing$label)
  
  # Remove errors
  training <- na.omit(training)
  testing <- na.omit(testing)
  
  # Combine all data
  all_data <- rbind(training,testing)
  
  # Normalize the data
  dimension <- length(features_vector)
  if (is_normalized) {
    # Calculate mean and SD for the first 'dimension' columns of all_data
    means <- sapply(all_data[1:dimension], mean)
    sds <- sapply(all_data[1:dimension], sd)
    tempdf <- all_data
    tempdf[1:dimension] <- mapply(function(x, mean, sd) (x - mean) / sd, all_data[1:dimension], means, sds)
    all_data <- tempdf
  }
  
  # ================================
  # LOOCV performance on all data
  # ================================
  # All drugs
  drugs <- unique(all_data$drug_name)
  
  # Fit the ordinal logistic regression model
  # Prepare the formula for logistic regression based on the number of features
  formula_string <- paste("label ~", paste(features_vector[1:dimension], collapse = " + "))
  formula <- as.formula(formula_string)
  
  # Specify the number of attempts
  # max_attempts <- 10000
  
  # Preallocated loocv data
  loocv_data <- data.frame()
  for (drug in drugs) {
    converged <- FALSE
    fit_results <- fitandgetTMS(training_data = subset(all_data,all_data$drug_name != drug),
                                testing_data = subset(all_data,all_data$drug_name == drug),
                                formula = formula,
                                dimension = dimension,
                                max_attempts = max_attempts,
                                logfile = logfile)
    converged <- fit_results$converged
    if (converged){
      testing_data <- fit_results$testing_data
      mod <- fit_results$mod
      probs <- predict(mod, newdata = testing_data, type = "probs")
      testing_data["1"] <- probs[1]
      testing_data["2"] <- probs[2]
      testing_data["3"] <- probs[3]
      testing_data["pred"] <- predict(mod, newdata = testing_data, type = "class")
      loocv_data <- rbind(loocv_data,testing_data)
    } else {
      write(sprintf("Cross validating with %s failed! Skipping...",drug),logfile)
    }
  }
  
  # Check whether the loocv_data contains three TdP risk samples
  if (length(unique(loocv_data$label))!=3 | length(unique(loocv_data$pred))!=3) {
    return(returnNA(log_filename))
  }
  
  # Preallocate the pmeasures dataframe
  pmeasures <- data.frame()
  
  # Calculate the classification error from the whole loocv data
  pred_err <- (abs(as.integer(loocv_data$pred) - as.integer(loocv_data$label)))
  
  # Calculate the pmeasures for the whole dataset
  training_drugs <- unique(training$drug_name)
  testing_drugs <- unique(testing$drug_name)
  loocv_data$is_training <- ifelse(loocv_data$drug_name %in% training_drugs, 1, 0)
  
  for (sample_id in training$Sample_ID[training$drug_name == training_drugs[1]]) {
    temp_data <- subset(loocv_data, loocv_data$Sample_ID == sample_id)
    temppmeasure <- pmeasuresfun_loocv(data = temp_data,
                                       label_values = values,
                                       all_or_loocv = TRUE,
                                       pairwise_test = TRUE)
    pmeasures <- rbind(pmeasures,temppmeasure)
  }
  all_or_loocv <- TRUE
  writepmeasuresfun(pmeasures = pmeasures,
                    pred_error = pred_err,
                    logfile = logfile,
                    all_or_loocv = all_or_loocv,
                    is_loocv = is_loocv,
                    th1 = NA,
                    th2 = NA)
  
  # Close the connection to the text file
  close(logfile)
  
  # Rank score for the OLR model
  rank_score <- rankscorefun(pmeasures = pmeasures,
                             pred_error = pred_err,
                             is_normalized = TRUE)
  
  # Store summary of pmeasures to summarydf
  feature_pair_name <- log_filename
  
  summarydf <- data.frame(
    Feature_Pair = feature_pair_name,
    Accuracy_cipa_1 = quantile(pmeasures$Accuracy_cipa_1, 0.025),
    Accuracy_cipa_2 = quantile(pmeasures$Accuracy_cipa_2, 0.025),
    AUC_cipa_1 = quantile(pmeasures$AUC_cipa_1, 0.025),
    AUC_cipa_2 = quantile(pmeasures$AUC_cipa_2, 0.025),
    Sensitivity_cipa_1 = quantile(pmeasures$Sensitivity_cipa_1, 0.025),
    Sensitivity_cipa_2 = quantile(pmeasures$Sensitivity_cipa_2, 0.025),
    Specificity_cipa_1 = quantile(pmeasures$Specificity_cipa_1, 0.025),
    Specificity_cipa_2 = quantile(pmeasures$Specificity_cipa_2, 0.025),
    LR_positive_cipa_1 = quantile(pmeasures$LR_positive_cipa_1, 0.025),
    LR_positive_cipa_2 = quantile(pmeasures$LR_positive_cipa_2, 0.025),
    LR_negative_cipa_1 = quantile(pmeasures$LR_negative_cipa_1, 0.975),
    LR_negative_cipa_2 = quantile(pmeasures$LR_negative_cipa_2, 0.975),
    F1score_cipa_1 = quantile(pmeasures$F1score_cipa_1, 0.025),
    F1score_cipa_2 = quantile(pmeasures$F1score_cipa_2, 0.025),
    Classification_error = mean(pred_err) + 1.96 * sd(pred_err) / sqrt(length(pred_err)),
    Pairwise_classification_accuracy = quantile(pmeasures$Pairwise, 0.025),
    Rank_score = rank_score
  )
  
  return(summarydf)
}

TMS <- function(alphas, betas, features_row){
  z1 <- alphas[1] - sum(betas * features_row)
  z2 <- alphas[2] - sum(betas * features_row)

  return((z1 + z2) / 2.0)
}

pairwisefun<- function(fulltable){
  cmb <- combn(seq_len(nrow(fulltable)), 2)
  mergedtable<-cbind(fulltable[cmb[1,],], fulltable[cmb[2,],])
  validpairidx <- (mergedtable[,4]!=mergedtable[,10])&(!mergedtable[,6]|!mergedtable[,12])
  correctidx1 <- ((mergedtable[,4]>mergedtable[,10])&(mergedtable[,5]>mergedtable[,11]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,5]<mergedtable[,11])) #when predicted class are different
  correctidx2 <- (mergedtable[,5]==3)&(mergedtable[,11]==3)&(((mergedtable[,4]>mergedtable[,10])&(mergedtable[,3]<mergedtable[,9]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,3]>mergedtable[,9]))) #when predicted class are both high
  correctidx3 <- (mergedtable[,5]==1)&(mergedtable[,11]==1)&(((mergedtable[,4]>mergedtable[,10])&(mergedtable[,3]<mergedtable[,9]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,3]>mergedtable[,9]))) #when predicted class are both low
  correctidx4 <- (mergedtable[,5]==2)&(mergedtable[,11]==2)&(((mergedtable[,4]>mergedtable[,10])&(mergedtable[,3]<mergedtable[,9]))|((mergedtable[,4]<mergedtable[,10])&(mergedtable[,3]>mergedtable[,9]))) #when predicted class are both intermediate
  correctidx <- correctidx1|correctidx2|correctidx3|correctidx4
  sum(validpairidx&correctidx)/sum(validpairidx)
}

aucrocfun_loocv <- function(data, label_values){
  auc_scores <- c()
  for (class_label in label_values) {
    if (class_label == 1) {
      actual <- as.integer(data$label != class_label)
    } else {
      actual <- as.integer(data$label == class_label)
    }
    predicted_prob <- data[,c("1","2","3")]
    if (class_label == 1) {
      predicted_prob <- predicted_prob[, 2] + predicted_prob[, 3]
    } else {
      predicted_prob <- predicted_prob[, class_label]
    }
    roc_obj <- roc(actual,
                   unlist(predicted_prob),
                   direction = "<",
                   quiet = TRUE)
    auc_score <- auc(roc_obj)
    auc_scores <- c(auc_scores, auc_score)
  }
  return(auc_scores)
}

pmeasuresfun_loocv <- function(data, label_values, tms_name = "TMS", all_or_loocv = TRUE, pairwise_test = TRUE){
  # Calculate performance measures
  if (all_or_loocv) {
    testing_data <- data
  } else {
    testing_data <- subset(data, data$is_training == 0)
  }
  # Create confusion matrix
  confusion_matrix <- table(testing_data$pred, testing_data$label)
  
  # Calculate the accuracy for each class
  tp_cipa_1 = confusion_matrix[2,2] + confusion_matrix[2,3] + confusion_matrix[3,2] + confusion_matrix[3,3]
  tn_cipa_1 = confusion_matrix[1,1]
  fp_cipa_1 = confusion_matrix[2,1] + confusion_matrix[3,1]
  fn_cipa_1 = confusion_matrix[1,2] + confusion_matrix[1,3]
  
  tp_cipa_2 = confusion_matrix[3,3]
  tn_cipa_2 = confusion_matrix[1,1] + confusion_matrix[1,2] + confusion_matrix[2,1] + confusion_matrix[2,2]
  fp_cipa_2 = confusion_matrix[3,1] + confusion_matrix[3,2]
  fn_cipa_2 = confusion_matrix[1,3] + confusion_matrix[2,3]
  
  f1score_cipa_1 <- 2.0 * tp_cipa_1 / (2.0 * tp_cipa_1 + fp_cipa_1 + fn_cipa_1)
  f1score_cipa_2 <- 2.0 * tp_cipa_2 / (2.0 * tp_cipa_2 + fp_cipa_2 + fn_cipa_2)
  
  accuracy_cipa_1 <- (tp_cipa_1 + tn_cipa_1) / sum(confusion_matrix)
  accuracy_cipa_2 <- (tp_cipa_2 + tn_cipa_2) / sum(confusion_matrix)
  
  sensitivity_cipa_1 <- tp_cipa_1 / (tp_cipa_1 + fn_cipa_1)
  sensitivity_cipa_2 <- tp_cipa_2 / (tp_cipa_2 + fn_cipa_2)
  
  specificity_cipa_1 <- tn_cipa_1 / (tn_cipa_1 + fp_cipa_1)
  specificity_cipa_2 <- tn_cipa_2 / (tn_cipa_2 + fp_cipa_2)
  
  # Add random number to LR+ and LR- to prevent zero devision
  u <- 1e-6
  sd <- 1e-12
  
  lr_positive_cipa_1 <- (sensitivity_cipa_1 + rnorm(1, mean = u, sd = sd)) / (1 - specificity_cipa_1 + rnorm(1, mean = u, sd = sd))
  lr_positive_cipa_2 <- (sensitivity_cipa_2 + rnorm(1, mean = u, sd = sd)) / (1 - specificity_cipa_2 + rnorm(1, mean = u, sd = sd))
  
  lr_negative_cipa_1 <- (1 - sensitivity_cipa_1 + rnorm(1, mean = u, sd = sd)) / (specificity_cipa_1 + rnorm(1, mean = u, sd = sd))
  lr_negative_cipa_2 <- (1 - sensitivity_cipa_2 + rnorm(1, mean = u, sd = sd)) / (specificity_cipa_2 + rnorm(1, mean = u, sd = sd))
  
  auc_scores <- aucrocfun_loocv(testing_data, label_values)
  
  # Calculate the pairwise classification error
  if (pairwise_test) {
    data["risk_label"] <- as.integer(data$label)
    data["risk_pred"] <- as.integer(data$pred)
    data <- data[,c("Sample_ID","drug_name",tms_name,"risk_label","risk_pred","is_training")]
    pairwise <- pairwisefun(data)
  } else {
    pairwise <- NA
  }
  
  # Fill the pmeasures row by row
  pmeasures <- data.frame(
    TP_cipa_1 = tp_cipa_1,
    TP_cipa_2 = tp_cipa_2,
    TN_cipa_1 = tn_cipa_1,
    TN_cipa_2 = tn_cipa_2,
    FP_cipa_1 = fp_cipa_1,
    FP_cipa_2 = fp_cipa_2,
    FN_cipa_1 = fn_cipa_1,
    FN_cipa_2 = fn_cipa_2,
    F1score_cipa_1 = f1score_cipa_1,
    F1score_cipa_2 = f1score_cipa_2,
    Accuracy_cipa_1 = accuracy_cipa_1,
    Accuracy_cipa_2 = accuracy_cipa_2,
    AUC_cipa_1 = auc_scores[1],
    AUC_cipa_2 = auc_scores[3],
    Sensitivity_cipa_1 = sensitivity_cipa_1,
    Sensitivity_cipa_2 = sensitivity_cipa_2,
    Specificity_cipa_1 = specificity_cipa_1,
    Specificity_cipa_2 = specificity_cipa_2,
    LR_positive_cipa_1 = lr_positive_cipa_1,
    LR_positive_cipa_2 = lr_positive_cipa_2,
    LR_negative_cipa_1 = lr_negative_cipa_1,
    LR_negative_cipa_2 = lr_negative_cipa_2,
    Pairwise = pairwise
  )
  
  return(pmeasures)
}

tmsplotfun <- function(data, th1, th2, label_colors, title, file_name, tms_name, tms_unit){
  data$drug_name <- factor(data$drug_name, levels = unique(data$drug_name[order(data$label)]))
  tms <- tms_name
  plot <- ggplot(data, aes_string(x = tms_name, y = "drug_name", fill = "risk")) +
    geom_boxplot(color = "black", width = 0.5, size = 0.2, outlier.size = 0.5, outlier.shape = NA) +
    labs(title = title, x = tms, y = "") +
    geom_vline(xintercept = th1, linetype = "dashed", color = "blue", size = 1)  +
    geom_vline(xintercept = th2, linetype = "dashed", color = "red", size = 1)  +
    scale_fill_manual(values = label_colors) + # Set the fill colors
    theme(plot.title = element_text(size = 20), # Title font size
          # Change axis title font sizes
          axis.title.x = element_text(size = 14), # X axis title font size
          axis.title.y = element_text(size = 14), # Y axis title font size
          # Change axis text font sizes
          axis.text.x = element_text(size = 12), # X axis text font size
          axis.text.y = element_text(size = 12), # Y axis text font size
          # Change legend title and text font sizes
          legend.title = element_text(size = 10), # Legend title font size
          legend.text = element_text(size = 8) # Legend text font size
    )
  ggsave(file_name, plot, width = 8, height = 6, dpi = 900)
}

scatterplotfun <- function(data, 
                           mod = NA, 
                           feature_1, 
                           feature_2, 
                           unit_1, 
                           unit_2, 
                           is_converged, 
                           is_training, 
                           is_legend, 
                           results_folder,
                           is_normalized, 
                           all_or_loocv = TRUE){
  # Check the column index
  data <- data.frame(data)
  idx_model_1 <- as.integer(which(colnames(data) == feature_1))
  idx_model_2 <- as.integer(which(colnames(data) == feature_2))
  idx_label <- as.integer(which(colnames(data) == "label"))
  
  if (is_converged) {
    # Variables from Ordinal Logistic Regression model
    alpha1 <- as.numeric(mod$zeta[1])
    alpha2 <- as.numeric(mod$zeta[2])
    beta1 <- as.numeric(mod$coefficients[1])
    beta2 <- as.numeric(mod$coefficients[2])
    
    # Some descriptions of decision boundaries
    m <- - beta1 / beta2
    c1 <- - 1.0 / beta2 * (- alpha1 + log(1.0 - 2.0 * exp(alpha1 - alpha2)))
    c2 <- 1.0 / beta2 * (alpha1 + log(exp(alpha2 - alpha1) - 2.0))
  }
  
  # Create a meshgrid for contour plotting testing dataset
  x <- seq(min(data[,idx_model_1]), max(data[,idx_model_1]), length.out = 100)
  y <- seq(min(data[,idx_model_2]), max(data[,idx_model_2]), length.out = 100)
  if (is_converged) {
    z1 <- outer(x, y, Vectorize(function(x, y) calculate_B1(alpha1, alpha2, beta1, beta2, x, y)))
    z2 <- outer(x, y, Vectorize(function(x, y) calculate_B2(alpha1, alpha2, beta1, beta2, x, y)))
  }
  if (is_training) {
    title <- "Training dataset"
    file_name <- "training"
  } else {
    title <- "Testing dataset"
    file_name <- "testing"
  }
  if (is_normalized) {
    unit_1 <- ""
    unit_2 <- ""
  }
  if (!all_or_loocv) {
    file_name <- paste(file_name,"development",sep = "_")
  }

  jpeg(paste(results_folder,"/",paste(feature_1,feature_2,file_name,"dataset.jpg",sep = "_"), sep = ""),quality = 100, units = "in", width = 5, height = 5, res = 900)
  plot(data[,idx_model_1], 
       data[,idx_model_2], 
       xlab = paste(feature_1, unit_1, sep = " "), 
       ylab = paste(feature_2, unit_2, sep = " "),
       main = title,
       cex.axis = 1.5, 
       cex.lab = 1.5, 
       cex.main = 1.5, 
       cex = 0.5)
  if (is_legend) {
    legend("bottomright", legend = c("Low", "Intermediate", "High"), fill = c("green", "blue", "red"))
  }
  points(data[,idx_model_1][data$label == "1"], data[,idx_model_2][data$label == "1"], col = "green", cex = 0.5)
  points(data[,idx_model_1][data$label == "2"], data[,idx_model_2][data$label == "2"], col = "blue", cex = 0.5)
  points(data[,idx_model_1][data$label == "3"], data[,idx_model_2][data$label == "3"], col = "red", cex = 0.5)
  if (is_converged) {
    abline(a = c1, b = m, col = "blue", lty = 2, lwd = 2)  # Replace 'a' with the intercept if needed
    abline(a = c2, b = m, col = "red", lty = 2, lwd = 2)  # Replace 'a' with the intercept if needed
  }
  dev.off()
}

pairsdfinitfun <- function(features, units, dimension) {
  if (dimension > length(features)) {
    stop("Dimension cannot be greater than the number of features.")
  }

  # Calculate all possible combinations of features and units
  feature_combinations <- combn(features, dimension, simplify = FALSE)
  unit_combinations <- combn(units, dimension, simplify = FALSE)

  # Initialize an empty list to store data frames for each combination
  pairsdf_list <- list()

  # Loop through each combination and create a data frame
  for (i in seq_along(feature_combinations)) {
    feature_combination <- feature_combinations[[i]]
    unit_combination <- unit_combinations[[i]]

    # Create a data frame with the feature names and units for this combination
    feature_names <- paste0("feature_", seq_len(dimension))
    unit_names <- paste0("unit_", seq_len(dimension))
    tempdf <- setNames(data.frame(t(feature_combination)), feature_names)
    unitdf <- setNames(data.frame(t(unit_combination)), unit_names)

    # Combine features and units into one data frame
    combinedf <- cbind(tempdf, unitdf)

    # Add the combined data frame to the list
    pairsdf_list[[i]] <- combinedf
  }

  # Combine all data frames in the list into a single data frame
  pairsdf <- do.call("rbind", pairsdf_list)

  return(pairsdf)
}

solvepair <- function(results_folder,
                      pair_id,
                      filepath_training,
                      filepath_testing,
                      features_vector,
                      units_vector,
                      is_normalized,
                      max_attempts = 1000,
                      is_loocv){
  # Print the current pair_id and features being processed
  print(paste(c(pair_id, features_vector), collapse = "_"))

  # Read in the training and testing datasets
  training <- read_csv(filepath_training, show_col_types = FALSE)
  testing <- read_csv(filepath_testing, show_col_types = FALSE)

  # Perform all
  if (is_loocv) {
    resultdf <- run_all_loocv(results_folder = results_folder,
                              training = training,
                              testing = testing,
                              features_vector = features_vector,
                              units_vector = units_vector,
                              is_normalized = is_normalized,
                              max_attempts = max_attempts)
  } else {
    resultdf <- run_all_training(results_folder = results_folder,
                              training = training,
                              testing = testing,
                              features_vector = features_vector,
                              units_vector = units_vector,
                              is_normalized = is_normalized,
                              max_attempts = max_attempts)
  }
  
  return(resultdf)
}

writepmeasuresfun <- function(pmeasures, 
                            pred_error,
                            logfile,
                            all_or_loocv = TRUE,
                            is_loocv = TRUE,
                            th1 = NA,
                            th2 = NA){
  
  # Print the Thresholds into logfile
  write("=======================",logfile)
  write(sprintf("TMS thresholds"),logfile)
  write("=======================",logfile)
  write(paste(sprintf('Threshold_1: %.4f ', th1)), logfile)
  write(paste(sprintf('Threshold_2: %.4f ', th2)), logfile)

  # Print the pmeasures dataframe into logfile
  if (is_loocv) {
    write("==============================",logfile)
    write("Performance measures for loocv data",logfile)
    write("==============================",logfile)
  } else if (all_or_loocv ) {
    write("==============================",logfile)
    write("Performance measures for whole training data",logfile)
    write("==============================",logfile)
  } else {
    write("==============================",logfile)
    write("Performance measures for testing data based on num_tests",logfile)
    write("==============================",logfile)
  }
  write(sprintf('AUC cipa 1: %.4f (%.4f, %.4f)',
                median(pmeasures$AUC_cipa_1),
                quantile(pmeasures$AUC_cipa_1, 0.025),
                quantile(pmeasures$AUC_cipa_1, 0.975)),
        logfile)
  write(sprintf('AUC cipa 2: %.4f (%.4f, %.4f)',
                median(pmeasures$AUC_cipa_2),
                quantile(pmeasures$AUC_cipa_2, 0.025),
                quantile(pmeasures$AUC_cipa_2, 0.975)),
        logfile)
  write(sprintf('Accuracy cipa 1: %.4f (%.4f, %.4f)',
                median(pmeasures$Accuracy_cipa_1),
                quantile(pmeasures$Accuracy_cipa_1, 0.025),
                quantile(pmeasures$Accuracy_cipa_1, 0.975)),
        logfile)
  write(sprintf('Accuracy cipa 2: %.4f (%.4f, %.4f)',
                median(pmeasures$Accuracy_cipa_2),
                quantile(pmeasures$Accuracy_cipa_2, 0.025),
                quantile(pmeasures$Accuracy_cipa_2, 0.975)),
        logfile)
  write(sprintf('Sensitivity cipa 1: %.4f (%.4f, %.4f)',
                median(pmeasures$Sensitivity_cipa_1),
                quantile(pmeasures$Sensitivity_cipa_1, 0.025),
                quantile(pmeasures$Sensitivity_cipa_1, 0.975)),
        logfile)
  write(sprintf('Sensitivity cipa 2: %.4f (%.4f, %.4f)',
                median(pmeasures$Sensitivity_cipa_2),
                quantile(pmeasures$Sensitivity_cipa_2, 0.025),
                quantile(pmeasures$Sensitivity_cipa_2, 0.975)),
        logfile)
  write(sprintf('Specificity cipa 1: %.4f (%.4f, %.4f)',
                median(pmeasures$Specificity_cipa_1),
                quantile(pmeasures$Specificity_cipa_1, 0.025),
                quantile(pmeasures$Specificity_cipa_1, 0.975)),
        logfile)
  write(sprintf('Specificity cipa 2: %.4f (%.4f, %.4f)',
                median(pmeasures$Specificity_cipa_2),
                quantile(pmeasures$Specificity_cipa_2, 0.025),
                quantile(pmeasures$Specificity_cipa_2, 0.975)),
        logfile)
  write(sprintf('LR+ cipa 1: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_positive_cipa_1),
                quantile(pmeasures$LR_positive_cipa_1, 0.025),
                quantile(pmeasures$LR_positive_cipa_1, 0.975)),
        logfile)
  write(sprintf('LR+ cipa 2: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_positive_cipa_2),
                quantile(pmeasures$LR_positive_cipa_2, 0.025),
                quantile(pmeasures$LR_positive_cipa_2, 0.975)),
        logfile)
  write(sprintf('LR- cipa 1: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_negative_cipa_1),
                quantile(pmeasures$LR_negative_cipa_1, 0.025),
                quantile(pmeasures$LR_negative_cipa_1, 0.975)),
        logfile)
  write(sprintf('LR- cipa 2: %.4f (%.4f, %.4f)',
                median(pmeasures$LR_negative_cipa_2),
                quantile(pmeasures$LR_negative_cipa_2, 0.025),
                quantile(pmeasures$LR_negative_cipa_2, 0.975)),
        logfile)
  write(sprintf('F1score cipa 1: %.4f (%.4f, %.4f)',
                median(pmeasures$F1score_cipa_1),
                quantile(pmeasures$F1score_cipa_1, 0.025),
                quantile(pmeasures$F1score_cipa_1, 0.975)),
        logfile)
  write(sprintf('F1score cipa 2: %.4f (%.4f, %.4f)',
                median(pmeasures$F1score_cipa_2),
                quantile(pmeasures$F1score_cipa_2, 0.025),
                quantile(pmeasures$F1score_cipa_2, 0.975)),
        logfile)
  write(sprintf('Classification error: %.4f (%.4f, %.4f)',
                mean(pred_error),
                mean(pred_error) - 1.96 * sd(pred_error) / sqrt(length(pred_error)),
                mean(pred_error) + 1.96 * sd(pred_error) / sqrt(length(pred_error))),
        logfile)
  if (all_or_loocv) {
    write(sprintf('Pairwise classification accuracy: NA'),
          logfile)
  } else {
    write(sprintf('Pairwise classification accuracy: %.4f (%.4f, %.4f)',
                  median(pmeasures$Pairwise),
                  quantile(pmeasures$Pairwise, 0.025),
                  quantile(pmeasures$Pairwise, 0.975)),
          logfile)
  }
    
}

rankscorefun <- function(pmeasures,
                         pred_error,
                         is_normalized = TRUE){
  Performance_measures <- c("AUC_cipa_1", "AUC_cipa_2",
                            "LR_positive_cipa_1", "LR_positive_cipa_2",
                            "LR_negative_cipa_1", "LR_negative_cipa_2",
                            "Pairwise","Classification_error")
  Performance_levels <- c("Excellent_performance", 
                          "Good_performance",
                          "Minimally_acceptable_performance",
                          "Not_acceptabel")
  
  # A dataframe for the weights of performance measures
  pm_df <- data.frame(
    Performance_measure = Performance_measures,
    Weight = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0)
  )

  # A dataframe for the weights of performance levels
  pl_df <- data.frame(
    Performance_level = Performance_levels,
    Weight = c(3.0, 2.0, 1.0, 0.0)
  )
  if (is_normalized) {
    pm_df$Weight <- pm_df$Weight / sum(pm_df$Weight)
    pl_df$Weight <- pl_df$Weight / max(pl_df$Weight)
  }
  pl_df$Weight[4] <- NA # Not acceptable performance is removed
  
  # Initialize the performance dataframe for the model
  model_df <- data.frame(
    Performance_measure = Performance_measures,
    Performance_level_weight = c(NA, NA, NA, NA, NA, NA, NA, NA)
  )
  
  # Check the performance level for each performance measure 
  # by looking at the 95% confidence interval that match the CiPA's criteria
  
  # AUC_cipa_1
  score_to_check <- quantile(pmeasures$AUC_cipa_1, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # AUC_cipa_2
  score_to_check <- quantile(pmeasures$AUC_cipa_2, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_positive_cipa_1
  score_to_check <- quantile(pmeasures$LR_positive_cipa_1, 0.025)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_positive_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_positive_cipa_2
  score_to_check <- quantile(pmeasures$LR_positive_cipa_2, 0.025)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_positive_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_negative_cipa_1
  score_to_check <- quantile(pmeasures$LR_negative_cipa_1, 0.975)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_negative_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_negative_cipa_2
  score_to_check <- quantile(pmeasures$LR_negative_cipa_2, 0.975)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_negative_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Pairwise
  score_to_check <- quantile(pmeasures$Pairwise, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Pairwise",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Classification_error
  score_to_check <- mean(pred_error) + 1.96 * sd(pred_error) / sqrt(length(pred_error))
  if (score_to_check > 1.0) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 1.0 & score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.3) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Calculate rank score
  rank_score <- model_df$Performance_level_weight %*% pm_df$Weight
  
  return(as.numeric(rank_score))
}

rankscorefun_training <- function(pmeasures,
                         pred_error,
                         is_normalized = TRUE){
  Performance_measures <- c("AUC_cipa_1", "AUC_cipa_2",
                            "LR_positive_cipa_1", "LR_positive_cipa_2",
                            "LR_negative_cipa_1", "LR_negative_cipa_2",
                            "Classification_error")
  Performance_levels <- c("Excellent_performance", 
                          "Good_performance",
                          "Minimally_acceptable_performance",
                          "Not_acceptabel")
  
  # A dataframe for the weights of performance measures
  pm_df <- data.frame(
    Performance_measure = Performance_measures,
    Weight = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0)
  )
  
  # A dataframe for the weights of performance levels
  pl_df <- data.frame(
    Performance_level = Performance_levels,
    Weight = c(3.0, 2.0, 1.0, 0.0)
  )
  if (is_normalized) {
    pm_df$Weight <- pm_df$Weight / sum(pm_df$Weight)
    pl_df$Weight <- pl_df$Weight / max(pl_df$Weight)
  }
  pl_df$Weight[4] <- NA # Not acceptable performance is removed
  
  # Initialize the performance dataframe for the model
  model_df <- data.frame(
    Performance_measure = Performance_measures,
    Performance_level_weight = c(NA, NA, NA, NA, NA, NA, NA)
  )
  
  # Check the performance level for each performance measure 
  # by looking at the 95% confidence interval that match the CiPA's criteria
  
  # AUC_cipa_1
  score_to_check <- quantile(pmeasures$AUC_cipa_1, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # AUC_cipa_2
  score_to_check <- quantile(pmeasures$AUC_cipa_2, 0.025)
  if (score_to_check < 0.7) {
    model_df[model_df$Performance_measure == "AUC_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 0.7 & score_to_check < 0.8) {
    model_df[model_df$Performance_measure == "AUC_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 0.8 & score_to_check < 0.9) {
    model_df[model_df$Performance_measure == "AUC_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "AUC_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_positive_cipa_1
  score_to_check <- quantile(pmeasures$LR_positive_cipa_1, 0.025)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_positive_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_positive_cipa_2
  score_to_check <- quantile(pmeasures$LR_positive_cipa_2, 0.025)
  if (score_to_check < 2.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check >= 2.0 & score_to_check < 5.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check >= 5.0 & score_to_check < 10.0) {
    model_df[model_df$Performance_measure == "LR_positive_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_positive_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_negative_cipa_1
  score_to_check <- quantile(pmeasures$LR_negative_cipa_1, 0.975)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_negative_cipa_1",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # LR_negative_cipa_2
  score_to_check <- quantile(pmeasures$LR_negative_cipa_2, 0.975)
  if (score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.2) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.2 & score_to_check > 0.1) {
    model_df[model_df$Performance_measure == "LR_negative_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "LR_negative_cipa_2",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Classification_error
  score_to_check <- mean(pred_error) + 1.96 * sd(pred_error) / sqrt(length(pred_error))
  if (score_to_check > 1.0) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Not_acceptabel",]$Weight
  } else if (score_to_check <= 1.0 & score_to_check > 0.5) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Minimally_acceptable_performance",]$Weight
  } else if (score_to_check <= 0.5 & score_to_check > 0.3) {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Good_performance",]$Weight
  } else {
    model_df[model_df$Performance_measure == "Classification_error",]$Performance_level_weight <- pl_df[pl_df$Performance_level == "Excellent_performance",]$Weight
  }
  
  # Calculate rank score
  rank_score <- model_df$Performance_level_weight %*% pm_df$Weight
  
  return(as.numeric(rank_score))
}


fitandgetTMS <- function(training_data,testing_data,formula,dimension,max_attempts,logfile){
  converged <- FALSE
  
  # First trial without start option
  tryCatch({
    # Attempt to fit the model without specifying start
    mod <- polr(formula, data = training_data, Hess = TRUE)
    converged <- TRUE
  }, error = function(e) {
    write(sprintf("First trial without 'start' failed. Retrying..."),logfile)
  })
  
  # If the first trial fails, attempt to fit with random starting values
  attempts <- as.integer(0)
  set.seed(1)
  if (!converged) {
    for (attemp_idx in 1:max_attempts) {
      # Generate random starting values
      start_values <- rnorm(dimension + 2, mean = 0.0, sd = 100.0)
      
      # Attempt to fit the model with random starting values
      tryCatch({
        mod <- polr(formula, data = training_data, start = start_values, Hess = TRUE)
        converged <- TRUE
        break
      }, error = function(e) {
        attempts <- attempts + 1
      })
    }
  }
  
  # If all attempts fail, print a message
  if (!converged) {
    write(sprintf("Could not fit the model after %d attempts",max_attempts),logfile)
    return(
      list(
        training_data = training_data,
        testing_data = testing_data,
        mod = NA,
        converged = converged
      )
    )
    
  } else {
    # Make predictions for training dataset
    predicted_labels <- predict(mod, newdata = training_data, type = "class")
    
    # Variables from Ordinal Logistic Regression model
    alphas <- as.numeric(mod$zeta)  # Alpha values
    betas <- as.numeric(mod$coefficients)  # Beta coefficients for all features
    
    # if (!is_single) {
      # Apply TMS function row-wise efficiently
      training_data$TMS <- apply(training_data[, 1:dimension], 1, function(row) TMS(alphas, betas, row))
      testing_data$TMS <- apply(testing_data[, 1:dimension], 1, function(row) TMS(alphas, betas, row))
    # }
    return(
      list(
        training_data = training_data,
        testing_data = testing_data,
        mod = mod,
        converged = converged
      )
    )
  }
}

returnNA <- function(log_filename,dimension = NA){
  summarydf <- data.frame(
    Feature_Pair = log_filename,
    Accuracy_cipa_1 = NA,
    Accuracy_cipa_2 = NA,
    AUC_cipa_1 = NA,
    AUC_cipa_2 = NA,
    Sensitivity_cipa_1 = NA,
    Sensitivity_cipa_2 = NA,
    Specificity_cipa_1 = NA,
    Specificity_cipa_2 = NA,
    LR_positive_cipa_1 = NA,
    LR_positive_cipa_2 = NA,
    LR_negative_cipa_1 = NA,
    LR_negative_cipa_2 = NA,
    F1score_cipa_1 = NA,
    F1score_cipa_2 = NA,
    Classification_error = NA,
    Pairwise_classification_accuracy = NA,
    Rank_score = NA
  )

  if (!is.na(dimension)) {
    # Add log likelihood value
    summarydf['logLik'] <- NA 
    
    # Add Alphas
    for (i in 1:2) {
      summarydf[[paste0("Alpha_", i)]] <- NA
    }
    
    # Add Betas
    for (i in 1:dimension) {
      summarydf[[paste0("Beta_", i)]] <- NA
    }
  }
  
  return (summarydf)
}
