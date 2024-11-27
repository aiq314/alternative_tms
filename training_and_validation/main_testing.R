library(readr)
library(MASS)
library(pROC)
library(ggplot2)
library(foreach)
library(doParallel)
source("functions.R")

# Declare features and units
features <- c("qNet",
              "dvdtmax",
              "vmax",
              "vrest",
              "APD50",
              "APD90",
              "max_dv",
              "camax",
              "carest",
              "CaTD50",
              "CaTD90")

units <- c("",
           "",
           "",
           "",
           "",
           "",
           "",
           "",
           "",
           "",
           "")

# Declare the file paths
filepath_training <- "data/nanion_training_scaled_w_tms_vrest-max_dv-camax.csv"
filepath_testing <- "data/nanion_testing_scaled_w_tms_vrest-max_dv-camax.csv"

# Set feature dimension
dimension <-11

# Create pairsdf with all unique combinations
pairsdf <- pairsdfinitfun(features = features, units = units, dimension = dimension)

# Create the results folder
results_folder <- "results_testing_nanion_11_latest"

# Choose whether data needs to be normalized
is_normalized <- FALSE

# Choose how many attemps required before skipping the fitting process (when divergence is found)
max_attempts <- 10000

# Choose how many tests required for evaluating model performance
num_tests <- 10000

# Register parallel backend
numCores <- 1

# Check if the folder exists
if (!dir.exists(results_folder)) {
  # The folder does not exist, so create it
  dir.create(results_folder)
  cat("Folder created:", results_folder, "\n")
} else {
  # The folder already exists
  # List all files in the folder
  files <- list.files(path = results_folder, full.names = TRUE)
  # Remove all files in the folder
  if (length(files) > 0) {
    file.remove(files)
    cat("All files in the folder have been removed.\n")
  } else {
    cat("The folder is already empty.\n")
  }
}

# Execute the tasks in parallel
cl <- makeCluster(numCores)
registerDoParallel(cl)

is_loocv <- FALSE
summarydf <- foreach( pair_id = 1:nrow(pairsdf),
                      .combine = 'rbind',
                      .packages = c("readr","MASS","pROC","ggplot2")) %dopar% {
                        # Select the row corresponding to pair_id
                        pair_row <- pairsdf[pair_id, ]
                        
                        # Extract vectors of features and units
                        features_vector <- pair_row[grep("feature_", names(pair_row))]
                        units_vector <- pair_row[grep("unit_", names(pair_row))]
                        
                        # Print the current pair_id and features being processed
                        print(paste(c(pair_id, features_vector), collapse = "_"))
                        
                        # Read in the training and testing datasets
                        training <- read_csv(filepath_training, show_col_types = FALSE)
                        testing <- read_csv(filepath_testing, show_col_types = FALSE)
                        
                        result <- run_all_testing (results_folder = results_folder,
                                                    training = training,
                                                    testing = testing,
                                                    features_vector = features_vector,
                                                    units_vector = units_vector,
                                                    is_normalized = is_normalized,
                                                    max_attempts = max_attempts,
                                                    num_tests = num_tests)
                        
                        # Return the result
                        result
                      }
write.csv(summarydf, paste(results_folder,"/","summary.csv", sep = ""), row.names = FALSE)

# Stop the parallel cluster
stopCluster(cl)