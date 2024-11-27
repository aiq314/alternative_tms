# File:         compute_TdP_error.R
# Author:       Kelly Chang
#               Zhihua Li
# Modifier:     AYD x LTF X ALI
# Date:         Nov 2017 (Nov 2023)
# Version:      1.0 (1.1)
# 
# Description:  R script to perform Torsade de Pointes (TdP) risk
#               classification using ordinal logistic regression and
#               leave-one-out cross validation (LOOCV).
#               For help, run this script with command line option "-h".
#

#--- specify command line arguments
library(optparse)
parser <- OptionParser()
#parser <- add_option(parser, c("-t", "--tdpfile"), default = "../chantest_AP_simulation/data/newCiPA.csv", help = "Filepath to drug list")
parser <- add_option(parser, c("-u", "--uncertainty"), default = TRUE, action = "store_true", help = "Flag to use simulations with uncertainty inputs for training and cross-validation")
parser <- add_option(parser, c("-o", "--omit"), default = FALSE, action = "store_true", help = "Flag to use drug simulations with ionic current effects omitted")
parser <- add_option(parser, c("-i", "--individual"), default = TRUE, action = "store_true", help = "Flag to output the classification for each individual")
parser <- add_option(parser, c("-m", "--metric"), default = "qNet", help = "the metric to use")
parser <- add_option(parser, c("-d", "--APpath"), default = "../chantest_AP_simulation/", help = "Filepath to the metrics.rds")
args <- parse_args(parser)

set.seed(100)

#--- load libraries
library(rms)
require(ROCR)
library(ggplot2)
library(comprehenr)
library(MASS)

source("drug_combination.R")
print(sessionInfo())

#--- get arguments
#tdpfile <- args$tdpfile
useUQ <- args$uncertainty
omit <- args$omit
outputI <- args$individual
# metric <- args$metric
metric <- "tms"
APpath <- args$APpath
scale <- "sara"
CL <- 2000
metricv <- c("tms")  
datasetvec <- c("chantest")  

# setup output directory
input <- ifelse(useUQ, "uncertainty", "fixed")
simdir <- sprintf("results_6/%s/", input)
if (omit) {
    outdir_vec <- Sys.glob(paste0(simdir, "drop_*/"))
} else {
    outdir_vec <- simdir
}

#--- function to try ordinal logistic regression
try_lrm <- function(datadf, tol = 1e-10, maxit = 1e6) {
  # Create and assign datadist object
  dd <- datadist(datadf)
  options(datadist = "dd")
  try({ lrm(CiPA ~ get(metric), data = datadf, penalty = 0, x = TRUE, y = TRUE, tol = tol, maxit = maxit) })
}

#--- initialize dictionary of drugs
drug_risks <- c(c("mexiletine"=0,"loratadine"=0,"metoprolol"=0,"ranolazine"=0,"diltiazem"=0,"tamoxifen"=0,"nitrendipine"=0,
                  "droperidol"=1,"cisapride"=1,"risperidone"=1,"astemizole"=1,"clozapine"=1,"domperidone"=1,
                  "sotalol"=2,"dofetilide"=2,"disopyramide"=2,"bepridil"=2,"vandetanib"=2))

# Define measures
measurevec <- c("AUC1","AUC2","Pairwise","LR1plus","LR1minus","LR2plus","LR2minus","Mean_error",
              "AUC1_LOOCV","AUC2_LOOCV","Pairwise_LOOCV","LR1plus_LOOCV","LR1minus_LOOCV","LR2plus_LOOCV","LR2minus_LOOCV","Mean_error_LOOCV",
              "AUC1_training","AUC2_training","Pairwise_training","LR1plus_training","LR1minus_training","LR2plus_training","LR2minus_training","Mean_error_training",
              "Threshold1", "Threshold2", "Normalized_logLik")

#--- creating drug combination
#TODO: Adjust combine_drug value based on number of drugs
drug_size = 6

columns <- c(to_vec(for(i in 1:drug_size) paste0("drug",i)), measurevec) 
result_df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(result_df) <- columns

drug_combination <- combine_drug(drug_size)
# for(drug_index in seq_len(nrow(drug_combination))) {
for(drug_index in 1:3) {
  print(paste0("-------------------------------------", drug_index))  
  drugtable <- data.frame(drug = c(drug_combination[drug_index,]),
                          CiPA = c(to_vec(for(i in drug_combination[drug_index,]) drug_risks[i])))
  drugnames <- as.character(drugtable$drug)
  
  #--- do logistic regression
  for (outdir in outdir_vec) {
      # read in dataset
      # if (metric == "qNet") {
          # infile <- paste0(APpath, outdir, "metrics.rds")
          infile <- paste0(APpath, outdir, "metrics_nanion_scaled_w_tms_dvdtmax-APD90-camax-carest.csv")
          df <- read.csv(infile)
          df <- df[df$drug != "control" & df$dose != 0, ]
          df <- df[!is.na(df$max_dv), ] # depolarization failures
          df <- df[(df$dose) %in% c("1", "2", "3", "4"), ]                  #very important cause some drugs may have additional doses in metrics.rds!!!!
          if (input == "uncertainty") {
              qNettable <- aggregate(df[, metric], list(df$drug, df$sample), mean)  #name is qNettable, but it could be any metric
              colnames(qNettable) <- c("drug", "sample", metric)
          } else {
              qNettable <- aggregate(df[, metric], list(df$drug), mean)
              colnames(qNettable) <- c("drug", "qNet")
              qNettable$sample <- NA
          }
          # print(head(qNettable))
          # print(class(qNettable$qNet))
          if (metric == "qNet")
              qNettable$qNet <- qNettable$qNet / 1000
          df <- merge(qNettable, drugtable[, c("drug", "CiPA")], by.x = "drug", by.y = "drug")
      # }
   
      #organize df
      df$drug <- factor(df$drug, levels = drugnames)
      df <- df[order(df$drug, df$sample), ]
      df$CiPA <- ordered(df$CiPA)
  
      # perform logistic regression and leave-one-out cross validation for each drug
      errdf <- data.frame()
      probdf <- data.frame()
      allprobdf <- data.frame()   #corresponding to cvdf but without LOOCV
      cverrdf <- data.frame()
      cvprobdf <- data.frame()
      cvdf <- data.frame()                              #used to store individual sample's probs
      is_training_all_error <- FALSE
  
      for (drug in c(NA, drugnames)) {
  
          # fit logistic regression model
          if (is.na(drug)) {
              print("training on all drugs")
              datadf <- df[, c("drug", "CiPA", "sample", metric)]
              testdf <- datadf
              traindf <- datadf
          } else {
              print(sprintf("cross validating with %s", drug))
              testdf <- datadf[datadf$drug == drug, ]
              traindf <- datadf[datadf$drug != drug, ]
              if (nrow(testdf) == 0) {
                  print("no data to test! skipping...")
                  next
              }
          }
  
          # lmod <- try_lrm(traindf)
        dd <- datadist(traindf)
        options(datadist = "dd")
        lmod <- try({ lrm(CiPA ~ get(metric), data = traindf, penalty = 0, x = TRUE, y = TRUE, tol = 1e-10, maxit = 1e6) })
        
          # print(lmod)
          if (inherits(lmod, "try-error")) {
              if (is.na(drug)) {
                  print("fitting all drugs failed? skipping...")
                  is_training_all_error <- TRUE
                  break
              } else {
                  print(sprintf("cross validating with %s failed! skipping...", drug))
                  next
              }
          }
  
          # save coefficients
          print(sprintf("Convergence failure: %s", lmod$fail))
          cf <- coefficients(lmod)
          cfvec <- c()
          for (kint in 1:(length(cf) - 1))
              cfvec[[paste0("intercept", kint)]] <- cf[[kint]]
          cfvec[["slope"]] <- cf[[length(cf)]]
          
          #use formal math for threshold 1 and 2
          ebeta1<-exp(-cfvec[["intercept1"]]); ebeta2<- exp(-cfvec[["intercept2"]])
          t1 <- log(ebeta1*ebeta2/-(2*ebeta1-ebeta2))/cfvec[["slope"]]
          t2 <- log(exp(-cfvec[["intercept2"]])-2*exp(-cfvec[["intercept1"]]))/cfvec[["slope"]]
          
          # get training/prediction error
          y0 <- testdf$CiPA
          if (is.na(drug)) {
              probs <- predict(lmod, type = "fitted.ind")
          } else {
              probs <- matrix(predict(lmod, newdata = testdf, type = "fitted.ind"), ncol = length(levels(y0)))
          }
          
          #sometimes APD90 is NA, giving rise to NA in probs
          idx <- apply(probs, 1, function(x) any(is.na(x)))
          probs[idx, ] <- matrix(rep(c(0, 0, 1), sum(idx)), nrow = sum(idx), byrow = T) #this is to assume all NAs are due to repolarization failure, and thus high risk
          
          yPred <- apply(probs, 1, function(x) which.max(x))
          pred_err <- mean(abs(yPred - as.integer(y0)))
          
          # Loglikelihood value of the model
          logLikVal <- numeric(length(y0))
          
          for (i in 1:length(y0)) {
            logLikVal[i] <- log(probs[i,y0[i]])
          }
          
          # Normalized loglikelihood value
          nrows_trainingdf <- nrow(traindf)
          normalized_logLik <- logLikVal / nrows_trainingdf
        
          # append to data frame
          if (is.na(drug)) {
              newrow <- data.frame(t(cfvec), error = pred_err)
              errdf <- rbind(errdf, newrow)
              
              # allprobdf <- cbind(traindf[, 1:3], probs)
              # colnames(allprobdf) <- c("drug", "CiPA", "sample", "low_prob", "inter_prob", "high_prob")
              allprobdf <- cbind(traindf[, 1:3], probs, abs(yPred - as.integer(y0)), normalized_logLik)
              colnames(allprobdf) <- c("drug", "CiPA", "sample", "low_prob", "inter_prob", "high_prob", "pred_err", "normalized_logLik")
             
          } else {
              newrow <- data.frame(drug = drug, t(cfvec), error = pred_err)
              cverrdf <- rbind(cverrdf, newrow)
              if (outputI) {                            #make detailed prob table here
                  testdf <- cbind(testdf, high_prob = probs[, 3], inter_prob = probs[, 2], low_prob = probs[, 1], pred_err = abs(yPred - as.integer(y0)))
                  cvdf <- rbind(cvdf, testdf)
              } #if outputI
          }
  
          # detailed results
          newcols <- sapply(levels(y0), function(s) paste0("predict_", s), USE.NAMES = F)
          predictions <- factor(newcols[yPred], levels = newcols)
          tmpdf <- as.data.frame(tapply(1:nrow(testdf), list(drug = as.character(testdf$drug), predict = predictions), FUN = length))
          tmpdf[is.na(tmpdf)] <- 0
          tmpdf$drug <- factor(rownames(tmpdf), levels = drugnames)
          
          tmpdf <- merge(drugtable[, c("drug", "CiPA")], tmpdf)
          tmpdf <- tmpdf[order(tmpdf$drug), ]
          rownames(tmpdf) <- NULL
          if (is.na(drug)) {
              probdf <- rbind(probdf, tmpdf)
             
          } else {
              cvprobdf <- rbind(cvprobdf, tmpdf)
          }
  
      } # for drug
   
      if (!is_training_all_error) {
        outfile <- paste0(outdir, metric, "_", drug_index, "_training_probs.csv")   #here training means "trained on all drugs"
        write.csv(probdf, outfile, row.names = F, quote = F)
        # print(head(probdf))
        
        outfile <- paste0(outdir, metric,  "_", drug_index, "_training_errors.csv")
        list_columns <- sapply(errdf, is.list) # Identify columns in errdf that are of type list
        if (any(list_columns)) { # Convert list columns to character (or other suitable type) if any
          errdf[list_columns] <- lapply(errdf[list_columns], as.numeric)
        }
        write.csv(errdf, outfile, row.names = F, quote = F)
        # print(head(errdf))
        
        outfile <- paste0(outdir, metric,  "_", drug_index, "_LOOCV_probs.csv")
        write.csv(cvprobdf, outfile, row.names = F, quote = F)
        # print(head(cvprobdf))
        
        outfile <- paste0(outdir, metric,  "_", drug_index, "_LOOCV_errors.csv")
        list_columns <- sapply(cverrdf, is.list) # Identify columns in errdf that are of type list
        if (any(list_columns)) { # Convert list columns to character (or other suitable type) if any
          cverrdf[list_columns] <- lapply(cverrdf[list_columns], as.numeric)
        }
        write.csv(cverrdf, outfile, row.names = F, quote = F)
        # print(head(cverrdf))
        
        if (outputI) {
          outfile <- paste0(outdir, metric,  "_", drug_index, "_LOOCV_allprobes.csv")
          write.csv(cvdf, outfile, row.names = F, quote = F)
          # print(head(cvdf))
          outfile <- paste0(outdir, metric,  "_", drug_index, "_training_allprobes.csv")
          write.csv(allprobdf, outfile, row.names = F, quote = F)
        }
      }
  } # for outdir
  
  if (!is_training_all_error) {
    
    # Initialize data frames
    # all28df <- data.frame(metric=c("AUC1","AUC2","Pairwise","LR1plus","LR1minus","LR2plus","LR2minus","Mean_error"))
    all28df <- data.frame(metric=c("AUC1","AUC2","Pairwise","LR1plus","LR1minus","LR2plus","LR2minus","Mean_error",
                                   "AUC1_LOOCV","AUC2_LOOCV","Pairwise_LOOCV","LR1plus_LOOCV","LR1minus_LOOCV","LR2plus_LOOCV","LR2minus_LOOCV","Mean_error_LOOCV",
                                   "AUC1_training","AUC2_training","Pairwise_training","LR1plus_training","LR1minus_training","LR2plus_training","LR2minus_training","Mean_error_training",
                                   "Threshold1", "Threshold2", "Normalized_logLik"))
    colordf <- data.frame(
      drug = c("dofetilide","bepridil","quinidine","sotalol","ibutilide", "azimilide", "disopyramide", "vandetanib",
               "cisapride","terfenadine","ondansetron","chlorpromazine","clarithromycin","risperidone","domperidone","astemizole", "pimozide","droperidol","clozapine",
               "ranolazine","mexiletine","diltiazem","verapamil","metoprolol","loratadine","tamoxifen","nifedipine","nitrendipine"),
      classidx=c(rep(2,8),rep(1,11),rep(0,9)),
      coloridx=c(rep(1,8),rep(3,11),rep(2,9)),
      isTraining=c(1,1,1,1,0,0,0,0,
                   1,1,1,1,0,0,0,0,0,0,0,
                   1,1,1,1,0,0,0,0,0)
    )
    
    # Loop over datasets and metrics
    for(dataset in datasetvec){
      for(metric in metricv){
        altstr <- paste0(dataset,"_",metric)
        altvec <- 0
        valdf <- read.delim(paste0("results_6/uncertainty/",metric,"_", drug_index, "_LOOCV_allprobes.csv"),sep=",",as.is=T)
        trainingdf <- read.delim(paste0("results_6/uncertainty/",metric,"_", drug_index, "_training_allprobes.csv"),sep=",",as.is=T)
        file.remove(paste0("results_6/uncertainty/",metric,"_", drug_index, "_LOOCV_allprobes.csv"))
        file.remove(paste0("results_6/uncertainty/",metric,"_", drug_index, "_training_allprobes.csv"))
        
        # Classify data
        valdf$classidx <- valdf$CiPA
        valdf$class <- 2 - valdf$classidx 
        
        trainingdf$classidx <- trainingdf$CiPA
        trainingdf$class <- 2 - trainingdf$classidx 
        
        # Prepare data for ROC analysis
        forROC1 <- valdf[,c("inter_prob","classidx","sample", "drug")]
        forROC1$inter_prob <- rowSums(valdf[,c("inter_prob","high_prob")])
        colnames(forROC1)[1] <- "positive_prob"
        forROC1$class <- forROC1$classidx > 0  # now high and intermediate are both "TRUE" or "positive" 
        
        forROC1_training <- trainingdf[,c("inter_prob","classidx","sample", "drug")]
        forROC1_training$inter_prob <- rowSums(trainingdf[,c("inter_prob","high_prob")])
        colnames(forROC1_training)[1] <- "positive_prob"
        forROC1_training$class <- forROC1_training$classidx > 0  # now high and intermediate are both "TRUE" or "positive" 
        
        # Generate predictions and labels
        forROC1predictions <- do.call(cbind,by(forROC1[,"positive_prob"],forROC1$sample,function(x) x))
        forROC1labels <- do.call(cbind,by(forROC1$class,forROC1$sample,function(x) x))
        
        forROC1predictions_training <- do.call(cbind,by(forROC1_training[,"positive_prob"],forROC1_training$sample,function(x) x))
        forROC1labels_training <- do.call(cbind,by(forROC1_training$class,forROC1_training$sample,function(x) x))
        
        # Perform ROC analysis
        ROC1obj <- prediction(forROC1predictions, forROC1labels)        
        AUC1 <- performance(ROC1obj,"auc")  #the value is in y.values
        AUC1dist <- formatC(quantile(unlist(AUC1@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(AUC1dist[2]," (",AUC1dist[1]," - ",AUC1dist[3],")"))
        
        ROC1obj_training <- prediction(forROC1predictions_training, forROC1labels_training)        
        AUC1_training <- performance(ROC1obj_training,"auc")  #the value is in y.values
        AUC1dist_training <- formatC(quantile(unlist(AUC1_training@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(AUC1dist_training[2]," (",AUC1dist_training[1]," - ",AUC1dist_training[3],")"))
        
        # Prepare data for second ROC analysis
        forROC2 <- valdf[,c("high_prob","classidx","sample","drug")]
        colnames(forROC2)[1] <- "positive_prob"
        forROC2$class <- forROC2$classidx > 1                         
        
        forROC2_training <- trainingdf[,c("high_prob","classidx","sample","drug")]
        colnames(forROC2_training)[1] <- "positive_prob"
        forROC2_training$class <- forROC2_training$classidx > 1
        
        # Generate predictions and labels
        forROC2predictions <- do.call(cbind,by(forROC2[,"positive_prob"],forROC2$sample,function(x) x))  
        forROC2labels <- do.call(cbind,by(forROC2$class,forROC2$sample,function(x) x))
        
        forROC2predictions_training <- do.call(cbind,by(forROC2_training[,"positive_prob"],forROC2_training$sample,function(x) x))  
        forROC2labels_training <- do.call(cbind,by(forROC2_training$class,forROC2_training$sample,function(x) x))
        
        # Perform second ROC analysis
        ROC2obj <- prediction(forROC2predictions, forROC2labels)        
        AUC2 <- performance(ROC2obj,"auc")  #the value is in y.values 
        AUC2dist <- formatC(quantile(unlist(AUC2@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(AUC2dist[2]," (",AUC2dist[1]," - ",AUC2dist[3],")"))
        
        ROC2obj_training <- prediction(forROC2predictions_training, forROC2labels_training)        
        AUC2_training <- performance(ROC2obj_training,"auc")  #the value is in y.values 
        AUC2dist_training <- formatC(quantile(unlist(AUC2_training@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(AUC2dist_training[2]," (",AUC2dist_training[1]," - ",AUC2dist_training[3],")"))
        
        # Pairwise ranking
        pairwisefun <- function(fulltable){
          cmb <- combn(seq_len(nrow(fulltable)), 2)
          mergedtable <- cbind(fulltable[cmb[1,],], fulltable[cmb[2,],])
          validpairidx <- (mergedtable[,7]!=mergedtable[,16])&(!mergedtable[,9]|!mergedtable[,18])
          correctidx1 <- ((mergedtable[,7]>mergedtable[,16])&(mergedtable[,6]<mergedtable[,15]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,6]>mergedtable[,15])) #when predicted class are different
          correctidx2 <- (mergedtable[,6]==1)&(mergedtable[,15]==1)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,3]>mergedtable[,12]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,3]<mergedtable[,12]))) #when predicted class are both high
          correctidx3 <- (mergedtable[,6]==3)&(mergedtable[,15]==3)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,5]<mergedtable[,14]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,5]>mergedtable[,14]))) #when predicted class are both low
          correctidx4 <- (mergedtable[,6]==2)&(mergedtable[,15]==2)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,5]<mergedtable[,14]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,5]>mergedtable[,14]))) #when predicted class are both intermediate
          correctidx <- correctidx1|correctidx2|correctidx3|correctidx4
          sum(validpairidx&correctidx)/sum(validpairidx)
        }
        tempdf <- valdf[,c("drug","sample","high_prob","inter_prob","low_prob")]; 
        tempdf$pred <- apply(valdf[,5:7],1,which.max) #note here 1 is high, 2 is int, 3 is low
        tempdf <- merge(tempdf,colordf,by="drug")
        distpairwise <- by(tempdf, valdf$sample, pairwisefun)
        pairwisedist <- formatC(quantile(unlist(distpairwise),prob=c(0.025,0.5,0.975),na.rm=TRUE),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(pairwisedist[2]," (",pairwisedist[1]," - ",pairwisedist[3],")"))
        
        tempdf_training <- trainingdf[,c("drug","sample","high_prob","inter_prob","low_prob")]; 
        tempdf_training$pred <- apply(trainingdf[,5:7],1,which.max) #note here 1 is high, 2 is int, 3 is low
        tempdf_training <- merge(tempdf_training,colordf,by="drug")
        distpairwise_training <- by(tempdf_training, trainingdf$sample, pairwisefun)
        pairwisedist_training <- formatC(quantile(unlist(distpairwise_training),prob=c(0.025,0.5,0.975),na.rm=TRUE),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(pairwisedist_training[2]," (",pairwisedist_training[1]," - ",pairwisedist_training[3],")"))
        
        # Calculate sensitivity and specificity
        sens1 <- with(valdf, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0))
        spec1 <-  with(valdf, sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0))
        LRplus1 <- sens1/(1-spec1)
        LRminus1 <- (1-sens1)/spec1
        
        sens2 <- with(valdf, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2))
        spec2 <-  with(valdf, sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2))
        LRplus2 <- sens2/(1-spec2)
        LRminus2 <- (1-sens2)/spec2
        
        sens1_training <- with(trainingdf, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0))
        spec1_training <-  with(trainingdf, sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0))
        LRplus1_training <- sens1_training/(1-spec1_training)
        LRminus1_training <- (1-sens1_training)/spec1_training
        
        sens2_training <- with(trainingdf, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2))
        spec2_training <-  with(trainingdf, sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2))
        LRplus2_training <- sens2_training/(1-spec2_training)
        LRminus2_training <- (1-sens2_training)/spec2_training
        
        # Calculate LR using random sampling
        sens1UQ <- by(valdf, valdf$sample, function(x) with(x, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0)))
        spec1UQ <- by(valdf, valdf$sample, function(x) with(x,sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0)))
        LRplus1UQ <- (sens1UQ+rnorm(length(sens1UQ),1e-6,1e-12))/(1-spec1UQ+rnorm(length(sens1UQ),1e-6,1e-12))
        LRplus1dist <- formatC(quantile(unlist(LRplus1UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRplus1dist[2]," (",LRplus1dist[1]," - ",LRplus1dist[3],")"))
        LRminus1UQ <- (1-sens1UQ+rnorm(length(sens1UQ),1e-6,1e-12))/(spec1UQ+rnorm(length(sens1UQ),1e-12,1e-12))
        LRminus1dist <- formatC(quantile(unlist(LRminus1UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRminus1dist[2]," (",LRminus1dist[1]," - ",LRminus1dist[3],")"))
        
        sens2UQ <- by(valdf, valdf$sample, function(x) with(x, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2)))
        spec2UQ <- by(valdf, valdf$sample, function(x) with(x,sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2)))
        LRplus2UQ <- (sens2UQ+rnorm(length(sens2UQ),1e-6,1e-12))/(1-spec2UQ+rnorm(length(sens2UQ),1e-6,1e-12))
        LRplus2dist <- formatC(quantile(unlist(LRplus2UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRplus2dist[2]," (",LRplus2dist[1]," - ",LRplus2dist[3],")"))
        LRminus2UQ <- (1-sens2UQ+rnorm(length(sens2UQ),1e-6,1e-12))/(spec2UQ+rnorm(length(sens2UQ),1e-6,1e-12))
        LRminus2dist <- formatC(quantile(unlist(LRminus2UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRminus2dist[2]," (",LRminus2dist[1]," - ",LRminus2dist[3],")"))
        
        sens1UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0)))
        spec1UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x,sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0)))
        LRplus1UQ_training <- (sens1UQ_training+rnorm(length(sens1UQ_training),1e-6,1e-12))/(1-spec1UQ_training+rnorm(length(sens1UQ_training),1e-6,1e-12))
        LRplus1dist_training <- formatC(quantile(unlist(LRplus1UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRplus1dist_training[2]," (",LRplus1dist_training[1]," - ",LRplus1dist_training[3],")"))
        LRminus1UQ_training <- (1-sens1UQ_training+rnorm(length(sens1UQ_training),1e-6,1e-12))/(spec1UQ_training+rnorm(length(sens1UQ_training),1e-12,1e-12))
        LRminus1dist_training <- formatC(quantile(unlist(LRminus1UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRminus1dist_training[2]," (",LRminus1dist_training[1]," - ",LRminus1dist_training[3],")"))
        
        sens2UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2)))
        spec2UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x,sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2)))
        LRplus2UQ_training <- (sens2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))/(1-spec2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))
        LRplus2dist_training <- formatC(quantile(unlist(LRplus2UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRplus2dist_training[2]," (",LRplus2dist_training[1]," - ",LRplus2dist_training[3],")"))
        LRminus2UQ_training <- (1-sens2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))/(spec2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))
        LRminus2dist_training <- formatC(quantile(unlist(LRminus2UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRminus2dist_training[2]," (",LRminus2dist_training[1]," - ",LRminus2dist_training[3],")"))
        
        # Calculate mean error
        meanerror <- mean(valdf$pred_err)
        lowererror <- meanerror - 1.96* sd(valdf$pred_err)/sqrt(dim(valdf)[1]); uppererror <- meanerror+1.96*sd(valdf$pred_err)/sqrt(dim(valdf)[1])
        errordist <- formatC(c(meanerror,lowererror,uppererror),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(errordist[1]," (",errordist[2]," - ",errordist[3],")"))
        
        meanerror_training <- mean(trainingdf$pred_err)
        lowererror_training <- meanerror_training - 1.96* sd(trainingdf$pred_err)/sqrt(dim(trainingdf)[1]); uppererror_training <- meanerror_training+1.96*sd(trainingdf$pred_err)/sqrt(dim(trainingdf)[1])
        errordist_training <- formatC(c(meanerror_training,lowererror_training,uppererror_training),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(errordist_training[1]," (",errordist_training[2]," - ",errordist_training[3],")"))
        
        # LOOCV results for ranking
        altvec <- c(altvec,paste0(AUC1dist[2]))
        altvec <- c(altvec,paste0(AUC2dist[2]))
        altvec <- c(altvec,paste0(pairwisedist[2]))
        altvec <- c(altvec,paste0(LRplus1dist[2]))
        altvec <- c(altvec,paste0(LRminus1dist[2]))
        altvec <- c(altvec,paste0(LRplus2dist[2]))
        altvec <- c(altvec,paste0(LRminus2dist[2]))
        altvec <- c(altvec,paste0(errordist[1]))
        
        # altvec <- c(altvec,paste0(AUC1dist[1]))
        # altvec <- c(altvec,paste0(AUC2dist[1]))
        # altvec <- c(altvec,paste0(pairwisedist[1]))
        # altvec <- c(altvec,paste0(LRplus1dist[1]))
        # altvec <- c(altvec,paste0(LRminus1dist[3]))
        # altvec <- c(altvec,paste0(LRplus2dist[1]))
        # altvec <- c(altvec,paste0(LRminus2dist[3]))
        # altvec <- c(altvec,paste0(errordist[3]))
        
        
        # Training results for ranking
        altvec <- c(altvec,paste0( AUC1dist_training[2]))
        altvec <- c(altvec,paste0( AUC2dist_training[2]))
        altvec <- c(altvec,paste0(pairwisedist_training[2]))
        altvec <- c(altvec,paste0(LRplus1dist_training[2]))
        altvec <- c(altvec,paste0(LRminus1dist_training[2]))
        altvec <- c(altvec,paste0(LRplus2dist_training[2]))
        altvec <- c(altvec,paste0(LRminus2dist_training[2]))
        altvec <- c(altvec,paste0(errordist_training[1]))
        
        # altvec <- c(altvec,paste0( AUC1dist_training[1]))
        # altvec <- c(altvec,paste0( AUC2dist_training[1]))
        # altvec <- c(altvec,paste0(pairwisedist_training[1]))
        # altvec <- c(altvec,paste0(LRplus1dist_training[1]))
        # altvec <- c(altvec,paste0(LRminus1dist_training[3]))
        # altvec <- c(altvec,paste0(LRplus2dist_training[1]))
        # altvec <- c(altvec,paste0(LRminus2dist_training[3]))
        # altvec <- c(altvec,paste0(errordist_training[3]))
        
        # Thresholds
        altvec <- c(altvec,paste0(t1))
        altvec <- c(altvec,paste0(t2))
        
        # The normalized loglikelihood value
        altvec <- c(altvec,paste0(mean(normalized_logLik)))
        
        # Update altvec and all28df
        altvec <- altvec[-1]
        all28df[,altstr] <- altvec
      }
    }
    
    # Write results to file
    
    table7 <- data.frame(metric=rep(measurevec,each=length(datasetvec)), 
                         dataset=rep(datasetvec, length(measurevec)),
                         qNet=NA)
    table7[,-1:-2] <- matrix(unlist(t(all28df[,-1])),nrow= length(measurevec)*length(datasetvec),byrow=T)
    outfile <- paste0("Result_PlotTableLOOCV_6/",drug_index,"_Table7.txt")
    write.table(table7, outfile, col.names=T,row.names=F,sep="\t",quote=F)
    table7 <- subset(table7, select=-c(dataset))
    #table7 <- subset(table7, select=-c(1))
    names(table7) <- NULL
    table7 <- setNames(data.frame(t(table7[,-1])), table7[,1])
    drugnames <- t(drugnames)
    colnames(drugnames) <- c(to_vec(for(i in 1:drug_size) paste0("drug_",i)))
    
    final_result <- cbind(drugnames, table7[1,])
    #print(final_result)
    result_df <- rbind(result_df, final_result)
    #print(result_df)
    outfile <- paste0(outdir, metric, "metrics.csv")
    write.csv(result_df, outfile ,row.names=F)
  } #if is_training_all_error
}

