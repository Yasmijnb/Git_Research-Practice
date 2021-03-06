###############################################################################

# Yasmijn Balder
# 02-03-2021

# Use different thresholds for young and old to find which gives the most
# accurate devision.

# Output
# The accuracy, sentivitiy, specificity, and AUC of the model

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/Lipids_age_sex.csv", check.names = TRUE)
data[,1] <- NULL

###############################################################################

# Packages
library(randomForest)
library(Rmisc)

# Functions
RForest = function(x.data, y.class, unbalance = TRUE, perc.unbal = 0.85, 
                   n.resamp = 100, max.perm = 1000, verbose = T){
  
  # Inputs:
  #
  #       x.data            data matrix
  #       y.class           vector of class labels
  #       perc.unbal        % of data to retain for each group in case of unblanced data. Default 85%
  #       n.resamp          Number of resampling for the unbalance. Default is 100.
  #       max.perm          Number of permutaion. Default is 1000
  #
  #       Edoardo Saccenti 10/2020 v 0.1
  
  #Number of samples
  n.samples =  nrow(x.data)
  
  # initialize result matrices
  ACCURACY.perm = matrix(NA, nrow = n.resamp, ncol = 1)
  SPECIFICITY.perm = matrix(NA, nrow = n.resamp, ncol = 1)
  SENSITIVITY.perm = matrix(NA, nrow = n.resamp, ncol = 1)
  AUC.perm = matrix(NA, nrow = n.resamp, ncol = 1)
  
  # UNBALANCE
  if(isTRUE(unbalance)){
    STAT = matrix(NA,4,4)
    print('Correction and resampling for unbalance performed')
    
    # Find unique class values
    class.values = unique(y.class)
    print(class.values)
    
    # Check number of classes. Must be 2
    if(length(class.values)> 2)
      stop('More than two classes found! Only two class discrimination problem is supported!')
    if(length(class.values)< 2)
      stop('Only one class found!')
    
    # Trick to take into account unbalance
    n1 = sum(y.class == class.values[1]); 
    n2 = sum(y.class == class.values[2]);
    nsize = round(min(n1,n2)*perc.unbal) #Use 85% of sample per classes
    STAT = matrix(NA,4,4)
    ACCURACY = matrix(NA, nrow = n.resamp, ncol = 1)
    SPECIFICITY = matrix(NA, nrow = n.resamp, ncol = 1)
    SENSITIVITY = matrix(NA, nrow = n.resamp, ncol = 1)
    AUC = matrix(NA, nrow = n.resamp, ncol = 1)
    
    # Random Forest on original data
    for(k in 1 : n.resamp){
      rf <- randomForest(as.factor(y.class) ~ ., data=x.data, strata = y.class, 
                         sampsize = c(nsize,nsize), replace = F)
      print('done')
      # Extract predicted y
      y.pred = as.factor(rf$predicted)
      
      # xtab <- table(as.numeric(y.pred), as.numeric(y.class))
      # conf.mat = confusionMatrix(xtab)
      
      conf.mat = confusionMatrix(as.factor(y.pred), as.factor(y.class))
      
      acc = conf.mat$overall["Accuracy"]
      sens = conf.mat$byClass["Sensitivity"]
      spec = conf.mat$byClass["Specificity"]
      
      # Calculate AUC
      AUC[k] = colAUC(rf$votes[,2],y.class,plotROC=F)
      ACCURACY[k] = acc
      SPECIFICITY[k] = spec
      SENSITIVITY[k] = sens
    }
    
    print(sprintf('RF discrimination model fitted'))
    
    # confidence intervals
    cis.acc = CI(ACCURACY, ci = 0.95);
    cis.sens = CI(SENSITIVITY, ci = 0.95);
    cis.spec = CI(SPECIFICITY, ci = 0.95);
    cis.auc = CI(AUC, ci = 0.95);
    
    # average specs
    average.acc = cis.acc[2]
    average.sens = cis.sens[2]
    average.spec = cis.spec[2]
    average.auc = cis.auc[2]
    
    #STAT
    STAT[1,1] = average.acc
    STAT[2,1] = average.spec
    STAT[3,1] = average.sens
    STAT[4,1] = average.auc
    
    STAT[1,2] = cis.acc[3]
    STAT[2,2] = cis.spec[3]
    STAT[3,2] = cis.sens[3]
    STAT[4,2] = cis.auc[3]
    
    STAT[1,3] = cis.acc[1]
    STAT[2,3] = cis.spec[1]
    STAT[3,3] = cis.sens[1]
    STAT[4,3] = cis.auc[1] 
    
    #### Run a permutation test ####
    
    for(i in 1 : max.perm ){
      if(isTRUE(verbose))
        print(sprintf('Permutation %i of %i',i,max.perm ))
      #Pemute class labels
      idx.p = sample(n.samples)
      y.class.perm = y.class[idx.p]
      
      acc.perm = NULL
      sens.perm = NULL
      spec.perm = NULL
      auc.perm = NULL
      
      for(k in 1 : n.resamp){
        rf.perm <- randomForest(as.factor(y.class.perm) ~ ., data= x.data, 
                                strata = y.class.perm, 
                                sampsize = c(nsize,nsize), replace=F)
        y.pred.perm = as.factor(rf.perm$predicted)
        
        #xtab.perm <- table(as.numeric(y.pred.perm), as.numeric(y.class.perm))
        
        conf.mat.perm = confusionMatrix(as.factor(y.pred.perm), 
                                        as.factor(y.class.perm))
        
        acc.perm[k] = conf.mat.perm$overall["Accuracy"]
        sens.perm[k] = conf.mat.perm$byClass["Sensitivity"]
        spec.perm[k] = conf.mat.perm$byClass["Specificity"]
        
        #Calculate AUC
        auc.perm[k] = colAUC(rf.perm$votes[,2],y.class.perm,plotROC=F)
      }
      ACCURACY.perm[i] = mean(acc.perm)
      SENSITIVITY.perm[i] = mean(sens.perm)
      SPECIFICITY.perm[i] = mean(spec.perm)
      AUC.perm[i] = mean(auc.perm)
    }
    #Calculate p-value
    #count the number of test-statistics as or more extreme than our initial test statistic
    Pval.acc = (length(which(ACCURACY.perm > average.acc)) + 1)/max.perm 
    Pval.spec = (length(which(SPECIFICITY.perm > average.spec)) + 1)/max.perm 
    Pval.sens = (length(which(SENSITIVITY.perm > average.sens)) + 1)/max.perm 
    Pval.auc = (length(which(AUC.perm > average.auc)) + 1)/max.perm 
    
    STAT[1,4] = Pval.acc
    STAT[2,4] = Pval.spec
    STAT[3,4] = Pval.sens
    STAT[4,4] = Pval.auc
    
    colnames(STAT) = c('Mean','95%CI Lower','95%CI Upper','P-val')
    rownames(STAT) = c('Accuracy','Specificity','Sensitivity','AUC')
    
    RF.models = list(ModelStatistics = STAT, RFaccuracy = ACCURACY, 
                     RFspecificity = SPECIFICITY, RFsensitivity = SENSITIVITY, 
                     RFaccuracyPerm = ACCURACY.perm, 
                     RFspecificityPerm = SPECIFICITY.perm, 
                     RFsensitivityPerm = SENSITIVITY.perm, RFaucPerm = AUC.perm)
  }
  
  ## Modelling if BALANCED
  if(isFALSE(unbalance)){
    STAT = matrix(NA,4,2)
    print('Correction and resampling for unblance NOT performed. Fitting a regular RF model')
    rf <- randomForest(as.factor(y.class) ~ ., data=x.data)
    print('done')
    #Extract predicted y
    y.pred = as.factor(rf$predicted)
    #  xtab <- table(as.numeric(y.pred), as.numeric(y.class))
    conf.mat = confusionMatrix(as.factor(y.pred), as.factor(y.class))
    
    model.acc = conf.mat$overall["Accuracy"]
    model.sens = conf.mat$byClass["Sensitivity"]
    model.spec = conf.mat$byClass["Specificity"]
    
    #Calculate AUC
    model.auc = colAUC(rf$votes[,2],y.class,plotROC=F)
    
    for(i in 1 : max.perm ){
      
      if(isTRUE(verbose))
        print(sprintf('Permutation %i of %i',i,max.perm ))
      
      #Pemute class labels
      idx.p = sample(n.samples)
      y.class.perm = y.class[idx.p]
      rf.perm <- randomForest(as.factor(y.class.perm) ~ ., data= x.data)
      y.pred.perm = as.factor(rf.perm$predicted)
      #    xtab.perm <- table(as.numeric(y.pred.perm), as.numeric(y.class.perm))
      # conf.mat.perm = confusionMatrix(xtab.perm)
      conf.mat.perm = confusionMatrix(as.factor(y.pred.perm), 
                                      as.factor(y.class.perm))
      acc.perm = conf.mat.perm$overall["Accuracy"]
      sens.perm = conf.mat.perm$byClass["Sensitivity"]
      spec.perm = conf.mat.perm$byClass["Specificity"]
      auc.perm = colAUC(rf.perm$votes[,2],y.class.perm,plotROC=F)
      
      ACCURACY.perm[i] = acc.perm
      SENSITIVITY.perm[i] = sens.perm
      SPECIFICITY.perm[i] = spec.perm
      AUC.perm[i] = auc.perm
    }
    
    # print(AUC.perm)
    #Calculate p-value
    
    # ACCURACY.perm =  read.table("C:\\Users\\Francesca\\Dropbox\\Forestone\\P_vs_C\\RFaccuracyPerm_PC.csv", sep = ",")
    # model.acc = read.table("C:\\Users\\Francesca\\Dropbox\\Forestone\\P_vs_C\\RFaccuracy_PC.csv", sep = ",")
    
    Pval.acc = (length(which(ACCURACY.perm > model.acc)) + 1)/max.perm 
    Pval.spec = (length(which(SPECIFICITY.perm > model.spec)) + 1)/max.perm 
    Pval.sens = (length(which(SENSITIVITY.perm > model.sens)) + 1)/max.perm 
    Pval.auc = (length(which(AUC.perm > model.auc[1])) + 1)/max.perm 
    
    STAT[1,1] = model.acc
    STAT[2,1] = model.spec
    STAT[3,1] = model.sens
    STAT[4,1] = model.auc
    
    STAT[1,2] = Pval.acc
    STAT[2,2] = Pval.spec
    STAT[3,2] = Pval.sens
    STAT[4,2] = Pval.auc
    
    colnames(STAT) = c('Mean','P-val')
    rownames(STAT) = c('Accuracy','Specificity','Sensitivity','AUC')
    
    RF.models = list(ModelStatistics = STAT, RFaccuracy = model.acc, 
                     RFspecificity = model.spec, RFsensitivity = model.sens, 
                     RFaccuracyPerm = ACCURACY.perm, 
                     RFspecificityPerm = SPECIFICITY.perm, 
                     RFsensitivityPerm = SENSITIVITY.perm, RFaucPerm = AUC.perm)
  }
  return(RF.models)
}
require(caret)        # confusionMatrix
require(caTools)      # colAUC

###############################################################################

women <- data[which(data$Gender == 'woman'),]
men <- data[which(data$Gender == 'man'),]

summary.table <- data.frame('Model', 'Threshold Young', 'Threshold Old', 
                              'Group size Young', 'Group size Old', 'Accuracy')
colnames(summary.table) <- summary.table[1,]

for (young in 32:38) {
  for (old in 43:49) {
    print(i)
    # Make age groups
    age.groups <- rep(0, nrow(data))
    age.groups[which(data$Age < young)] <- 'young'
    age.groups[which(data$Age > old)] <- 'old'
    age.women <- age.groups[which(data$Gender == 'woman')]
    age.men <- age.groups[which(data$Gender == 'man')]
    
    # All
    all.model <- RForest(x.data = data[which(data$Age < young | data$Age > old), 23:43], 
                   y.class = age.groups[which(age.groups == 'young' | age.groups == 'old')], 
                   unbalance = T, verbose = F, max.perm = 1)
    summary.table <- rbind(summary.table, c('All', young, old, length(age.groups[which(age.groups == 'young')]),
                           length(age.groups[which(age.groups == 'old')]), round(all.model$ModelStatistics[1,1],3)))
    # print('All data')
    # print(paste('Young =', young, 'Old =', old, 'Accuracy =', round(all.model$ModelStatistics[1,1],3)))
    # print('Group sizes')
    # print(paste('Young =', length(age.groups[which(age.groups == 'young')]), 'Old =', 
    #             length(age.groups[which(age.groups == 'old')])))
    
    # Women
    women.model <- RForest(x.data = women[which(women$Age < young | women$Age > old), 23:43], 
                   y.class = age.women[which(age.women == 'young' | age.women == 'old')], 
                   unbalance = T, verbose = F, max.perm = 1)
    summary.table <- rbind(summary.table, c('Women', young, old, length(age.women[which(age.women == 'young')]),
                           length(age.women[which(age.women == 'old')]), round(women.model$ModelStatistics[1,1],3)))
    # print('Women')
    # print(paste('Young =', young, 'Old =', old, 'Accuracy =', round(women.model$ModelStatistics[1,1],3)))
    # print('Group sizes')
    # print(paste('Young =', length(age.women[which(age.women == 'young')]), 'Old =', 
    #             length(age.women[which(age.women == 'old')])))
    
    # Men
    men.model <- RForest(x.data = men[which(men$Age < young | men$Age > old), 23:43], 
                   y.class = age.men[which(age.men == 'young' | age.men == 'old')], 
                   unbalance = T, verbose = F, max.perm = 1)
    summary.table <- rbind(summary.table, c('Men', young, old, length(age.men[which(age.men == 'young')]),
                           length(age.men[which(age.men == 'old')]), round(men.model$ModelStatistics[1,1],3)))
    # print('Men')
    # print(paste('Young =', young, 'Old =', old, 'Accuracy =', round(men.model$ModelStatistics[1,1],3)))
    # print('Group sizes')
    # print(paste('Young =', length(age.men[which(age.men == 'young')]), 'Old =', 
    #             length(age.men[which(age.men == 'old')])))
  }
}

###############################################################################

# Make 3D graph
library("plot3D")
x <- as.numeric(summary.table$`Threshold Young`[2:148])
y <- as.numeric(summary.table$`Threshold Old`[2:148])
z <- as.numeric(summary.table$Accuracy[2:148])

scatter3D(x, y, z, bty = "f", colvar = NULL, col = 'blue', main = "Accuracy for different age thresholds",
          pch = 19, cex = 0.5, xlab = "Threshold Young",
          ylab ="Threshold Old", zlab = "Accuracy", ticktype = "detailed")
