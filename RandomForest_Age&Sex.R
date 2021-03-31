###############################################################################

# Yasmijn Balder
# 12-02-2021

# Perform random forest to distinguish sex and age using main fractions. Four
# models are created: men vs. women, young vs. old, young men vs. old men, young
# women vs. old women

# Output
# The accuracy, sensitivity, specificity, and AUC of each model

###############################################################################

# Load data
data <- read.csv("../Data/LipidsAgeSex_SqrtNormalization.csv", check.names = T)
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

## Sex
sex.forest <- RForest(x.data = data[,23:43], y.class = as.factor(data$Gender))
sex.forest$ModelStatistics

###############################################################################

## Age

# Make young and old groups
# age_groups <- rep(0, nrow(data))
# age_groups[which(data$Age < quantile(data$Age, probs = 1/3))] <- 'young'
# age_groups[which(data$Age > quantile(data$Age, probs = 2/3))] <- 'old'
# 
# # Use only young and old, not the middle
# age.data <- data[which(age_groups == 'young' | age_groups == 'old'),]
# group.data <- age_groups[which(age_groups == 'young' | age_groups == 'old')]
# 
# age.forest <- RForest(x.data = age.data[,23:43], y.class = group.data, 
#                       unbalance = FALSE)
# age.forest$ModelStatistics

###############################################################################

## Women

# Use only women
women <- data[which(data$Gender == 'woman'),]

# Make young and old groups
women$AgeGroup <- rep(0, nrow(women))
women$AgeGroup[which(women$Age < quantile(women$Age, probs = 1/3))] <- 'young'
women$AgeGroup[which(women$Age > quantile(women$Age, probs = 2/3))] <- 'old'

# Use only young and old, not the middle
women <- women[which(women$AgeGroup == 'young' | women$AgeGroup == 'old'),]
# women.group.data <- women.age_groups[which(women.age_groups == 'young' | women.age_groups == 'old')]

women.forest <- RForest(x.data = women[,23:43], y.class = women$AgeGroup, 
                        unbalance = FALSE)
women.forest$ModelStatistics

###############################################################################

## Men

# Use only men
men <- data[which(data$Gender == 'man'),]

# Make young and old groups
men$AgeGroup <- rep(0, nrow(men))
men$AgeGroup[which(men$Age < quantile(men$Age, probs = 1/3))] <- 'young'
men$AgeGroup[which(men$Age > quantile(men$Age, probs = 2/3))] <- 'old'

# Use only young and old, not the middle
men <- men[which(men$AgeGroup == 'young' | men$AgeGroup == 'old'),]
# men.group.data <- men.age_groups[which(men.age_groups == 'young' | men.age_groups == 'old')]

men.forest <- RForest(x.data = men[,23:43], y.class = men$AgeGroup, 
                      unbalance = FALSE)
men.forest$ModelStatistics

###############################################################################

# Print all results again
sex.forest$ModelStatistics
# age.forest$ModelStatistics
women.forest$ModelStatistics
men.forest$ModelStatistics

###############################################################################

## Run a single random forest to get an idea of the importance

# Shorten colnames
library(stringr)            # Used to split and shorten the lipid names
# Shorten the names of the lipids
data <- read.csv("../Data/LipidsAgeSex_SqrtNormalization.csv", check.names = F)
data[,1] <- NULL
short.names <- NULL
long.names <- colnames(data)[23:43]
for (index in long.names) {
  short.name <- str_split(index, ', ')[[1]][2:3]
  combined <- paste(short.name, collapse = ' & ')
  short.names <- c(short.names, combined)
}
colnames(data)[23:43] <- short.names

# Take care of the unbalance in gender
n1 = sum(as.factor(data$Gender) == 'man')
n2 = sum(as.factor(data$Gender) == 'woman')
nsize = round(min(n1,n2)*0.85)
# Prepare for women
women <- data[which(data$Gender == 'woman'),]
women$AgeGroup <- rep(0, nrow(women))
women$AgeGroup[which(women$Age < quantile(women$Age, probs = 1/3))] <- 'young'
women$AgeGroup[which(women$Age > quantile(women$Age, probs = 2/3))] <- 'old'
women <- women[which(women$AgeGroup == 'young' | women$AgeGroup == 'old'),]
# Prepare for men
men <- data[which(data$Gender == 'man'),]
men$AgeGroup <- rep(0, nrow(men))
men$AgeGroup[which(men$Age < quantile(men$Age, probs = 1/3))] <- 'young'
men$AgeGroup[which(men$Age > quantile(men$Age, probs = 2/3))] <- 'old'
men <- men[which(men$AgeGroup == 'young' | men$AgeGroup == 'old'),]

# Run single random forests
single.sex.forest <- randomForest(data[,23:43], as.factor(data$Gender), 
                                  importance = T, strata = as.factor(data$Gender), 
                                  sampsize = c(nsize,nsize), replace = F,
                                  proximity = T)
single.women.forest <- randomForest(women[,23:43], as.factor(women$AgeGroup), 
                                    importance = T, proximity = T)
single.men.forest <- randomForest(men[,23:43], as.factor(men$AgeGroup), 
                                  importance = T, proximity = T)

# View importance
varImpPlot(single.sex.forest)
varImpPlot(single.women.forest)
varImpPlot(single.men.forest)

# View importance
varImpPlot(single.sex.forest, sort = F)
varImpPlot(single.women.forest, sort = F)
varImpPlot(single.men.forest, sort = F)

###############################################################################

# Proximity distance plots

MDSplot(single.sex.forest, as.factor(data$Gender), palette = c('blue', 'red'))
MDSplot(single.women.forest, as.factor(women$AgeGroup), palette = c('darkred','pink'))
MDSplot(single.men.forest, as.factor(men$AgeGroup), palette = c('darkblue', 'cyan'))
