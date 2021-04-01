###############################################################################

# Yasmijn Balder
# 02-03-2021

# Use different thresholds for young and old to find which gives the most
# accurate devision.

# Output
# The accuracy, sentivitiy, specificity, and AUC of the model

###############################################################################

# Load data
data <- read.csv("../Data/LipidsAgeSex_SqrtNormalization.csv", check.names = TRUE)
data[,1] <- NULL

###############################################################################

# Packages
library(randomForest)
library(Rmisc)

# Functions
source('../AllNetworkFunctions12-Jan-2021 10.56.R')
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
