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
source('../AllNetworkFunctions12-Jan-2021 10.56.R')
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
