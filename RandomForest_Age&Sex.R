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

# Print all results again
sex.forest$ModelStatistics
men.forest$ModelStatistics
women.forest$ModelStatistics

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
# Prepare for men
men <- data[which(data$Gender == 'man'),]
men$AgeGroup <- rep(0, nrow(men))
men$AgeGroup[which(men$Age < quantile(men$Age, probs = 1/3))] <- 'young'
men$AgeGroup[which(men$Age > quantile(men$Age, probs = 2/3))] <- 'old'
men <- men[which(men$AgeGroup == 'young' | men$AgeGroup == 'old'),]
# Prepare for women
women <- data[which(data$Gender == 'woman'),]
women$AgeGroup <- rep(0, nrow(women))
women$AgeGroup[which(women$Age < quantile(women$Age, probs = 1/3))] <- 'young'
women$AgeGroup[which(women$Age > quantile(women$Age, probs = 2/3))] <- 'old'
women <- women[which(women$AgeGroup == 'young' | women$AgeGroup == 'old'),]

# # Run single random forests
# single.sex.forest <- randomForest(data[,23:43], as.factor(data$Gender), 
#                                   importance = T, strata = as.factor(data$Gender), 
#                                   sampsize = c(nsize,nsize), replace = F,
#                                   proximity = T)
# single.men.forest <- randomForest(men[,23:43], as.factor(men$AgeGroup), 
#                                   importance = T, proximity = T)
# single.women.forest <- randomForest(women[,23:43], as.factor(women$AgeGroup), 
#                                     importance = T, proximity = T)
# 
# # View importance
# varImpPlot(single.sex.forest)
# varImpPlot(single.men.forest)
# varImpPlot(single.women.forest)
# 
# # View importance
# varImpPlot(single.sex.forest, sort = F)
# varImpPlot(single.men.forest, sort = F)
# varImpPlot(single.women.forest, sort = F)

###############################################################################

# Proximity distance plots
# 
# MDSplot(single.sex.forest, as.factor(data$Gender), palette = c('blue', 'red'))
# MDSplot(single.men.forest, as.factor(men$AgeGroup), palette = c('darkblue', 'cyan'))
# MDSplot(single.women.forest, as.factor(women$AgeGroup), palette = c('darkred','pink'))

###############################################################################

## Run rfPermuta to get the importance
library(rfPermute)

# Run random forests with permutations
permute.sex.forest <- rfPermute(
  # Use same random forest arguments
  data[,23:43], as.factor(data$Gender), importance = T, proximity = T,
  strata = as.factor(data$Gender), sampsize = c(nsize,nsize), replace = F,
  # Use 1000 permutations 
  nrep = 100,
  # Use multiple cores
  num.cores = 5)
permute.men.forest <- rfPermute(men[,23:43], as.factor(men$AgeGroup), 
                                importance = T, proximity = T, nrep = 100, 
                                num.cores = 5)
permute.women.forest <- rfPermute(women[,23:43], as.factor(women$AgeGroup), 
                                  importance = T, proximity = T, nrep = 100,
                                  num.cores = 5)

# Prepare for gg plot
imp.sex.forest <- as.data.frame(rp.importance(permute.sex.forest))
imp.sex.forest$Lipids <- rownames(imp.sex.forest)
imp.sex.forest$sig <- rep("Not significant", 21)
imp.sex.forest$sig[which(imp.sex.forest$MeanDecreaseGini.pval < 0.05)] <- 'Significant'

imp.men.forest <- as.data.frame(rp.importance(permute.men.forest))
imp.men.forest$Lipids <- rownames(imp.men.forest)
imp.men.forest$sig <- rep("Not significant", 21)
imp.men.forest$sig[which(imp.men.forest$MeanDecreaseGini.pval < 0.05)] <- 'Significant'

imp.women.forest <- as.data.frame(rp.importance(permute.women.forest))
imp.women.forest$Lipids <- rownames(imp.women.forest)
imp.women.forest$sig <- rep("Not significant", 21)
imp.women.forest$sig[which(imp.women.forest$MeanDecreaseGini.pval < 0.05)] <- 'Significant'

# Plot importance
ggplot(data = imp.sex.forest, aes(x = MeanDecreaseGini, y = Lipids, 
                                  fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("steelblue", "orange")) + 
  # Change the axis labels
  xlab('Mean Decrease Gini') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank())

ggplot(data = imp.men.forest, aes(x = MeanDecreaseGini, y = Lipids, 
                                  fill = sig)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("steelblue", "orange")) + 
  xlab('Mean Decrease Gini') + ylab('Lipoprotein Main Fractions') + 
  theme_bw(base_size = 17) +
  theme(legend.title = element_blank())

ggplot(data = imp.women.forest, aes(x = MeanDecreaseGini, y = Lipids, 
                                  fill = sig)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("steelblue", "orange")) + 
  xlab('Mean Decrease Gini') + ylab('Lipoprotein Main Fractions') +
  theme_bw(base_size = 17) +
  theme(legend.title = element_blank()) +
  xlim(0, 16)

