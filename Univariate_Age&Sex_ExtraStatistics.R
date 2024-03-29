###############################################################################

# Yasmijn Balder
# 01-04-2021

# More stats are compiled from the three comparisons per lipid. 
# Comparisons:
  # Men vs. Women
  # Young vs. old men
  # Young vs. old women
# Stats:
  # AUC
  # CI.AUC
  # threshold
  # accuracy
  # specificity
  # sensitivity
  # p-value
  # Adjusted p-value

# Output
# Three tables are saved as .csv 

###############################################################################

# Load data
data <- read.csv("../Data/LipidsAgeSex_NoNormalization.csv", check.names = F)
data[,1] <- NULL

###############################################################################

# Packages
library(pROC)

###############################################################################

# Get short names
short.names <- NULL
long.names <- colnames(data)[23:43]
for (index in long.names) {
  short.name <- str_split(index, ', ')[[1]][2:3]
  combined <- paste(short.name, collapse = ' & ')
  short.names <- c(short.names, combined)
}

###############################################################################

# Make the gender column a factor 
data$Gender <- as.factor(data$Gender)

# Splits the data into men and women
men <- data[which(data$Gender=='man'),]
women <- data[which(data$Gender=='woman'),]
# Splits the data into young men and old men
young.men <- men[men$Age < quantile(men$Age, probs = 1/3),]
old.men <- men[men$Age > quantile(men$Age, probs = 2/3),]
# Splits the data into young women and old women
young.women <- women[women$Age < quantile(women$Age, probs = 1/3),]
old.women <- women[women$Age > quantile(women$Age, probs = 2/3),]

# Make age groups
men.age <- NULL
for (sample in 1:nrow(men)) {
  if (men$Age[sample] < 35) {
    men.age <- c(men.age, 'young')
  } 
  if (men$Age[sample] > 45) {
    men.age <- c(men.age, 'old')
  } 
}
men.age <- as.factor(men.age)
women.age <- NULL
for (sample in 1:nrow(women)) {
  if (women$Age[sample] < 37) {
    women.age <- c(women.age, 'young')
  } 
    if (women$Age[sample] > 48) {
      women.age <- c(women.age, 'old')
    } 
}
women.age <- as.factor(women.age)

# Remove middle group
youngold.men <- men[which(men$Age < 35 | men$Age > 45),]
youngold.women <- women[which(women$Age < 37 | women$Age > 48),]

###############################################################################

## Sex 

# Create an empty dataframe to save the stats
sex.stats <- NULL

# Perform a Wilcoxon rank-sum test for every lipid comparing men and women
for (lipid in 23:43) {
  RES = roc(data$Gender, data[,lipid], levels = levels(data$Gender))
  wil.test <- wilcox.test(men[,lipid], women[,lipid], alternative = "two.sided")
  pval <- wil.test$p.value
  auc = RES$auc
  CO = coords(RES, x="best", input="threshold", ret=c("threshold","accuracy","specificity", "sensitivity"))
  threshold= CO$threshold
  accuracy = CO$accuracy
  specificity = CO$specificity
  sensitivity = CO$sensitivity
  CI.AUC = ci.auc(RES, x="best")
  # Collect the stats
  stats <- c(auc, CI.AUC[1], CI.AUC[3], threshold, accuracy, specificity, sensitivity, pval)
  sex.stats <- rbind(sex.stats, stats)
}

sex.stats <- as.data.frame(sex.stats)
colnames(sex.stats) <- c('AUC', 'CI.AUC Lower', 'CI.AUC Upper', 'Threshold', 
                         'Accuracy', 'Specificity', 'Sensitivity', 'P-value')
rownames(sex.stats) <- short.names

###############################################################################

## Men 

# Create an empty dataframe to save the stats
men.stats <- NULL

# Perform a Wilcoxon rank-sum test for every lipid comparing men and women
for (lipid in 23:43) {
  RES = roc(men.age, youngold.men[,lipid], levels = levels(men.age))
  wil.test <- wilcox.test(young.men[,lipid], old.men[,lipid], alternative = "two.sided")
  pval <- wil.test$p.value
  auc = RES$auc
  CO = coords(RES, x="best", input="threshold", ret=c("threshold","accuracy","specificity", "sensitivity"))
  threshold= CO$threshold
  accuracy = CO$accuracy
  specificity = CO$specificity
  sensitivity = CO$sensitivity
  CI.AUC = ci.auc(RES, x="best")
  # Collect the stats
  stats <- c(auc, CI.AUC[1], CI.AUC[3], threshold, accuracy, specificity, sensitivity, pval)
  men.stats <- rbind(men.stats, stats)
}

men.stats <- as.data.frame(men.stats)
colnames(men.stats) <- c('AUC', 'CI.AUC Lower', 'CI.AUC Upper', 'Threshold', 
                         'Accuracy', 'Specificity', 'Sensitivity', 'P-value')
rownames(men.stats) <- short.names

###############################################################################

## Women 

# Create an empty dataframe to save the stats
women.stats <- NULL

# Perform a Wilcoxon rank-sum test for every lipid comparing men and women
for (lipid in 23:43) {
  RES = roc(women.age, youngold.women[,lipid], levels = levels(women.age))
  wil.test <- wilcox.test(young.women[,lipid], old.women[,lipid], alternative = "two.sided")
  pval <- wil.test$p.value
  auc = RES$auc
  CO = coords(RES, x="best", input="threshold", ret=c("threshold","accuracy","specificity", "sensitivity"))
  threshold= CO$threshold
  accuracy = CO$accuracy
  specificity = CO$specificity
  sensitivity = CO$sensitivity
  CI.AUC = ci.auc(RES, x="best")
  # Collect the stats
  stats <- c(auc, CI.AUC[1], CI.AUC[3], threshold, accuracy, specificity, sensitivity, pval)
  women.stats <- rbind(women.stats, stats)
}

women.stats <- as.data.frame(women.stats)
colnames(women.stats) <- c('AUC', 'CI.AUC Lower', 'CI.AUC Upper', 'Threshold', 
                         'Accuracy', 'Specificity', 'Sensitivity', 'P-value')
rownames(women.stats) <- short.names

###############################################################################

# Adjust the p-values
pvalues <- c(sex.stats$`P-value`, men.stats$`P-value`, women.stats$`P-value`)
adjusted <- p.adjust(pvalues, method = 'bonferroni')
sex.stats$`Adjusted P-value` <- adjusted[1:21]
men.stats$`Adjusted P-value` <- adjusted[22:42]
women.stats$`Adjusted P-value` <- adjusted[43:63]

###############################################################################

sex.stats <- round(sex.stats, 3)
men.stats <- round(men.stats, 3)
women.stats <- round(women.stats, 3)

write.csv(sex.stats, 'Results/Wilcoxon_stats_Sex.csv')
write.csv(men.stats, 'Results/Wilcoxon_stats_Men.csv')
write.csv(women.stats, 'Results/Wilcoxon_stats_Women.csv')
