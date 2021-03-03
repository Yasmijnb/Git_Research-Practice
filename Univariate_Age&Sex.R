###############################################################################

# Yasmijn Balder
# 11-02-2021

# Univariate data analysis using only main fractions
    # Men and Women
    # Young and Old
    # Young Women and Old women
    # Young Men and Old men

# Output
# Table
# This is also written to a csv file

###############################################################################

# Load packages

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/Lipids_age_sex.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

# Initiate a vector to store all p.values
P.values <- NULL
# Initiate a vector to store the shifts
shifts <- NULL

###############################################################################

## Sex

# Splits the data into men and women
men <- data[which(data$Gender=='man'),]
women <- data[which(data$Gender=='woman'),]

# Perform a Wilcoxon rank-sum test for every lipid comparing men and women
# for (lipid in 4:ncol(data)) { # All lipid measurements
for (lipid in 23:43) {
  wil.test <- wilcox.test(men[,lipid], women[,lipid], 
                          alternative = "two.sided", conf.int = T)
  P.values <- c(P.values, wil.test$p.value)
  shifts <- c(shifts, wil.test$estimate)
}

###############################################################################

## Age

# Splits the data into young and old
young <- data[which(data$Age < 33),]
old <- data[which(data$Age > 48),]

# Perform a Wilcoxon rank-sum test for every lipid comparing young and old
for (lipid in 23:43) {
# for (lipid in 4:ncol(data)) { # All lipid measurements
  wil.test <- wilcox.test(young[,lipid], old[,lipid], 
                          alternative = "two.sided", conf.int = T)
  P.values <- c(P.values, wil.test$p.value)
  shifts <- c(shifts, wil.test$estimate)
}

###############################################################################

## Age + men

# Splits the data into young men and old men
young.men <- men[men$Age < 32,]
old.men <- men[men$Age > 49,]

# Perform a Wilcoxon rank-sum test for every lipid comparing young and old
# for (lipid in 4:ncol(data)) { # All lipid measurements
for (lipid in 23:43) {
  wil.test <- wilcox.test(young.men[,lipid], old.men[,lipid], 
                          alternative = "two.sided", conf.int = T)
  P.values <- c(P.values, wil.test$p.value)
  shifts <- c(shifts, wil.test$estimate)
}

###############################################################################

## Age + women

# Splits the data into young women and old women
young.women <- women[women$Age < 37,]
old.women <- women[women$Age > 49,]

# Perform a Wilcoxon rank-sum test for every lipid comparing young and old
# for (lipid in 4:ncol(data)) { # All lipid measurements
for (lipid in 23:43) {
  # One sample contains only zero's which the test cannot handle
  if (sum(young.women[,lipid]) == 0) {
    P.values <- c(P.values, 1)
    shifts <- c(shifts, 0)
  } else{
  wil.test <- wilcox.test(young.women[,lipid], old.women[,lipid], 
                          alternative = "two.sided", conf.int = T)
  P.values <- c(P.values, wil.test$p.value)
  shifts <- c(shifts, wil.test$estimate)
  }
}

###############################################################################

# Compile the results

# Adjust the p-values
adjusted.P.values <- p.adjust(P.values, method = 'bonferroni')

# Make a vector with the directions of the shifts
directions <- rep('', length(adjusted.P.values))
# + means a higher value for women (Sex), old (Age, Women, Men)
directions[which(adjusted.P.values <= 0.01 & shifts < 0)] <- '+' # Reverse is intentional!!!
directions[which(adjusted.P.values <= 0.01 & shifts > 0)] <- '-' # Reverse is intentional!!!

# Make an empty dataframe to store the results
# wilcoxon.summary <- matrix(ncol = 4, nrow = 114) # All lipid measurements
wilcoxon.summary <- matrix(ncol = 4, nrow = 21)
colnames(wilcoxon.summary) <- c('Sex', 'Age', 'Women', 'Men')
# rownames(wilcoxon.summary) <- colnames(data)[4:117] # All lipid measurements
rownames(wilcoxon.summary) <- colnames(data)[23:43]

# Fill in the directions in the summary table
wilcoxon.summary[1:21,1] <- directions[1:21]
wilcoxon.summary[1:21,2] <- directions[22:42]
wilcoxon.summary[1:21,3] <- directions[43:63]
wilcoxon.summary[1:21,4] <- directions[64:84]

# How many lipids are significant before and after?
length(which(P.values <= 0.01))
length(which(adjusted.P.values <= 0.01))
# Plot the original and adjusted p.values to see the effect of the adjustment
plot(P.values, adjusted.P.values)
abline(h = 0.01, col = 'red'); abline(v = 0.01, col = 'red')

# Set different working directory results
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results")
# Write the summary to csv file
write.csv(wilcoxon.summary, file = 'Wilcoxon_results.csv')
