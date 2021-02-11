###############################################################################

# Yasmijn Balder
# 11-02-2021

# Univariate data analysis between age and sex

# Output
# Lists of significantly different lipids between:
    # Men and Women
    # Young and Old
    # Young Women and Old women
    # Young Men and Old men
# These lists are also written to csv files

###############################################################################

# Load packages

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/Lipids_age_sex.csv", check.names = FALSE)
data[,1] <- NULL

# Set different working directory to save results
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results")

###############################################################################

# Create an empty dataframe to store the results
wilcoxon.summary <- rbind(c('Lipid','Sex','','Age','', 'Men', '', 'Women', ''),
                          c('', 'Men', 'Women', 'Young', 'Old', 'Young men',
                            'Old men', 'Young women', 'Old women'))
for (lipid in 4:ncol(data)) {
  wilcoxon.summary <- rbind(wilcoxon.summary, c(colnames(data)[lipid], rep('', 8)))
}

###############################################################################

## Sex

# Splits the data into men and women
men <- data[which(data$Gender=='man'),]
women <- data[which(data$Gender=='woman'),]

# Initiate a vector to store all p.values
P.values <- NULL

# Perform a Wilcoxon rank-sum test for every lipid comparing men and women
for (lipid in 4:ncol(data)) {
  wil.test <- wilcox.test(men[,lipid], women[,lipid], alternative = "two.sided")
  P.values <- c(P.values, wil.test$p.value)
}

# Adjust the p.values with Bonferroni 
adjusted.P.values <- p.adjust(P.values, method = "bonferroni")

# # How many lipids are significant before and after?
# length(which(P.values <= 0.05))
# length(which(adjusted.P.values <= 0.05))
# # Plot the original and adjusted p.values to see the effect of the adjustment
# plot(adjusted.P.values, P.values)
# abline(h = 0.05, col = 'red'); abline(v = 0.05, col = 'red')

# Which lipids are significant different between men and women?
colnames(data)[which(adjusted.P.values <= 0.05) + 3]
# Write to csv file
write.csv(colnames(data)[which(adjusted.P.values <= 0.05) + 3], 
          file = 'Wilcoxon_Men_Women.csv')

i = 3
for (result in adjusted.P.values) {
  if (result <= 0.05) {
    wilcoxon.summary[i, 2] <- '+'
    wilcoxon.summary[i, 3] <- '+'
  } else {
    wilcoxon.summary[i, 2] <- '.'
    wilcoxon.summary[i, 3] <- '.'
  }
  i = i + 1
}

###############################################################################

## Age

# Splits the data into young and old
young <- data[which(data$Age < quantile(data$Age, probs = 1/3)),]
old <- data[which(data$Age > quantile(data$Age, probs = 2/3)),]

# Initiate a vector to store all p.values
P.values <- NULL

# Perform a Wilcoxon rank-sum test for every lipid comparing young and old
for (lipid in 4:ncol(data)) {
  wil.test <- wilcox.test(young[,lipid], old[,lipid], alternative = "two.sided")
  P.values <- c(P.values, wil.test$p.value)
}

# Adjust the p.values with Bonferroni 
adjusted.P.values <- p.adjust(P.values, method = "bonferroni")

# Which lipids are significant different between young and old?
colnames(data)[which(adjusted.P.values <= 0.05) + 3]
# Write to csv file
write.csv(colnames(data)[which(adjusted.P.values <= 0.05) + 3], 
          file = 'Wilcoxon_Young_Old.csv')

i = 3
for (result in adjusted.P.values) {
  if (result <= 0.05) {
    wilcoxon.summary[i, 4] <- '+'
    wilcoxon.summary[i, 5] <- '+'
  } else {
    wilcoxon.summary[i, 4] <- '.'
    wilcoxon.summary[i, 5] <- '.'
  }
  i = i + 1
}

###############################################################################

## Age + men

# Splits the data into young men and old men
young.men <- men[men$Age < quantile(men$Age, probs = 1/3),]
old.men <- men[men$Age > quantile(men$Age, probs = 2/3),]

# Initiate a vector to store all p.values
P.values <- NULL

# Perform a Wilcoxon rank-sum test for every lipid comparing young and old
for (lipid in 4:ncol(data)) {
  wil.test <- wilcox.test(young.men[,lipid], old.men[,lipid], alternative = "two.sided")
  P.values <- c(P.values, wil.test$p.value)
}

# Adjust the p.values with Bonferroni 
adjusted.P.values <- p.adjust(P.values, method = "bonferroni")

# Which lipids are significant different between young men and old men?
colnames(data)[which(adjusted.P.values <= 0.05) + 3]
# Write to csv file
write.csv(colnames(data)[which(adjusted.P.values <= 0.05) + 3], 
          file = 'Wilcoxon_YoungMen_OldMen.csv')

i = 3
for (result in adjusted.P.values) {
  if (result <= 0.05) {
    wilcoxon.summary[i, 6] <- '+'
    wilcoxon.summary[i, 7] <- '+'
  } else {
    wilcoxon.summary[i, 6] <- '.'
    wilcoxon.summary[i, 7] <- '.'
  }
  i = i + 1
}

###############################################################################

## Age + women

# Splits the data into young women and old women
young.women <- women[women$Age < quantile(women$Age, probs = 1/3),]
old.women <- women[women$Age > quantile(women$Age, probs = 2/3),]

# Initiate a vector to store all p.values
P.values <- NULL

# Perform a Wilcoxon rank-sum test for every lipid comparing young and old
for (lipid in 4:ncol(data)) {
  wil.test <- wilcox.test(young.women[,lipid], old.women[,lipid], alternative = "two.sided")
  P.values <- c(P.values, wil.test$p.value)
  print(paste(lipid, lipid-3, wil.test$p.value))
}

P.values[65] <- 0

# Adjust the p.values with Bonferroni 
adjusted.P.values <- p.adjust(P.values, method = "bonferroni")

# Which lipids are significant different between young women and old women?
colnames(data)[which(adjusted.P.values <= 0.05) + 3]
# Write to csv file
write.csv(colnames(data)[which(adjusted.P.values <= 0.05) + 3], 
          file = 'Wilcoxon_YoungWomen_OldWomen.csv')

i = 3
for (result in adjusted.P.values) {
  if (result <= 0.05) {
    wilcoxon.summary[i, 8] <- '+'
    wilcoxon.summary[i, 9] <- '+'
  } else {
    wilcoxon.summary[i, 8] <- '.'
    wilcoxon.summary[i, 9] <- '.'
  }
  i = i + 1
}

