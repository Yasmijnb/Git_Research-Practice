###############################################################################

# Yasmijn Balder
# 15-02-2021

# PLS-DA analysis between sex and age using main fractions

# Output
# 

###############################################################################

# Load packages
library(caret)

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/Lipids_age_sex.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

sex.plsda <- plsda(data[,23:43], as.factor(data$Gender))
