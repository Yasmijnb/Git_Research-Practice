###############################################################################

# Yasmijn Balder
# 20-01-2021

# Check the correlation between the measured and estimated lipids

# Output
  # Plots between estimated and measured values with the correlation value in
  # the title of the plot

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/Lipids_all.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

# Check the correlation between certain estimations
# Ignore NAs

Chol.cor <- cor(data$Chol, data$`Main Parameters, Cholesterol, Chol`, use = 'complete.obs')

plot(data$Chol, data$`Main Parameters, Cholesterol, Chol`,
     main = paste0('Correlation between measured and estimated Cholesterol \n(',
                   round(Chol.cor, 2), ')'),
     xlab = 'Estimated', ylab = 'Measured')
# abline(0, Chol.cor, col = 'red')

########################################

# TG
TG.cor <- cor(data$TG, data$`Main Parameters, Triglycerides, TG`, use = 'complete.obs')

plot(data$TG, data$`Main Parameters, Triglycerides, TG`,
     main = paste0('Correlation between measured and estimated TG \n(',
                   round(TG.cor, 2), ')'),
     xlab = 'Estimated', ylab = 'Measured')
# abline(0, 1, col = 'red')

########################################

# HDL-Chol	
HDL.cor <- cor(data$`HDL-Chol`, 
               data$`Main Parameters, HDL Cholesterol, HDL-Chol`, use = 'complete.obs')

plot(data$`HDL-Chol`, 
     data$`Main Parameters, HDL Cholesterol, HDL-Chol`,
     main = paste0('Correlation between measured and estimated HDL-Chol \n(',
                   round(HDL.cor, 2), ')'),
     xlab = 'Estimated', ylab = 'Measured')
# abline(0, 1, col = 'red')

########################################

# LDL-Chol
LDL.cor <- cor(data$`LDL-Chol`, 
               data$`Main Parameters, LDL Cholesterol, LDL-Chol`, use = 'complete.obs')

plot(data$`LDL-Chol`, 
     data$`Main Parameters, LDL Cholesterol, LDL-Chol`,
     main = paste0('Correlation between measured and estimated LDL-Chol \n(',
                   round(LDL.cor, 2), ')'),
     xlab = 'Estimated', ylab = 'Measured')
# abline(0, 1, col = 'red')

###############################################################################

# If all NAs are 0

data[which(is.na(data), arr.ind = T)] <- 0

# Chol
Chol.cor <- cor(data$Chol, data$`Main Parameters, Cholesterol, Chol`)

plot(data$Chol, data$`Main Parameters, Cholesterol, Chol`,
     main = paste0('Correlation between measured and estimated Cholesterol \n(',
                   round(Chol.cor, 2), ')'),
     xlab = 'Estimated', ylab = 'Measured')

########################################

# TG
TG.cor <- cor(data$TG, data$`Main Parameters, Triglycerides, TG`)

plot(data$TG, data$`Main Parameters, Triglycerides, TG`,
     main = paste0('Correlation between measured and estimated TG \n(',
                   round(TG.cor, 2), ')'),
     xlab = 'Estimated', ylab = 'Measured')

########################################

# HDL-Chol	
HDL.cor <- cor(data$`HDL-Chol`, 
               data$`Main Parameters, HDL Cholesterol, HDL-Chol`)

plot(data$`HDL-Chol`, 
     data$`Main Parameters, HDL Cholesterol, HDL-Chol`,
     main = paste0('Correlation between measured and estimated HDL-Chol \n(',
                   round(HDL.cor, 2), ')'),
     xlab = 'Estimated', ylab = 'Measured')

########################################

# LDL-Chol
LDL.cor <- cor(data$`LDL-Chol`, 
               data$`Main Parameters, LDL Cholesterol, LDL-Chol`)

plot(data$`LDL-Chol`, 
     data$`Main Parameters, LDL Cholesterol, LDL-Chol`,
     main = paste0('Correlation between measured and estimated LDL-Chol \n(',
                   round(LDL.cor, 2), ')'),
     xlab = 'Estimated', ylab = 'Measured')
