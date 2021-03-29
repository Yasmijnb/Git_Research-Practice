###############################################################################

# Yasmijn Balder
# 25-01-2021

# Investigate patient characteristics/demographics

# Output
# Plots showing number of individuals per age (one plot for men, one for women)

###############################################################################

# Load packages

library(ggplot2)

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/LipidsAgeSex_NoNormalization.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

# Create histograms of the age distribution per gender 
# (figure 1 of Vignoli et. al. 2017)

# Women
hist(data$Age[which(data$Gender == 'woman')], 
     breaks = length(unique(data$Age[which(data$Gender == 'woman')])), 
     xlab = 'Age (years)', ylab = 'Number of individuals', col = 'magenta',
     main = 'Age distribution of women')
abline(v = 37, col = 'blue')
abline(v = 48, col = 'blue')

ggplot(data[which(data$Gender == 'woman'),], aes(x = Age)) + 
  geom_histogram(binwidth = 1, color = 'white', fill = 'magenta') + 
  labs(title = "Age distribution of women",
       x = "Age (years)", y = "Number of individuals") + 
  scale_x_continuous(breaks = seq(range(data$Age)[1], range(data$Age)[2], 1)) + 
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Men
hist(data$Age[which(data$Gender == 'man')], 
     breaks = length(unique(data$Age[which(data$Gender == 'man')])), 
     xlab = 'Age (years)', ylab = 'Number of individuals', col = 'blue',
     main = 'Age distribution of men')
abline(v = 35, col = 'red')
abline(v = 45, col = 'red')

ggplot(data[which(data$Gender == 'man'),], aes(x = Age)) + 
  geom_histogram(binwidth = 1, color = 'white', fill = 'cyan') + 
  labs(title = "Age distribution of men",
       x = "Age (years)", y = "Number of individuals") + 
  scale_x_continuous(breaks = seq(range(data$Age)[1], range(data$Age)[2], 1)) + 
  scale_y_continuous(breaks = seq(0, 30, 5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

###############################################################################

# Find groups for young and old using quantiles

# Women
quantile(data$Age[which(data$Gender == 'woman')], probs = seq(0, 1, 1/3))
# Young < 37
# Old   > 48

# Men
quantile(data$Age[which(data$Gender == 'man')], probs = seq(0, 1, 1/3))
# Young < 35
# Old   > 45

###############################################################################

# Women total
women.number <- length(which(data$Gender == 'woman'))
women.age <- median(data[which(data$Gender == 'woman'),'Age'])

# Women young
y.women.number <- length(which(data$Gender == 'woman' & data$Age < 37))
y.women.age <- median(data[which(data$Gender == 'woman' & data$Age < 37),'Age'])

# Women old
o.women.number <- length(which(data$Gender == 'woman' & data$Age > 48))
o.women.age <- median(data[which(data$Gender == 'woman' & data$Age > 48),'Age'])

# Men total
men.number <- length(which(data$Gender == 'man'))
men.age <- median(data[which(data$Gender == 'man'),'Age'])

# Men young
y.men.number <- length(which(data$Gender == 'man' & data$Age < 35))
y.men.age <- median(data[which(data$Gender == 'man' & data$Age < 35),'Age'])

# Men old
o.men.number <- length(which(data$Gender == 'man' & data$Age > 45))
o.men.age <- median(data[which(data$Gender == 'man' & data$Age > 45),'Age'])

women.number; women.age
y.women.number; y.women.age
o.women.number; o.women.age
men.number; men.age
y.men.number; y.men.age
o.men.number; o.men.age



