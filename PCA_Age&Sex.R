###############################################################################

# Yasmijn Balder
# 02-02-2021

# Investigate difference between age and sex

# Output
# PCA coloured by age and sex

###############################################################################

# Load packages

library(ggfortify)

###############################################################################

# Load data
data <- read.csv("../Data/LipidsAgeSex_SqrtNormalization.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

# Make groups of age and sex
group <- NULL
for (sample in 1:nrow(data)) {
  # Men
  if (data$Gender[sample] == 'man') {
    if (data$Age[sample] < 35) {
      group <- c(group, 'young man')
    } else 
      if (data$Age[sample] > 45) {
        group <- c(group, 'old man')
      } else
        group <- c(group, 'man')
  }
  # Women
  if (data$Gender[sample] == 'woman') {
    if (data$Age[sample] < 37) {
      group <- c(group, 'young woman')
    } else 
      if (data$Age[sample] > 48) {
        group <- c(group, 'old woman')
      } else
        group <- c(group, 'woman')
  }
}

data$Groups <- group

###############################################################################

# Perform PCA analysis
pca <- prcomp(data[,4:117])

# Percentage of variance explained
summary(pca)$importance[2:3,1:3]*100

# Loadings of first two components
pca$rotation[,1]^2
pca$rotation[,2]^2

# Create plots
autoplot(pca, data = data, colour = 'Gender') +
  # Make the plot with a white background
  theme_bw()
autoplot(pca, data = data, colour = 'Age') +
  # Make the plot with a white background
  theme_bw()
autoplot(pca, data = data, colour = 'Groups') +
  # Make the plot with a white background
  theme_bw()

###############################################################################

# Use only main fractions

# Perform PCA analysis
pca <- prcomp(data[,23:43])

# Percentage of variance explained
summary(pca)$importance[2:3,1:4]*100

# Loadings of first two components
pca$rotation[,1]^2
pca$rotation[,2]^2

# Create plots 
autoplot(pca, data = data, colour = 'Gender') +
  scale_colour_manual(values=c('blue', 'red')) +
  # Make the plot with a white background
  theme_bw()
# autoplot(pca, data = data, colour = 'Age') +
  # scale_colour_manual(values=c('darkblue', 'cyan')) +
  # Make the plot with a white background
  # theme_bw()
autoplot(pca, data = data, colour = 'Groups') +
  scale_colour_manual(values=c('blue', 'darkblue', 'darkred', 'red', 'cyan', 'pink')) +
  # Make the plot with a white background
  theme_bw()

###############################################################################

# Prepare data frames
men <- data[which(data$Gender == 'man'),]
men <- men[which(men$Age < 35 | men$Age > 45),]
women <- data[which(data$Gender == 'woman'),]
women <- women[which(women$Age < 37 | women$Age > 48),]

# Make three separate PCAs
sex.pca <- prcomp(data[,23:43])
men.pca <- prcomp(men[,23:43])
women.pca <- prcomp(women[,23:43])

# Create plots 
autoplot(sex.pca, data = data, colour = 'Gender') + 
  scale_colour_manual(values=c('blue', 'red')) +
  # Make the plot with a white background and make font size bigger
  theme_bw(base_size = 17) + 
  # Remove legend title
  theme(legend.title = element_blank())
  
autoplot(men.pca, data = men, colour = 'Groups') +
  scale_colour_manual(values=c('darkblue', 'cyan')) +
  # Make the plot with a white background and make font size bigger
  theme_bw(base_size = 17) +
  # Remove legend title
  theme(legend.title = element_blank())

autoplot(women.pca, data = women, colour = 'Groups') +
  scale_colour_manual(values=c('darkred', 'orange')) +
  # Make the plot with a white background and make font size bigger
  theme_bw(base_size = 17) +
  # Remove legend title
  theme(legend.title = element_blank())
