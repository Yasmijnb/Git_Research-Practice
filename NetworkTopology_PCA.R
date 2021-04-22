###############################################################################

# Yasmijn Balder
# 10-03-2021

# PCA of network topology

# Output
# Two plot types for each PCA:
  # A plot that colours the two networks
  # A plot that colours the lipid groups and shapes the networks

###############################################################################

# Load packages
library(ggfortify)      # Used to plot the PCA
library(RNOmni)         # Used to add confidence interval ellipses

###############################################################################

# Load the data
men <- read.csv('Results/NetworkTopology_Men.csv')
women <- read.csv('Results/NetworkTopology_Women.csv')

young.men <- read.csv('Results/NetworkTopology_YoungMen.csv')
old.men <- read.csv('Results/NetworkTopology_OldMen.csv')

young.women <- read.csv('Results/NetworkTopology_YoungWomen.csv')
old.women <- read.csv('Results/NetworkTopology_OldWomen.csv')

###############################################################################

# Make a PCA for sex
SEX <- rbind(men, women)
rownames(SEX) <- c(1:21, '1 ', '2 ', '3 ', '4 ', '5 ', '6 ', '7 ', '8 ', '9 ', 
                   '10 ', '11 ', '12 ', '13 ', '14 ', '15 ', '16 ', '17 ', 
                   '18 ', '19 ', '20 ', '21 ')
SEX$network <- c(rep('men', nrow(men)), rep('women', nrow(women)))
SEX.pca <- prcomp(SEX[,c(2:7, 12:13, 16, 20:21)], scale = T)

# Make a PCA for men
MEN <- rbind(young.men, old.men)
rownames(MEN) <- c(1:21, '1 ', '2 ', '3 ', '4 ', '5 ', '6 ', '7 ', '8 ', '9 ', 
                   '10 ', '11 ', '12 ', '13 ', '14 ', '15 ', '16 ', '17 ', 
                   '18 ', '19 ', '20 ', '21 ')
MEN$network <- c(rep('young men', nrow(young.men)), rep('old men', nrow(old.men)))
MEN.pca <- prcomp(MEN[,c(2:7, 12:13, 16, 20:21)], scale = T)

# Make a PCA for women
WOMEN <- rbind(young.women, old.women)
rownames(WOMEN) <- c(1:21, '1 ', '2 ', '3 ', '4 ', '5 ', '6 ', '7 ', '8 ', '9 ', 
                   '10 ', '11 ', '12 ', '13 ', '14 ', '15 ', '16 ', '17 ', 
                   '18 ', '19 ', '20 ', '21 ')
WOMEN$network <- c(rep('young women', nrow(young.women)), rep('old women', nrow(old.women)))
WOMEN.pca <- prcomp(WOMEN[,c(2:7, 12:13, 16, 20:21)], scale = T)

###############################################################################

# Create figures

# Make figures that colour network and number lipids
autoplot(SEX.pca, data = SEX, colour = 'network', label = T, shape = F) + 
  # Add a confidence interval ellipse
  stat_ellipse(level = 0.95, aes(group = SEX$network, color = SEX$network), type = "norm") +
  # Set the colours
  scale_colour_manual(values=c('blue', 'red')) +
  # Add circles for clarity
  geom_point(shape=1, size = 8, colour = rep(c("blue","red"), each=21), stroke=1.5) +
  # Change number of decimals of PVE
  xlab(paste0('PC1 (', round(summary(SEX.pca)$importance[2,1]*100, 1), '%)')) +
  ylab(paste0('PC1 (', round(summary(SEX.pca)$importance[2,2]*100, 1), '%)')) +
  # Make the plot with a white background
  theme_bw(base_size = 17) + 
  # Remove the legend
  theme(legend.position = 'none') + 
  # Make square plots
  theme(aspect.ratio = 1)

autoplot(MEN.pca, data = MEN, colour = 'network', label = T, shape = F) + 
  # Add a confidence interval ellipse
  stat_ellipse(level = 0.95, aes(group = MEN$network, color = MEN$network), type = "norm") +
  # Set the colours
  scale_colour_manual(values=c('darkblue', 'cyan')) +
  # Add circles for clarity
  geom_point(shape=1, size = 8, colour = rep(c("cyan","darkblue"), each=21), stroke=1.5) +
  # Change number of decimals of PVE
  xlab(paste0('PC1 (', round(summary(MEN.pca)$importance[2,1]*100, 1), '%)')) +
  ylab(paste0('PC1 (', round(summary(MEN.pca)$importance[2,2]*100, 1), '%)')) +
  # Make the plot with a white background
  theme_bw(base_size = 17) + 
  # Remove legend
  theme(legend.position = 'none') + 
  # Make square plots
  theme(aspect.ratio = 1)

autoplot(WOMEN.pca, data = WOMEN, colour = 'network', label = T, shape = F) + 
  # Add a confidence interval ellipse
  stat_ellipse(level = 0.95, aes(group = WOMEN$network, color = WOMEN$network), type = "norm") +
  # Set the colours
  scale_colour_manual(values=c('darkred', 'orange')) +
  # Add circles for clarity
  geom_point(shape=1, size = 8, colour = rep(c("orange","darkred"), each=21), stroke=1.5) +
  # Change number of decimals of PVE
  xlab(paste0('PC1 (', round(summary(WOMEN.pca)$importance[2,1]*100, 1), '%)')) +
  ylab(paste0('PC1 (', format(round(summary(WOMEN.pca)$importance[2,2]*100, 1), nsmall = 1), '%)')) +
  # Make the plot with a white background
  theme_bw(base_size = 17) + 
  # Remove legend
  theme(legend.position = 'none') + 
  # Make square plots
  theme(aspect.ratio = 1)

###############################################################################
