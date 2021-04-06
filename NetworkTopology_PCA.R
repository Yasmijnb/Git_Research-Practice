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

young <- read.csv('Results/NetworkTopology_Young.csv')
old <- read.csv('Results/NetworkTopology_Old.csv')

young.men <- read.csv('Results/NetworkTopology_YoungMen.csv')
old.men <- read.csv('Results/NetworkTopology_OldMen.csv')

young.women <- read.csv('Results/NetworkTopology_YoungWomen.csv')
old.women <- read.csv('Results/NetworkTopology_OldWomen.csv')

###############################################################################

# Make a PCA for age
AGE <- rbind(young, old)
AGE$network <- c(rep('young', nrow(young)), rep('old', nrow(old)))
AGE.pca <- prcomp(AGE[,c(2:7, 12:13, 16, 20:21)], scale = T)

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
WOMEN <- rbind(young, old)
rownames(WOMEN) <- c(1:21, '1 ', '2 ', '3 ', '4 ', '5 ', '6 ', '7 ', '8 ', '9 ', 
                   '10 ', '11 ', '12 ', '13 ', '14 ', '15 ', '16 ', '17 ', 
                   '18 ', '19 ', '20 ', '21 ')
WOMEN$network <- c(rep('young women', nrow(young.women)), rep('old women', nrow(old.women)))
WOMEN.pca <- prcomp(WOMEN[,c(2:7, 12:13, 16, 20:21)], scale = T)

###############################################################################

# Create figures

# Make figures that shape network and colour lipid type
# autoplot(SEX.pca, data = SEX, colour = 'group', main = 'Sex', shape = 'network', size = 3) + 
#   scale_colour_manual(values=c("#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC")) +
# theme_bw() + 
#   theme(legend.title = element_blank())
# autoplot(MEN.pca, data = MEN, colour = 'group', main = 'Men', shape = 'network', size = 3) + 
#   scale_colour_manual(values=c("#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC")) +
# theme_bw() + 
#   theme(legend.title = element_blank())
# autoplot(WOMEN.pca, data = WOMEN, colour = 'group', main = 'Women', shape = 'network', size = 3) + 
#   scale_colour_manual(values=c("#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC")) +
# theme_bw() + 
#   theme(legend.title = element_blank())

# Make figures that colour network and number lipids
autoplot(SEX.pca, data = SEX, colour = 'network', label = T, shape = F) + 
  # Add a confidence interval ellipse
  stat_ellipse(level = 0.95, aes(group = SEX$network, color = SEX$network), type = "norm") +
  # Set the colours
  scale_colour_manual(values=c('blue', 'red')) +
  # Make the plot with a white background
  theme_bw(base_size = 17) + 
  # Remove the legend title
  theme(legend.title = element_blank())

autoplot(MEN.pca, data = MEN, colour = 'network', label = T, shape = F) + 
  # Add a confidence interval ellipse
  stat_ellipse(level = 0.95, aes(group = MEN$network, color = MEN$network), type = "norm") +
  # Set the colours
  scale_colour_manual(values=c('darkblue', 'cyan')) +
  # Make the plot with a white background
  theme_bw(base_size = 17) + 
  # Remove legend title
  theme(legend.title = element_blank())

autoplot(WOMEN.pca, data = WOMEN, colour = 'network', label = T, shape = F) + 
  # Add a confidence interval ellipse
  stat_ellipse(level = 0.95, aes(group = WOMEN$network, color = WOMEN$network), type = "norm") +
  # Set the colours
  scale_colour_manual(values=c('red', 'orange')) +
  # Make the plot with a white background
  theme_bw(base_size = 17) + 
  # Remove legend title
  theme(legend.title = element_blank())

###############################################################################
