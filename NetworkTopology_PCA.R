###############################################################################

# Yasmijn Balder
# 10-03-2021

# PCA of network topology

# Output
# ...

###############################################################################

# Load packages
library(ggfortify)      # Used to plot the PCA

###############################################################################

# Load the data
men <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_Men.csv')
women <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_Women.csv')

young <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_Young.csv')
old <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_Old.csv')

young.men <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_YoungMen.csv')
old.men <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_OldMen.csv')

young.women <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_YoungWomen.csv')
old.women <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_OldWomen.csv')

###############################################################################

# Make a PCA for age
AGE <- rbind(young, old)
AGE$network <- c(rep('young', nrow(young)), rep('old', nrow(old)))
AGE.pca <- prcomp(AGE[,c(2:7, 12:13, 16, 20:21)], scale = T)
autoplot(AGE.pca, data = AGE, colour = 'network', main = 'Age', label = T, shape = F)

# Make a PCA for sex
SEX <- rbind(men, women)
SEX$network <- c(rep('men', nrow(men)), rep('women', nrow(women)))
SEX.pca <- prcomp(SEX[,c(2:7, 12:13, 16, 20:21)], scale = T)
autoplot(SEX.pca, data = SEX, colour = 'network', main = 'Sex', label = T, shape = F)

# Make a PCA for men
MEN <- rbind(young.men, old.men)
MEN$network <- c(rep('young men', nrow(young.men)), rep('old men', nrow(old.men)))
MEN.pca <- prcomp(MEN[,c(2:7, 12:13, 16, 20:21)], scale = T)
autoplot(MEN.pca, data = MEN, colour = 'network', main = 'Men', label = T, shape = F)

# Make a PCA for women
WOMEN <- rbind(young, old)
WOMEN$network <- c(rep('young women', nrow(young.women)), rep('old women', nrow(old.women)))
WOMEN.pca <- prcomp(WOMEN[,c(2:7, 12:13, 16, 20:21)], scale = T)
autoplot(WOMEN.pca, data = WOMEN, colour = 'network', main = 'Women', label = T, shape = F)

###############################################################################

