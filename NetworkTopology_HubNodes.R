###############################################################################

# Yasmijn Balder
# 10-03-2021

# Finds hub nodes in the networks

# Output
# Shows a table with hub nodes given the name of the network and the node

###############################################################################

# Load packages

###############################################################################

# Load the data
men <- read.csv('Results/NetworkTopology_Men.csv')
women <- read.csv('Results/NetworkTopology_Women.csv')

young.men <- read.csv('Results/NetworkTopology_YoungMen.csv')
old.men <- read.csv('Results/NetworkTopology_OldMen.csv')

young.women <- read.csv('Results/NetworkTopology_YoungWomen.csv')
old.women <- read.csv('Results/NetworkTopology_OldWomen.csv')

###############################################################################

# Save hub nodes in a data frame
hub.nodes <- data.frame()

# Make list of all networks
networks <- list('men' = men, 'women' = women, 
                 'young.men' = young.men, 'old.men' = old.men, 
                 'young.women' = young.women, 'old.women' = old.women)

###############################################################################

# Which change from < 0.03 to > 0.03?
cc.thrs <- data.frame(row.names = women$name)
cc.thrs$sex <- rep('', 21)
cc.thrs$sex[which(women$ClusteringCoefficient < 0.03 & men$ClusteringCoefficient > 0.03)] <- 'women'
# cc.thrs$sex[which(women$ClusteringCoefficient > 0.03 & men$ClusteringCoefficient < 0.03)] <- 'men'
cc.thrs$men <- rep('', 21)
# cc.thrs$men[which(young.men$ClusteringCoefficient < 0.03 & old.men$ClusteringCoefficient > 0.03)] <- 'young'
cc.thrs$men[which(young.men$ClusteringCoefficient > 0.03 & old.men$ClusteringCoefficient < 0.03)] <- 'old'
cc.thrs$women <- rep('', 21)
# cc.thrs$women[which(young.women$ClusteringCoefficient < 0.03 & old.women$ClusteringCoefficient > 0.03)] <- 'young'
cc.thrs$women[which(young.women$ClusteringCoefficient > 0.03 & old.women$ClusteringCoefficient < 0.03)] <- 'old'

View(cc.thrs)

