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
men <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_Men.csv')
women <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_Women.csv')

young <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_Young.csv')
old <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_Old.csv')

young.men <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_YoungMen.csv')
old.men <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_OldMen.csv')

young.women <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_YoungWomen.csv')
old.women <- read.csv('C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/NetworkTopology_OldWomen.csv')

###############################################################################

# Save hub nodes in a data frame
hub.nodes <- data.frame()

# Make list of all networks
networks <- list('men' = men, 'women' = women, 'young' = young, 'old' = old, 
                 'young.men' = young.men, 'old.men' = old.men, 
                 'young.women' = young.women, 'old.women' = old.women)

# Go through all networks
name <- 1
for (network in networks) {
  for (lipid in 1:nrow(network)) {
    # Get the mean and sd of the clustering coefficient in this network
    mean.clust <- mean(network$ClusteringCoefficient)
    sd.clust <- sd(network$ClusteringCoefficient)
    # Get the threshold for a hub by adding 1 sd to the mean
    thr.clust <- mean.clust - sd.clust
    # Get the mean and sd of the degree in this network
    mean.degree <- mean(network$Degree)
    sd.degree <- sd(network$Degree)
    # Get the threshold for a hub by adding 1 sd to the mean
    thr.degree <- mean.degree + sd.degree
    # Find the hub nodes
    if (network$ClusteringCoefficient[lipid] < thr.clust & network$Degree[lipid] > thr.degree) {
      print('Found a hub node!')
      # Add the hub node to the table
      hub.nodes <- rbind(hub.nodes, c(names(networks)[name], network$name[lipid]))
    }
  }
  name <- name+ 1
}

# Give the hub node table column names
colnames(hub.nodes) <- c('Network','Lipid')

View(hub.nodes)

###############################################################################

# What is the difference between the cc's per network?
cc.diff <- data.frame(row.names = women$name)
cc.diff$sex <- women$ClusteringCoefficient - men$ClusteringCoefficient
cc.diff$age <- young$ClusteringCoefficient - old$ClusteringCoefficient
cc.diff$men <- young.men$ClusteringCoefficient - old.men$ClusteringCoefficient
cc.diff$women <- young.women$ClusteringCoefficient - old.women$ClusteringCoefficient

View(cc.diff)

# Which change from < 0.03 to > 0.03?
cc.thrs <- data.frame(row.names = women$name)
cc.thrs$sex <- rep('', 21)
cc.thrs$sex[which(women$ClusteringCoefficient < 0.03 & men$ClusteringCoefficient > 0.03)] <- 'women'
cc.thrs$sex[which(women$ClusteringCoefficient > 0.03 & men$ClusteringCoefficient < 0.03)] <- 'men'
cc.thrs$age <- rep('', 21)
cc.thrs$age[which(young$ClusteringCoefficient < 0.03 & old$ClusteringCoefficient > 0.03)] <- 'young'
cc.thrs$age[which(young$ClusteringCoefficient > 0.03 & old$ClusteringCoefficient < 0.03)] <- 'old'
cc.thrs$men <- rep('', 21)
cc.thrs$men[which(young.men$ClusteringCoefficient < 0.03 & old.men$ClusteringCoefficient > 0.03)] <- 'young'
cc.thrs$men[which(young.men$ClusteringCoefficient > 0.03 & old.men$ClusteringCoefficient < 0.03)] <- 'old'
cc.thrs$women <- rep('', 21)
cc.thrs$women[which(young.women$ClusteringCoefficient < 0.03 & old.women$ClusteringCoefficient > 0.03)] <- 'young'
cc.thrs$women[which(young.women$ClusteringCoefficient > 0.03 & old.women$ClusteringCoefficient < 0.03)] <- 'old'

View(cc.thrs)
