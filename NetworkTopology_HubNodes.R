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


