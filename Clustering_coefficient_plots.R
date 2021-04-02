###############################################################################

# Yasmijn Balder
# 01-04-2021

# Clustering coefficient plots

# Output
# A plot for each group of interest with the clustering coefficient per lipid

###############################################################################

# Load the data
men <- read.csv('Results/NetworkTopology_Men.csv')
women <- read.csv('Results/NetworkTopology_Women.csv')

young.men <- read.csv('Results/NetworkTopology_YoungMen.csv')
old.men <- read.csv('Results/NetworkTopology_OldMen.csv')

young.women <- read.csv('Results/NetworkTopology_YoungWomen.csv')
old.women <- read.csv('Results/NetworkTopology_OldWomen.csv')

###############################################################################

men$sig <- rep('Not significant', 21)
men$sig[which(men$ClusteringCoefficient < 0.03)] <- 'Significant'

women$sig <- rep('Not significant', 21)
women$sig[which(women$ClusteringCoefficient < 0.03)] <- 'Significant'

young.men$sig <- rep('Not significant', 21)
young.men$sig[which(young.men$ClusteringCoefficient < 0.03)] <- 'Significant'

old.men$sig <- rep('Not significant', 21)
old.men$sig[which(old.men$ClusteringCoefficient < 0.03)] <- 'Significant'

young.women$sig <- rep('Not significant', 21)
young.women$sig[which(young.women$ClusteringCoefficient < 0.03)] <- 'Significant'

old.women$sig <- rep('Not significant', 21)
old.women$sig[which(old.women$ClusteringCoefficient < 0.03)] <- 'Significant'

###############################################################################

# Create a bar plots

ggplot(data = men, aes(x = ClusteringCoefficient, y = name, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("steelblue", "orange")) + 
  # Change the axis labels
  xlab('Clustering coefficient') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  # Set the x lim, so they are the same for each plot
  xlim(0, 1)

ggplot(data = women, aes(x = ClusteringCoefficient, y = name, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("steelblue", "orange")) + 
  # Change the axis labels
  xlab('Clustering coefficient') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  # Set the x lim, so they are the same for each plot
  xlim(0, 1)

ggplot(data = young.men, aes(x = ClusteringCoefficient, y = name, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("steelblue", "orange")) + 
  # Change the axis labels
  xlab('Clustering coefficient') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  # Set the x lim, so they are the same for each plot
  xlim(0, 1)

ggplot(data = old.men, aes(x = ClusteringCoefficient, y = name, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("steelblue", "orange")) + 
  # Change the axis labels
  xlab('Clustering coefficient') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  # Set the x lim, so they are the same for each plot
  xlim(0, 1)

ggplot(data = young.women, aes(x = ClusteringCoefficient, y = name, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("steelblue", "orange")) + 
  # Change the axis labels
  xlab('Clustering coefficient') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  # Set the x lim, so they are the same for each plot
  xlim(0, 1)

ggplot(data = old.women, aes(x = ClusteringCoefficient, y = name, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("steelblue", "orange")) + 
  # Change the axis labels
  xlab('Clustering coefficient') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  # Set the x lim, so they are the same for each plot
  xlim(0, 1)
