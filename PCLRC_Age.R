###############################################################################

# Yasmijn Balder
# 11-02-2021

# PCLRC and GGM comparing young to old using main fractions

# Output
# Opens a network in cytoscape

###############################################################################

# Load packages
library(GeneNet)
library(corpcor)
library(longitudinal)
library(fdrtool)
library(minet)
library(qgraph)             # Used to create the network
library(igraph)             # Used to create the network
library(RCy3)               # Used to open the network in cytoscape
library(openxlsx)           # Used to write excel files
library(stringr)            # Used to split and shorten the lipid names
library(ggplot2)            # Used to make pretty plots
library(dplyr)
library(qvalue)             # Used to calculate storey's q-value

###############################################################################

# Functions from Edoardo, Maria and, Sanjeevan
source('../AllNetworkFunctions28-Mar-2021 16.49.R')
source('../VisualiseNetwork.R')

###############################################################################

# Load data
data <- read.csv("../Data/LipidsAgeSex_SqrtNormalization.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

# Split the data based on sex
young <- data[which(data$Age < quantile(data$Age, probs = 1/3)),]
old <- data[which(data$Age > quantile(data$Age, probs = 2/3)),]

# Perform PCLRC
age.pclrc <- Diff.Conn.PCLRC.gmm(young[,23:43], old[,23:43], verbose = TRUE, 
                                 adjust.diff = 'BH',
                                 prob.threshold = 0.99)

###############################################################################

# Save the PCLRC output
saveRDS(age.pclrc, 'Results/age.pclrc.rds')

# Save as xlsx for COVSCA
write.xlsx(age.pclrc$AdjMat1, 'COVSCA/Adjacency_matrix_young.xlsx')
write.xlsx(age.pclrc$AdjMat2, 'COVSCA/Adjacency_matrix_old.xlsx')

###############################################################################

# Make groups of lipoprotein main fractions
groups <- as.vector(c(rep('Triglycerides', 4), rep('Cholesterol', 4), 
            rep('FreeCholesterol', 4), rep ('Phospholipids', 4), rep ('Apo', 5)))

# Retrieve the adjacency matrix
# young_adj <- as.data.frame(age.pclrc$AdjMat1)
# old_adj <- as.data.frame(age.pclrc$AdjMat2)
young_adj <- read.csv('Results/Adjacency_matrix_young.csv')
young_adj <- young_adj[,-1]
old_adj <- read.csv('Results/Adjacency_matrix_old.csv')
old_adj <- old_adj[,-1]
rownames(old_adj) <- colnames(old_adj) <- 
  rownames(young_adj) <- colnames(young_adj) <- c("Triglycerides_VLDL", 
                                                "Triglycerides_IDL", 
                                                "Triglycerides_LDL", 
                                                "Triglycerides_HDL", 
                                                "Cholesterol_VLDL", 
                                                "Cholesterol_IDL", 
                                                "Cholesterol_LDL",
                                                "Cholesterol_HDL",
                                                "FreeCholesterol_VLDL",
                                                "FreeCholesterol_IDL",
                                                "FreeCholesterol_LDL",
                                                "FreeCholesterol_HDL",
                                                "Phospholipids_VLDL",
                                                "Phospholipids_IDL",
                                                "Phospholipids_LDL",
                                                "Phospholipids_HDL",
                                                "ApoA1_HDL","ApoA2_HDL",
                                                "ApoB_VLDL","ApoB_IDL",        
                                                "ApoB_LDL")

# Visualise the networks in cytoscape (make sure it is opened)
young_network <- VisualiseNetwork(A = young_adj, Group = TRUE, G = groups, type = 3)
old_network <- VisualiseNetwork(A = old_adj, Group = TRUE, G = groups, type = 3)

# Save adjacency matrices
# write.csv(young_adj, file = 'Adjacency_matrix_young.csv')
# write.csv(old_adj, file = 'Adjacency_matrix_old.csv')

###############################################################################

# Create figure

# Create a dataframe with the p-values
pvalues <- as.data.frame(age.pclrc$Pval)
# Change the name of the p-values column
colnames(pvalues)[1] <- 'P.values'
# Shorten the names of the lipids
short.names <- NULL
for (long.name in rownames(pvalues)) {
  short.name <- str_split(long.name, ', ')[[1]][2:3]
  combined <- paste(short.name, collapse = ' & ')
  short.names <- c(short.names, combined)
}
# Add the lipid names as a column
pvalues$lipids <- short.names
# Add the differential connectivity as a column
pvalues$diffcon <- age.pclrc$Diff_Conn

# Use various multiple testing corrections
pvalues$BH <- p.adjust(pvalues$P.values, method = 'BH')
pvalues$bonferroni <- p.adjust(pvalues$P.values, method = 'bonferroni')

# Decide significance based on the adjusted p-values
pvalues$BH1 <- rep("Not significant", 21)
pvalues$BH1[which(pvalues$BH <= 0.01)] <- 'Significant'
pvalues$BH5 <- rep("Not significant", 21)
pvalues$BH5[which(pvalues$BH <= 0.05)] <- 'Significant'
pvalues$bonferroni1 <- rep("Not significant", 21)
pvalues$bonferroni1[which(pvalues$bonferroni <= 0.01)] <- 'Significant'
pvalues$bonferroni5 <- rep("Not significant", 21)
pvalues$bonferroni5[which(pvalues$bonferroni <= 0.05)] <- 'Significant'
pvalues$storey1 <- qvalue(pvalues$P.values, fdr.level = 0.01, lambda = 0)$significant
pvalues$storey1[which(pvalues$storey1 == 'TRUE')] <- 'Significant'
pvalues$storey1[which(pvalues$storey1 == 'FALSE')] <- 'Not significant'
pvalues$storey5 <- qvalue(pvalues$P.values, fdr.level = 0.05, lambda = 0)$significant
pvalues$storey5[which(pvalues$storey5 == 'TRUE')] <- 'Significant'
pvalues$storey5[which(pvalues$storey5 == 'FALSE')] <- 'Not significant'
pvalues$storey05 <- qvalue(pvalues$P.values, fdr.level = 0.005, lambda = 0)$significant
pvalues$storey05[which(pvalues$storey05 == 'TRUE')] <- 'Significant'
pvalues$storey05[which(pvalues$storey05 == 'FALSE')] <- 'Not significant'

# For loop for plots
testing <- c('BH', 'BH', 'bonferroni', 'bonferroni', 'storey', 'storey', 'storey')
thresholds <- c('0.01', '0.05', '0.01', '0.05', '0.01', '0.05', '0.005')
nm <- names(pvalues)
# Change working directory to save the plots
setwd("Results")
for (i in 1:7) {
  g <- ggplot(data = pvalues, aes_string(x = nm[3], y = nm[2], fill = nm[5+i])) + 
    geom_bar(stat = "identity") +
    # Add a title
    ggtitle(paste('Young vs Old\nMT =', testing[i], '\np-value threshold =', thresholds[i])) +
    # Use nicer colours
    scale_fill_manual(values = c("orange", "steelblue")) + 
    # Change the axis labels
    xlab('Differential connectivity') + ylab('Lipoprotein Main Fractions') + 
    # Remove the legend title
    theme(legend.title = element_blank())
  print(g)
  ggsave(filename = paste0('PCLRC_Age_',testing[i], '_p-threshold=', thresholds[i], '.png'), 
         plot = g)
}

# First: CHOOSE MULTIPLE TESTING
# Create table; 
# summary.table <- data.frame(row.names = short.names)
# summary.table$Age <- rep('', 21)
# summary.table$Age[which(pvalues$sig == 'Significant' & pvalues$diffcon > 0)] <- '+'
# summary.table$Age[which(pvalues$sig == 'Significant' & pvalues$diffcon < 0)] <- '-'
