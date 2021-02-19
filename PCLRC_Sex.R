###############################################################################

# Yasmijn Balder
# 11-02-2021

# PCLRC and GGM comparing (all) men to (all) women using main fractions

# Output
# Opens a network in cytoscape

###############################################################################

# Load packages
library(GeneNet)
library(cropcor)
library(longitudinal)
library(fdrtool)
library(minet)
library(qgraph)             # Used to create the network
library(igraph)             # Used to create the network
library(RCy3)               # Used to open the network in cytoscape
library(openxlsx)           # Used to write excel files
library(stringr)            # Used to split and shorten the lipid names
library(ggplot2)            # Used to make pretty plots

###############################################################################

# Functions from Edoardo and Maria
getHighestLinks <- function(matrix, rank.thr = 0.3, verbose = FALSE) {
  # Obtains the links with the hightes weight from an weighted adjacency matrix
  #
  # Arguments:
  #   weighedAdjacencyMatrix:  a weighted adjacency matrix (symmetric matrix)
  
  #   rank.thre: fraction of interactions to keep. Default = 0.3 
  #
  # Returns: 
  #   binary matrix with only the higher links having non null values
  
  if(max(matrix, na.rm = TRUE) == 0){
    th = 0.01
    if(verbose) cat("null clr matrix \n")
  } else{                                      # get the threshold
    if(is.null(rank.thr)) {
      th <- min(matrix[matrix > 0], na.rm = TRUE) 
    } else{
      th <- quantile((matrix[upper.tri(matrix)]), 1 - rank.thr,na.rm = TRUE) 
      if(th == 0){
        if(verbose) 
          cat("threshold too low, min of the (non-null) matrix chosen instead\n")
        th <- min(matrix[matrix>0])
      }
    }
  }
  # select the ones that are above the threshold...
  net <- matrix
  net[which(net < th)] = 0
  net[which(net >= th)] = 1
  return(net)
}   


PCLRC.gmm <- function(datamatrix, prob.threshold = 0.95, Niter = 1000, 
                      frac = 0.75, rank.thr = 0.3){
  
  # performs the PCLRC to infer probabilistic associations from a data matrix 
  # using a Gaussian Graphical model to estimate partial correlations
  #
  # Arguments:
  #   datamatrix:   numeric matrix with the objects measurements.
  #                 objects in rows, samples in cols.
  #   Niter:        Integer value with the number of iterations. Default is 1000
  #   frac:         Fraction of the samples to be considered at each iteration.
  #                 Default is 0.75 (75%)            
  #   rank.thr:     Fraction of the total predicted interactions to be kept at
  #                 each iteration. Default is 0.3 (30% highest scoring interations kept).     
  # Returns:  
  #
  #     Correlation matrix of datamatrix
  #     Filtered Correlation matrix of datamatrix
  #     Matrix of probabilities obtained using PCLRC
  #
  #     Requires GeneNet package
  #
  #     Saccenti et al 2015 J. Proteome Res., 2015, 14 (2)
  
  require(GeneNet)
  
  if(prob.threshold > 1 || prob.threshold < 0)
    stop('prob.threshold must be between 0 and 1')
  
  if(frac > 1 || frac < 0)
    stop('frac must be between 0 and 1. Default 0.75')
  
  if(rank.thr > 1 || rank.thr  < 0)
    stop('frac must be between 0 and 1. Default 0.35')
  
  nsamp = round(dim(datamatrix)[1] * frac) # number of samples per iteration
  Nvaliter = 0  # number of valid iterations (no NAs generated)
  
  # this table will store the number of times an interaction was selected
  table <- mat.or.vec(dim(datamatrix)[2], dim(datamatrix)[2])
  
  
  for(it in (1:Niter)){
    # randomly select the samples
    samples = sample(dim(datamatrix)[1], nsamp)
    
    # build similarity matrix based on GMM from GeneNet
    
    similarityMatrixSubset <- ggm.estimate.pcor(as.matrix(datamatrix[samples,]), 
                                                method = "static",verbose = F)
    
    similarityMatrixSubset = as.matrix(unclass(similarityMatrixSubset))
    similarityMatrixSubset <- similarityMatrixSubset^2
    
    # infer network for this iteration using CLR
    
    adjSubset = clr(similarityMatrixSubset)
    
    if(!is.na(sum(adjSubset))){ #valid iteration
      #extract highest links
      out = getHighestLinks(adjSubset, rank.thr)
      #collect output
      table <- table+out
      Nvaliter = Nvaliter+1
    }
  }
  
  ProbMat  <- table / Nvaliter
  
  # Corr.mat = cor(datamatrix, method = "spearman")
  
  CorrMat = ggm.estimate.pcor(as.matrix(datamatrix), 
                              method = "static", verbose = F)
  CorrMat = as.matrix(unclass(CorrMat))
  
  CorrMat.Filtered  <- CorrMat
  CorrMat.Filtered[which(ProbMat < prob.threshold)] = 0
  diag(CorrMat.Filtered) = 1 #Force diagonal element equal to 1
  
  networks = list(CorrMatFiltered = CorrMat.Filtered, CorrMat = CorrMat, 
                  ProbMat = ProbMat, prob.threshold = prob.threshold)
  
  return(networks)
} 

Diff.Conn.PCLRC.gmm <- function(X1, X2, corr.type = 'pearson', 
                                prob.threshold = 0.95, adjust.diff = 'BH', 
                                Niter = 1000, MaxPerm = 1000, frac = 0.75, 
                                rank.thr = 0.3, verbose = FALSE){
  
  # Calculate differential connectivity between the two networks calculated
  # with PCLRC using partial correlation form data matrix X1 and X2
  #
  # Peform permutation test to assess statistical significance of the differential connectivity
  
  # Arguments:
  #   X1:             numeric matrix with the objects measurements.
  #   X2:             numeric matrix with the objects measurements.
  #   corr.type:      Type of correlation to be used
  #   prob.threshold: Probality treshold calculated by PCLRC on which to filter correlations. Default 0.95
  #   Niter:          Integer value with the number of iterations. Default is 1000
  #   frac:           Fraction of the samples to be considered at each iteration.
  #                   Default is 0.75 (75%)            
  #   rank.thr:       Fraction of the total predicted interactions to be kept at
  #                   each iteration. Default is 0.3 (30% highest scoring interations kept).
  #   Nperm:          Number of permutation for permutation test. Default 1000
  #   adjust.diff:    Pvalue correction for multiple testing. It uses p.adjust from stats package
  #                   Options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  #                   "fdr", "none"). Benjamini-Hochberg is default metohd.
  #  
  #   Returns:  
  #
  #   Connectivity C1 of network 1 from X1
  #   Connectivity C2 of network 2 from X2
  #   Differential Connectivity C1 - C2
  #   Adjacency matrix AdjMat1 for data X1
  #   Adjacency matrix AdjMat2 for data Xw
  #   P-values for each connectivity
  #   Matrix of permutations
  
  #
  
  if(ncol(X1) != ncol(X2)) 
    stop('X1 and X2 must have the same number of columns!')
  
  RES1 =   PCLRC.gmm(datamatrix = X1,  prob.threshold = prob.threshold, 
                     Niter = Niter, frac = frac, rank.thr = rank.thr)
  RES2 =   PCLRC.gmm(datamatrix = X2,  prob.threshold = prob.threshold, 
                     Niter = Niter, frac = frac, rank.thr = rank.thr)
  
  
  #CONNECTIVITY based on PCLRC
  W1C = RES1$CorrMatFiltered
  W2C = RES2$CorrMatFiltered
  
  G1C = abs(W1C) ##W1 is the correlation or mutual info for data set 1
  G2C = abs(W2C) ##W2 is the correlation or mutual info for data set 2
  
  Conn1_0C = rowSums(G1C) - 1
  Conn2_0C = rowSums(G2C) - 1
  
  
  #Calculate Differential connectivity
  diff0C = abs(Conn1_0C - Conn2_0C)
  
  
  
  ## PERMUTATION STARTS
  
  print('Entering permutation test')
  
  DIFFP_C = NULL
  DIFFP_M = NULL
  
  
  n1 = dim(X1)[1]
  n2 = dim(X2)[1]
  
  for (ii in 1:MaxPerm) {
    if(verbose)
      print(sprintf('Permutation %i of %i', ii, MaxPerm))
    X1p = NULL
    X2p = NULL
    for (jj in 1:dim(X1)[2])  {
      idx1p = sample(n1,n1)
      idx2p = sample(n2,n2)
      X1p = cbind(X1p,X1[idx1p,jj])
      X2p = cbind(X2p,X2[idx2p,jj])
    }
    
    # Calculate Corr
    
    RES1p = PCLRC.gmm(X1p, prob.threshold = prob.threshold, Niter = Niter, 
                      frac = frac, rank.thr = rank.thr)
    RES2p = PCLRC.gmm(X2p, prob.threshold = prob.threshold, Niter = Niter, 
                      frac = frac, rank.thr = rank.thr)
    
    
    W1pC = RES1p$CorrMatFiltered 
    W2pC = RES2p$CorrMatFiltered 
    
    
    # Calculate Connectivity based on CORR
    Conn_1C = rowSums(abs(W1pC)) - 1 
    Conn_2C = rowSums(abs(W2pC)) - 1 
    
    diffp_C = abs(Conn_1C - Conn_2C)
    DIFFP_C = cbind(DIFFP_C, diffp_C)
    
  }
  
  
  # Calculate P-value
  PVAL_C = NULL
  
  for(k in 1:dim(DIFFP_C)[1]) {
    
    LLC = length(which(DIFFP_C[k,] > diff0C[k]))
    PVAL_C[k] = (1 + LLC) / MaxPerm
  }
  
  # Adjust p-value of differential connectivity
  PVAL_C_ADJ = p.adjust(PVAL_C, method = adjust.diff)
  
  
  names(Conn1_0C) = colnames(X1)
  names(Conn2_0C) = colnames(X1)
  names(diff0C) = colnames(X1)
  names(PVAL_C) = colnames(X1)
  names(PVAL_C) = colnames(X1)
  names(PVAL_C_ADJ) = colnames(X1)
  
  results = list(Conn_X1 = Conn1_0C, Conn_X2 = Conn2_0C, Diff_Conn = diff0C, 
                 Pval = PVAL_C, Pval_adj = PVAL_C_ADJ, AdjMat1 = W1C, 
                 AdjMat2 = W2C, Permutations = DIFFP_C, 
                 adjust.diff =  adjust.diff)
  
  
  return(results)
  
}

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/Lipids_age_sex.csv", check.names = FALSE)
data[,1] <- NULL

# Split the data based on sex
men <- data[which(data$Gender=='man'),]
women <- data[which(data$Gender=='woman'),]

# Perform PCLRC + ggm
sex.pclrc <- Diff.Conn.PCLRC.gmm(men[,23:43], women[,23:43], verbose = TRUE)

###############################################################################

# Create figure

# Create a dataframe with the adjusted p-values
adjusted.pvalues <- as.data.frame(sex.pclrc$Pval_adj)
# Change the name of the p-values column
colnames(adjusted.pvalues)[1] <- 'P.values'
# Shorten the names of the lipids
short.names <- NULL
for (long.name in rownames(adjusted.pvalues)) {
  short.name <- str_split(long.name, ', ')[[1]][2:3]
  combined <- paste(short.name, collapse = ' & ')
  short.names <- c(short.names, combined)
}
# Add the lipid names as a column
adjusted.pvalues$lipids <- short.names
# Add the differential connectivity as a column
adjusted.pvalues$diffcon <- sex.pclrc$Diff_Conn
# Create a column with significance based on the adjusted p-values
adjusted.pvalues$sig <- rep("Not significant", 21)
adjusted.pvalues$sig[which(adjusted.pvalues$P.values <= 0.05)] <- 'Significant'

# Create a bar plot
ggplot(data = adjusted.pvalues, aes(x = diffcon, y = lipids, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("orange", "steelblue")) + 
  # Change the axis labels
  xlab('Differenrial connectivity') + ylab('Lipoprotein Main Fractions') + 
  # Remove the legend title
  theme(legend.title = element_blank())

###############################################################################

# Make network for men

# Retrieve the adjacency matrix
men_adj <- sex.pclrc$AdjMat1

# Save the adjacency matrix for COVSCA
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/COVSCA/")
write.xlsx(men_adj, 'Adjacency_matrix_men.xlsx')

# Remove disconnected nodes
disconnected.nodes <- which(apply(men_adj, 1, function(x){all(x==0)}))
if (length(disconnected.nodes)!=0) {
  men_adj <- men_adj[-disconnected.nodes,-disconnected.nodes]
  # groups <- groups[-disconnected.nodes]
  # shapes <- shapes[-disconnected.nodes]
  # colours <- colours[-disconnected.nodes]
}

# Create a qgraph with layout options
men_qgraph <- qgraph(input=men_adj,
                         labels=colnames(men_adj),
                         # groups=groups,
                         DoNotPlot=TRUE,
                         borders=FALSE,
                         palette="colorblind",
                         label.font='sans',
                         posCol="#009E73",  # colour of positive edges
                         negCol="#D55E00",  # colour of negative edges
                         # color=colours,     # colour of groups
                         # shape=shapes,      # shapes of groups
                         fade=FALSE,        # edge transparency based on weight
                         esize=2)

# Convert qgraph to igraph object
men_igraph <- as.igraph(men_qgraph, attributes = TRUE)
V(men_igraph)$name <- colnames(men_adj)

# Connect to cytoscape (Make sure cytoscape is opened)
cytoscapePing()
# Create the network
createNetworkFromIgraph(igraph = men_igraph,
                        title=paste0("Men"),
                        collection="PCLRC Men vs. Women")

###############################################################################

# Make network for women

# Retrieve the adjancy matrix
women_adj <- sex.pclrc$AdjMat2

# Save the adjacency matrix for COVSCA
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/COVSCA/")
write.xlsx(women_adj, 'Adjacency_matrix_women.xlsx')

disconnected.nodes <- which(apply(women_adj, 1, function(x){all(x==0)}))
if (length(disconnected.nodes)!=0) {
  women_adj <- women_adj[-disconnected.nodes,-disconnected.nodes]
  # groups <- groups[-disconnected.nodes]
  # shapes <- shapes[-disconnected.nodes]
  # colours <- colours[-disconnected.nodes]
}

women_qgraph <- qgraph(input=women_adj,
                       labels=colnames(women_adj),
                       # groups=groups,
                       DoNotPlot=TRUE,
                       borders=FALSE,
                       palette="colorblind",
                       label.font='sans',
                       posCol="#009E73",  # colour of positive edges
                       negCol="#D55E00",  # colour of negative edges
                       # color=colours,     # colour of groups
                       # shape=shapes,      # shapes of groups
                       fade=FALSE,        # edge transparency based on weight
                       esize=2)

women_igraph <- as.igraph(women_qgraph, attributes = TRUE)
V(women_igraph)$name <- colnames(women_adj)

createNetworkFromIgraph(igraph = women_igraph,
                        title=paste0("Women"),
                        collection="PCLRC Men vs. Women")

###############################################################################
