###############################################################################

# Yasmijn Balder
# 11-02-2021

# PCLRC and GGM comparing (all) men to (all) women using main fractions

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
# Function from Sanjee
VisualiseNetwork <- function(A, Group = TRUE, G, type = 1) {
  #Required packages:
  #   library(dplyr)
  #   library(RCy3)
  #Keep cytoscape open in the background
  #Inputs
  #   A is an adjacency matrix with the 1st column containing the column names and 1st row containing the row names and 
  #       the matrix containing numerical values corresponding to the variable relationships in the columns and rows
  #   A can also be a list of matrices.
  #   Group is boolean asking if certain variables are to be clustered together as neighbours in the figure. 
  #   Incase certain variables are to be clustered together in groups in the visualisation,
  #   G is a vector containing the label names of the "Groups to be clustered together" 
  #For visualising the Edge Weights,
  #   Choose "type = 1" for grading the edges based on partial correlation values
  #   Choose "type = 2" for grading the edges based on Pearson or Spearman correlation values
  #   Choose "type = 3" for grading the edges on a ranked percentile system (such that the edges are ranked and on an exponential scale
  #   the gradient of the width and colour of edges are assigned. for eg. 98th percentile with the highest width, 95th percentile with the next width etc.)
  #
  #
  #Example:
  #	source("VisualiseNetwork.R")
  #	library(dplyr)
  #	library(RCy3)
  #
  #	setwd("full path to working directory")
  #	Mat1 <- read.table("Adjacency_Matrix1.txt", header = TRUE, row.names = 1)
  #       Mat2 <- read.table("Adjacency_Matrix2.txt", header = TRUE, row.names = 1)
  #       mat3 <- read.table("Adjacency_Matrix3.txt", header = TRUE, row.names = 1)
  #
  #	A <- list(Mat1, Mat2, Mat3)
  #	G <- as.vector(c(label1, label2, label3...., labeln))
  #
  #	Visual <- VisualiseNetwork(A, Group = TRUE, G, type = 3)
  
  require(dplyr)
  while(inherits(A, "data.frame") == TRUE || inherits(A, "matrix") == TRUE){A <- list(A)}
  
  AdjMatrix = NULL
  NodesNetwork = NULL
  EdgesNetwork = NULL
  
  for (matrix in 1:length(A)) {
    
    Adjacency = NULL
    NodeTable = NULL
    EdgeTable = NULL    
    
    Adjacency = as.data.frame(A[matrix])
    
    if(ncol(Adjacency) != nrow(Adjacency)) stop("Adjacency matrix should be a square matrix with equal number of rows and columns")
    if(nrow(Adjacency) != length(G)) stop("The number of nodes/variables in the groups table should be the same as in the adjacency matrix")
    
    
    
    if(Group == FALSE){
      Adj <- Adjacency
      Node <- as.vector(colnames(Adj))
      Node <- as.data.frame(Node)
      Groups <- as.vector(rep("A", nrow(Adj)))
      NodeTable <- as.data.frame(cbind(Node, Groups))
    } else {
      Groups = as.vector(G)
      Adj <- Adjacency
      Node <- as.vector(colnames(Adj))
      Node <- as.data.frame(Node)
      NodeTable <- as.data.frame(cbind(Node, Groups))
    }
    #NodeTable <- left_join(Node, Node_Group, by = NULL,  copy = FALSE, suffix=c(".n",".g"))
    NodeTable <- with(NodeTable, NodeTable[order(NodeTable$Groups) , ])
    
    diag(Adj) = 0
    Adj[lower.tri(Adj, diag=TRUE)] <- 0
    
    Source = NULL
    Target = NULL
    Weight = NULL
    for (row in 1:nrow(Adj)) {
      for (col in 1:ncol(Adj)) {
        if (Adj[row, col] != 0) {
          Source <- as.vector(append(Source, rownames(Adj[row, ])))
          Target <- as.vector(append(Target, colnames(Adj[col])))
          Weight <- as.vector(append(Weight, as.numeric(Adj[row, col])))
        } else {}
      }
    }
    Interaction <- as.vector(rep("interacts", length(Weight)))
    EdgeTable <- as.data.frame(cbind(Source, Target, Weight, Interaction))
    EdgeTable$Weight <- as.numeric(as.character(EdgeTable$Weight))
    
    
    X = NULL
    Y = NULL
    R = round(nrow(Adj)/10, 0) * (100)
    
    for (i in 0:(nrow(Adj) - 1)) {
      print(i)
      x = R*cos((i*2*3.14159265359)/(nrow(Adj)))
      X <- as.vector(append(X, x)) 
      y = R*sin((i*2*3.14159265359)/(nrow(Adj)))
      Y <- as.vector(append(Y, y))
    }
    pos <- as.data.frame(cbind(X,Y))
    
    NodeTable <- cbind(NodeTable, pos)
    
    frac = as.vector(c(2, 3, 4, 6, 10, 15, 24, 36))
    E = nrow(EdgeTable)
    M = NULL
    for (i in 1:length(frac)) {
      f = frac[i]
      mes = NULL
      mes = round((f * E)/100, 0)
      M <- as.vector(append(M, mes))
    }
    diff = (sum(M)) - (nrow(EdgeTable))
    ifelse(diff == 0, print("perfect!"), M[8] <- M[8] - diff)
    
    wids <- as.vector(c(10, 8, 4, 2, 1, 0.5, 0.25, 0.25))
    wid = NULL
    for (j in 1:length(M)) {
      times = M[j]
      value = wids[j]
      wid <- as.vector(append(wid, c(rep(value, times))))
    }
    
    ifelse(type == 1, EdgeTable <- mutate(EdgeTable, width = ifelse(Weight > 0.5, 10, ifelse(Weight < -0.5, 10, ifelse(Weight < 0.5 & Weight > 0.4, 8, ifelse(Weight > -0.5 & Weight < -0.4, 8, ifelse(Weight < 0.4 & Weight > 0.3, 4, ifelse(Weight > -0.4 & Weight < -0.3, 4, ifelse(Weight < 0.3 & Weight > 0.2, 2, ifelse(Weight > -0.3 & Weight < -0.2, 2, ifelse(Weight < 0.2 & Weight > 0.15, 1, ifelse(Weight > -0.2 & Weight < -0.15, 1, ifelse(Weight < 0.15 & Weight > 0.1, 0.5, ifelse(Weight > -0.15 & Weight < -0.1, 0.5, 0.25))))))))))))), ifelse(type == 2, EdgeTable <- mutate(EdgeTable, width = ifelse(Weight > 0.9, 10, ifelse(Weight < -0.9, 10, ifelse(Weight < 0.9 & Weight > 0.8, 8, ifelse(Weight > -0.9 & Weight < -0.8, 8, ifelse(Weight < 0.8 & Weight > 0.7, 4, ifelse(Weight > -0.8 & Weight < -0.7, 4, ifelse(Weight < 0.7 & Weight > 0.6, 2, ifelse(Weight > -0.7 & Weight < -0.6, 2, ifelse(Weight < 0.6 & Weight > 0.5, 1, ifelse(Weight > -0.6 & Weight < -0.15, 1, ifelse(Weight < 0.5 & Weight > 0.4, 0.5, ifelse(Weight > -0.5 & Weight < -0.4, 0.5, 0.25))))))))))))), if(type == 3){
      wid <- as.data.frame(wid)
      EdgeTable <- EdgeTable[sort(abs(EdgeTable$Weight), decreasing=T, index.return=T)[[2]],]
      EdgeTable <- cbind(EdgeTable, wid)
      colnames(EdgeTable)[5] <- "width"
    } else{
      print("type not selected")}))
    
    ifelse(type == 1, EdgeTable <- mutate(EdgeTable, Stroke = ifelse(Weight > 0.5, "#DC1C13", ifelse(Weight < -0.5, "#1F1FFF", ifelse(Weight < 0.5 & Weight > 0.4, "#EA4C46", ifelse(Weight > -0.5 & Weight < -0.4, "#4949FF", ifelse(Weight < 0.4 & Weight > 0.3, "#F07470", ifelse(Weight > -0.4 & Weight < -0.3, "#7879FF", ifelse(Weight < 0.3 & Weight > 0.2, "#F1959B", ifelse(Weight > -0.3 & Weight < -0.2, "#A3A3FF", ifelse(Weight < 0.2 & Weight > 0.15, "#F6BDC0", ifelse(Weight > -0.2 & Weight < -0.15, "#BFBFFF", ifelse(Weight < 0.15 & Weight > 0.1, "#F6BDC0", ifelse(Weight > -0.15 & Weight < -0.1, "#BFBFFF", "#A6A6A6"))))))))))))), ifelse(type == 2, EdgeTable <- mutate(EdgeTable, Stroke = ifelse(Weight > 0.9, "#DC1C13", ifelse(Weight < -0.9, "#1F1FFF", ifelse(Weight < 0.9 & Weight > 0.8, "#EA4C46", ifelse(Weight > -0.9 & Weight < -0.8, "#4949FF", ifelse(Weight < 0.8 & Weight > 0.7, "#F07470", ifelse(Weight > -0.8 & Weight < -0.7, "#7879FF", ifelse(Weight < 0.7 & Weight > 0.6, "#F1959B", ifelse(Weight > -0.7 & Weight < -0.6, "#A3A3FF", ifelse(Weight < 0.6 & Weight > 0.5, "#F6BDC0", ifelse(Weight > -0.6 & Weight < -0.5, "#BFBFFF", ifelse(Weight < 0.5 & Weight > 0.4, "#F6BDC0", ifelse(Weight > -0.5 & Weight < -0.4, "#BFBFFF", "#A6A6A6"))))))))))))), ifelse(type == 3, EdgeTable <- mutate(EdgeTable, Stroke = ifelse(width == 10 & Weight > 0, "#DC1C13", ifelse(width == 10 & Weight < 0, "#1F1FFF", ifelse(width == 8 & Weight > 0, "#EA4C46", ifelse(width == 8 & Weight < 0, "#4949FF", ifelse(width == 4 & Weight > 0, "#F07470", ifelse(width == 4 & Weight < 0, "#7879FF", ifelse(width == 2 & Weight > 0, "#F1959B", ifelse(width == 2 & Weight < 0, "#A3A3FF", ifelse(width == 1 & Weight > 0, "#F6BDC0", ifelse(width == 1 & Weight < 0, "#BFBFFF", ifelse(width == 0.5 & Weight > 0, "#F6BDC0", ifelse(width == 0.5 & Weight < 0, "#BFBFFF", "#A6A6A6"))))))))))))), "type not selected")))
    
    
    EdgeTable$sharedname <- paste(EdgeTable$Source, "(interacts)", EdgeTable$Target)
    
    Network_name = sprintf("Visual_Network_%i", matrix)
    Network_Collection = sprintf("Visual_Networks_%i", matrix)
    
    nodes = NULL
    edges = NULL
    
    nodes <- data.frame(id=as.vector(NodeTable$Node), group=as.vector(NodeTable$Groups), stringsAsFactors = FALSE)
    edges <- data.frame(source=as.vector(EdgeTable$Source), target=as.vector(EdgeTable$Target), interaction=as.vector(EdgeTable$Interaction), weight=as.vector(EdgeTable$Weight), stringsAsFactors = FALSE)
    
    createNetworkFromDataFrames(nodes, edges, title=Network_name, collection=Network_Collection, style.name = "SanjeeNetworkStyle")
    
    
    Colour_palette <- as.vector(c("#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC", "#003C6799", "#8F7700", "#3B3B3B", "#A73030", "#4A6990"))
    
    style.name = "SanjeeNetworkStyle"
    defaults <- list(NODE_SHAPE="Ellipse", NODE_SIZE=25.0, EDGE_TRANSPARENCY=255 , NODE_LABEL_POSITION="W,E,c,0.00,0.00", NODE_BORDER_PAINT="#FFFFFF")
    nodeLabels <- mapVisualProperty("Node Label", "id", "p")
    nodecolour <- mapVisualProperty("Node Fill Color", "group", "d", as.vector(unique(NodeTable$Groups)), as.vector(Colour_palette[1:length(unique(NodeTable$Groups))]))
    nodeXlocation <- mapVisualProperty("Node X Location", "id", "d", as.vector(NodeTable$Node), as.vector(NodeTable$X))
    nodeYlocation <- mapVisualProperty("Node Y Location", "id", "d", as.vector(NodeTable$Node), as.vector(NodeTable$Y))
    edgeline <- mapVisualProperty("Edge Line Type", "interaction", "d", as.vector(unique(EdgeTable$Interaction)), as.vector(c("Solid")))
    edgewidth <- mapVisualProperty("Edge Width", "shared name", "d", as.vector(EdgeTable$sharedname), as.vector(EdgeTable$width))
    edgestroke <- mapVisualProperty("Edge Stroke Unselected Paint", "shared name", "d", as.vector(EdgeTable$sharedname), as.vector(EdgeTable$Stroke))
    
    createVisualStyle(style.name, defaults, list(nodeLabels, nodecolour, nodeXlocation, nodeYlocation, edgeline, edgewidth, edgestroke))
    setVisualStyle("SanjeeNetworkStyle")
    
    fitContent(selected.only = FALSE)
    fitContent(selected.only = FALSE)
    fitContent(selected.only = FALSE)
    
    Network_out = sprintf("Network_Image_%i", matrix)
    
    full.path = paste(getwd(), Network_out, sep="/")
    exportImage(full.path, "PNG", units="pixels", width=3600, height=1771)
    
    #Network files for building network using some other software
    AdjMatrix <- list(AdjMatrix, Adjacency)
    NodesNetwork <- list(NodesNetwork, NodeTable)
    EdgesNetwork <- list(EdgesNetwork, EdgeTable)
    
    Network_save = sprintf("Cytoscape_Network_%i", matrix)
    full.path.cps = paste(getwd(), Network_save, sep="/")
    closeSession(save.before.closing = TRUE, filename = full.path.cps)
    
  }
  
  
  Network = list(AdjMatrix, NodesNetwork, EdgesNetwork)
  return(Network)
  
}

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/Lipids_age_sex.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

# Split the data based on sex
men <- data[which(data$Gender=='man'),]
women <- data[which(data$Gender=='woman'),]

# Perform PCLRC
sex.pclrc <- Diff.Conn.PCLRC.gmm(men[,23:43], women[,23:43], verbose = TRUE, 
                                 adjust.diff = 'bonferroni',
                                 prob.threshold = 0.99)

###############################################################################

# Make groups of lipoprotein main fractions
groups <- as.vector(c(rep('Triglycerides', 4), rep('Cholesterol', 4), 
                      rep('FreeCholesterol', 4), rep ('Phospholipids', 4), rep ('Apo', 5)))

# Retrieve the adjacency matrix
# men_adj <- as.data.frame(sex.pclrc$AdjMat1)
# women_adj <- as.data.frame(sex.pclrc$AdjMat2)
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/")
young_adj <- read.csv('Adjacency_matrix_young.csv')
young_adj <- young_adj[,-1]
old_adj <- read.csv('Adjacency_matrix_old.csv')
old_adj <- old_adj[,-1]
rownames(men_adj) <- colnames(men_adj) <- 
  rownames(women_adj) <- colnames(women_adj) <- c("Triglycerides_VLDL", 
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
men_network <- VisualiseNetwork(A = men_adj, Group = TRUE, G = groups)
women_network <- VisualiseNetwork(A = women_adj, Group = TRUE, G = groups)

# Save adjacency matrices
# write.csv(men_adj, file = 'Adjacency_matrix_men.csv')
# write.csv(women_adj, file = 'Adjacency_matrix_women.csv')

###############################################################################

# Create figure

# Create a dataframe with the p-values
pvalues <- as.data.frame(sex.pclrc$Pval)
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
pvalues$diffcon <- sex.pclrc$Diff_Conn

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
setwd("C:/Users/Yasmijn/Pictures/Research practice")
for (i in 1:7) {
  g <- ggplot(data = pvalues, aes_string(x = nm[3], y = nm[2], fill = nm[5+i])) + 
    geom_bar(stat = "identity") +
    # Add a title
    ggtitle(paste('Men vs Women\nMT =', testing[i], '\np-value threshold =', thresholds[i])) +
    # Use nicer colours
    scale_fill_manual(values = c("orange", "steelblue")) + 
    # Change the axis labels
    xlab('Differential connectivity') + ylab('Lipoprotein Main Fractions') + 
    # Remove the legend title
    theme(legend.title = element_blank())
  print(g)
  ggsave(filename = paste0('PCLRC_Sex_',testing[i], '_p-threshold=', thresholds[i], '.png'), 
         plot = g)
}

# First: CHOOSE MULTIPLE TESTING
# Create table; 
# summary.table <- data.frame(row.names = short.names)
# summary.table$Age <- rep('', 21)
# summary.table$Age[which(pvalues$sig == 'Significant' & pvalues$diffcon > 0)] <- '+'
# summary.table$Age[which(pvalues$sig == 'Significant' & pvalues$diffcon < 0)] <- '-'
