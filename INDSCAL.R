###############################################################################

# Yasmijn Balder
# 22-02-2021

# INDSCAL

# Output
# 

###############################################################################

# Load packages
library(smacof)             # Used to perform INDSCAL
library(openxlsx)           # Used to open excel files
library(stringr)            # Used to split the names

###############################################################################

# Load data
men <- read.xlsx("COVSCA/Adjacency_matrix_men.xlsx")
old <- read.xlsx("COVSCA/Adjacency_matrix_old.xlsx")
old.men <- read.xlsx("COVSCA/Adjacency_matrix_old_men.xlsx")
old.women <- read.xlsx("COVSCA/Adjacency_matrix_old_women.xlsx")
women <- read.xlsx("COVSCA/Adjacency_matrix_women.xlsx")
young <- read.xlsx("COVSCA/Adjacency_matrix_young.xlsx")
young.men <- read.xlsx("COVSCA/Adjacency_matrix_young_men.xlsx")
young.women <- read.xlsx("COVSCA/Adjacency_matrix_young_women.xlsx")

# Load data
data <- read.csv("../Data/Lipids_age_sex.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

# Edit the colnames to be shorter and add them as rownames

short_names <- NULL
for (name in colnames(men)) {
  split <- str_split(name, ',')
  short <- paste(split[[1]][2], split[[1]][3], sep = '&')
  short <- str_replace_all(short, '\\.', '')
  short_names <- c(short_names, short)
}
short_names

rownames(men) <- short_names
rownames(old) <- short_names
rownames(old.men) <- short_names
rownames(old.women) <- short_names
rownames(women) <- short_names
rownames(young) <- short_names
rownames(young.men) <- short_names
rownames(young.women) <- short_names

###############################################################################

# Make dissimilarity matrices with dist?
?dist

# Make a list of dissimilarity matrices (use absolute values)
diss.list <- list(abs(men), abs(old), abs(old.men), abs(old.women), abs(women), 
                  abs(young), abs(young.men), abs(young.women))

# Set a seed to get the same results
set.seed(1)

###############################################################################

# Perform INDSCAL analysis
indscal <- smacofIndDiff(diss.list, constraint = 'indscal', verbose = TRUE)

plot(indscal)
barplot(sort(indscal$sps, decreasing = TRUE), main = "Stress per Subject", cex.names = 2)
plot(indscal, plot.type = "bubbleplot")
plot(indscal, plot.type = "stressplot")
plot(indscal, plot.type = "Shepard")

###############################################################################

# Perform INDSCAL analysis with different number of dimensions
for (dimensions in 2:25) {
  indscal <- smacofIndDiff(diss.list, ndim = dimensions, constraint = 'indscal', verbose = TRUE)
  assignment.name <- paste0('indscal.', dimensions)
  assign(value = indscal, x = assignment.name)
}

plot(indscal.2, main = 'Dimensions = 2')
