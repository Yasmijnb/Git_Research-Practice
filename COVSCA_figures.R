###############################################################################

# Yasmijn Balder
# 19-03-2021

# Make figures of COVSCA results

# Output
# Three figures. 
  # 1. Scores plot of COVSCA
  # 2. Loadings plot of first COVSCA prototype
  # 3. Loadings plot of second COVSCA prototype

###############################################################################

# Load the packages
library(stringr)            # Used to split and shorten the lipid names
library(ggplot2)            # Used to make pretty plots

###############################################################################

# Load the data
scores <- read.csv('COVSCA/COVSCA_scores.csv', header = FALSE)
loadings <- read.csv('COVSCA/COVSCA_loadings.csv', header = FALSE)

###############################################################################

# Change the colnames of the scores
colnames(scores) <- c('1st COVSCA prototype', '2nd COVSCA prototype')
scores$age <- c('all','all','young','old','young','old')
scores$gender <- c('man', 'woman', 'man', 'woman', 'woman', 'woman')
scores$age_gender <- c('men','women','young men','old men','young women','old women')

# COVSCA scores plot
ggplot(scores, aes(x = `1st COVSCA prototype`, y = `2nd COVSCA prototype`, 
                   color = age_gender)) + 
  geom_point(size = 3) +
  # Use manual colours
  scale_color_manual(values = c('blue','darkblue','darkred','red','cyan','orange')) + 
  # Add labels of age
  # geom_text(label=scores$age) +
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank())

###############################################################################

# Use short names for the loadings plots
loadings$X <- c("Triglycerides VLDL", "Triglycerides IDL", 
                        "Triglycerides LDL", "Triglycerides HDL", 
                        "Cholesterol VLDL", "Cholesterol IDL", 
                        "Cholesterol LDL", "Cholesterol HDL",
                        "FreeCholesterol VLDL", "FreeCholesterol IDL",
                        "FreeCholesterol LDL", "FreeCholesterol HDL",
                        "Phospholipids VLDL", "Phospholipids IDL",
                        "Phospholipids LDL", "Phospholipids HDL",
                        "ApoA1 HDL", "ApoA2 HDL", "ApoB VLDL", 
                        "ApoB IDL", "ApoB LDL")

# Split the loadings and change format for gg plot
first_loadings <- data.frame(c(loadings$X, loadings$X), 
                             c(rep('First', 21), rep('Second', 21)), 
                             c(loadings$V1, loadings$V2))
colnames(first_loadings) <- c('Lipids','Loadings','Value')
second_loadings <- data.frame(c(loadings$X, loadings$X), 
                             c(rep('First', 21), rep('Second', 21)), 
                             c(loadings$V3, loadings$V4))
colnames(second_loadings) <- c('Lipids','Loadings','Value')

# COVSCA first loadings plot
ggplot(first_loadings, aes(x = Lipids, y = Value, fill = Loadings)) + 
  geom_bar(position = "dodge", stat="identity") + 
  # Change label names
  labs(x = "Lipoprotein Main Fractions", y = "Loadings") +
  # Use costum colours
  scale_fill_manual(values = c('orange','steelblue')) +
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Change x axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Remove the legend title
  theme(legend.title = element_blank())

# COVSCA second loadings plot
ggplot(second_loadings, aes(x = Lipids, y = Value, fill = Loadings)) + 
  geom_bar(position = "dodge", stat="identity") + 
  # Change label names
  labs(x = "Lipoprotein Main Fractions", y = "Loadings") +
  # Use costum colours
  scale_fill_manual(values = c('orange','steelblue')) +
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Change x axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Remove the legend title
  theme(legend.title = element_blank())

