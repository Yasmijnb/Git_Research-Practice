###############################################################################

# Yasmijn Balder
# 16-03-2021

# Make differential connectivity plot with -log10(p-value)

# Output:
# Graph per differential connectivity result

###############################################################################

# Load function from Tuuli
source('../diffcon_plot.R')

###############################################################################

# Move to results folder
setwd("Results/")
# Perform function
diffcon_plot()

###############################################################################