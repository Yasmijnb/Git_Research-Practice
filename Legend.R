###############################################################################

# Yasmijn Balder
# 02-02-2021

# Legend

# Output
# Plot with legend

###############################################################################

# Load packages

library(ggplot2)

###############################################################################

# Load data
data <- read.csv("../Data/LipidsAgeSex_SqrtNormalization.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

mtcars$cyl <- as.factor(c('Men','Women','Young men','Old men','Young women','Old women','',''))#, levels = c('Men','Women','Young men','Old men','Young women','Old women'))

ggplot(mtcars, aes(x=wt, y=mpg, color=cyl)) + 
  geom_point() +
  # scale_x_discrete(limits=c("Men", "Women",'Young men','Old men','Young women','Old women')) + 
  scale_colour_manual(breaks = c("Men", "Women",'Young men','Old men','Young women','Old women','',''),
                      values=c('white','blue','darkblue', 'darkred','red','cyan','orange')) + 
  theme_bw() +
  theme(legend.position="bottom")
                     
  