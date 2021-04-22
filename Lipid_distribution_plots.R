###############################################################################

# Yasmijn Balder
# 21-04-2021

# Lipid distribution

# Output
# Plots showing the distribution of the lipidome

###############################################################################

# Load packages

library(ggplot2)
library(scales)
library(dplyr)

###############################################################################

# Load data
data <- read.csv("../Data/LipidsAgeSex_NoNormalization.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

dist <- colMeans(data[,23:43])
dist <- as.data.frame(dist)
dist$lipid <- c(rep('Triglyceride', 4), rep('Cholesterol', 4), 
                rep('Free cholesterol', 4), rep('Phospholipid', 4), rep('Apo', 5))
dist$lipo <- c('VLDL', 'IDL', 'LDL', 'HDL', 'VLDL', 'IDL', 'LDL', 'HDL', 'VLDL',
               'IDL', 'LDL', 'HDL', 'VLDL', 'IDL', 'LDL', 'HDL', 'HDL', 'HDL', 
               'VLDL', 'IDL', 'LDL')
dist <- dist[order(dist$lipo),]

###############################################################################

lipid_dist <- data.frame(sum(dist$dist[which(dist$lipid == 'Apo')]),
                        sum(dist$dist[which(dist$lipid == 'Cholesterol')]),
                        sum(dist$dist[which(dist$lipid == 'Free cholesterol')]),
                        sum(dist$dist[which(dist$lipid == 'Triglyceride')]),
                        sum(dist$dist[which(dist$lipid == 'Phospholipid')]))
lipid_dist <- rbind(lipid_dist, c('Apo', 'Cholesterol', 'Free cholesterol', 
                                  'Triglyceride','Phospholipid'))
lipid_dist <- t(lipid_dist)
colnames(lipid_dist) <- c('dist','lipid')
lipid_dist <- as.data.frame(lipid_dist)
lipid_dist$dist <- as.numeric(lipid_dist$dist)

# Make percentage labels
lipid_percent <- lipid_dist %>%
  mutate(lipid = factor(lipid, 
                       levels = rev(unique(lipid))),
         cumulative = cumsum(dist),
         midpoint = cumulative - dist / 2,
         labels = paste0(round((dist/ sum(dist)) * 100, 1), "%"))

ggplot(lipid_percent, aes(x = "", y = dist, fill = lipid)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = rev(c("#0073C2", "#EFC000", "#868686", "#7AA6DC", "#CD534C"))) + 
  theme_void() +
  geom_text(aes(x = 1.3, y = midpoint , label = labels), color="black",
            fontface = "bold") + 
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())

###############################################################################

lipo_dist <- data.frame(sum(dist$dist[which(dist$lipo == 'HDL')]),
                        sum(dist$dist[which(dist$lipo == 'IDL')]),
                        sum(dist$dist[which(dist$lipo == 'LDL')]),
                        sum(dist$dist[which(dist$lipo == 'VLDL')]))
lipo_dist <- rbind(lipo_dist, c('HDL', 'IDL', 'LDL', 'VLDL'))
lipo_dist <- t(lipo_dist)
colnames(lipo_dist) <- c('dist','lipo')
lipo_dist <- as.data.frame(lipo_dist)
lipo_dist$dist <- as.numeric(lipo_dist$dist)

# Make percentage labels
lipo_percent <- lipo_dist %>%
  mutate(lipo = factor(lipo, 
                       levels = rev(unique(lipo))),
         cumulative = cumsum(dist),
         midpoint = cumulative - dist / 2,
         labels = paste0(round((dist/ sum(dist)) * 100, 1), "%"))

ggplot(lipo_percent, aes(x = "", y = dist, fill = lipo)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "GnBu", direction = 1) + 
  theme_void() +
  geom_text(aes(x = 1.3, y = midpoint , label = labels), color="black",
            fontface = "bold") + 
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_blank())

