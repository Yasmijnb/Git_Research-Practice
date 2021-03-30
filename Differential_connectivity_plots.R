###############################################################################

# Yasmijn Balder
# 19-03-2021

# Make differential connectivity plots

# Output
# GGplots for sex, men, and women, with differential connectivity per lipid

###############################################################################

# Load the data
sex.pclrc <- readRDS("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/sex.pclrc.rds")
men.pclrc <- readRDS("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/men.pclrc.rds")
women.pclrc <- readRDS("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/women.pclrc.rds")

###############################################################################

# Create figure for sex

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
  scale_fill_manual(values = "steelblue") + 
  # Change the axis labels
  xlab('Differential connectivity') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw() +
  # Remove the legend title
  theme(legend.title = element_blank())

###############################################################################

# Create figure for men

# Create a dataframe with the adjusted p-values
adjusted.pvalues <- as.data.frame(men.pclrc$Pval_adj)
# Change the name of the p-values column
colnames(adjusted.pvalues)[1] <- 'P.values'
# Add the lipid names as a column
adjusted.pvalues$lipids <- short.names
# Add the differential connectivity as a column
adjusted.pvalues$diffcon <- men.pclrc$Diff_Conn
# Create a column with significance based on the adjusted p-values
adjusted.pvalues$sig <- rep("Not significant", 21)
adjusted.pvalues$sig[which(adjusted.pvalues$P.values <= 0.05)] <- 'Significant'

# Create a bar plot
ggplot(data = adjusted.pvalues, aes(x = diffcon, y = lipids, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("orange", "steelblue")) + 
  # Change the axis labels
  xlab('Differential connectivity') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw() +
  # Remove the legend title
  theme(legend.title = element_blank())

###############################################################################

# Create figure for women

# Create a dataframe with the adjusted p-values
adjusted.pvalues <- as.data.frame(women.pclrc$Pval_adj)
# Change the name of the p-values column
colnames(adjusted.pvalues)[1] <- 'P.values'
# Add the lipid names as a column
adjusted.pvalues$lipids <- short.names
# Add the differential connectivity as a column
adjusted.pvalues$diffcon <- women.pclrc$Diff_Conn
# Create a column with significance based on the adjusted p-values
adjusted.pvalues$sig <- rep("Not significant", 21)
adjusted.pvalues$sig[which(adjusted.pvalues$P.values <= 0.05)] <- 'Significant'

# Create a bar plot
ggplot(data = adjusted.pvalues, aes(x = diffcon, y = lipids, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("orange", "steelblue")) + 
  # Change the axis labels
  xlab('Differential connectivity') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw() +
  # Remove the legend title
  theme(legend.title = element_blank())

