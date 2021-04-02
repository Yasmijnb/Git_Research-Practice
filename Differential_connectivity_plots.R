###############################################################################

# Yasmijn Balder
# 19-03-2021

# Make differential connectivity plots

# Output
# GGplots for sex, men, and women, with differential connectivity per lipid

###############################################################################

# Load the data
sex.pclrc <- readRDS("Results/sex.pclrc.rds")
men.pclrc <- readRDS("Results/men.pclrc.rds")
women.pclrc <- readRDS("Results/women.pclrc.rds")

###############################################################################

# Create figure for sex

# Create a dataframe with the diff conn values
diff.conn <- as.data.frame(sex.pclrc$Diff_Conn)
# Change the name of the column
colnames(diff.conn)[1] <- 'Diff_Conn'
# Add the lipid names as a column
diff.conn$lipids <- short.names
# Create a column with significance based on the zscore
dd <- sex.pclrc$Diff_Conn
zscore <- abs(dd-mean(dd))/sd(dd)
diff.conn$sig <- rep("Not significant", 21)
diff.conn$sig[which(zscore > 1)] <- 'Significant'

# Create a bar plot
ggplot(data = diff.conn, aes(x = Diff_Conn, y = lipids, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("orange", "steelblue")) + 
  # Change the axis labels
  xlab('Differential connectivity') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  xlim(0, 1)

###############################################################################

# Create figure for men

# Create a dataframe with the diff conn values
diff.conn <- as.data.frame(men.pclrc$Diff_Conn)
# Change the name of the column
colnames(diff.conn)[1] <- 'Diff_Conn'
# Add the lipid names as a column
diff.conn$lipids <- short.names
# Create a column with significance based on the zscore
dd <- men.pclrc$Diff_Conn
zscore <- abs(dd-mean(dd))/sd(dd)
diff.conn$sig <- rep("Not significant", 21)
diff.conn$sig[which(zscore > 1)] <- 'Significant'

# Create a bar plot
ggplot(data = diff.conn, aes(x = Diff_Conn, y = lipids, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("orange", "steelblue")) + 
  # Change the axis labels
  xlab('Differential connectivity') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  xlim(0, 1)

###############################################################################

# Create figure for women

# Create a dataframe with the diff conn values
diff.conn <- as.data.frame(women.pclrc$Diff_Conn)
# Change the name of the column
colnames(diff.conn)[1] <- 'Diff_Conn'
# Add the lipid names as a column
diff.conn$lipids <- short.names
# Create a column with significance based on the zscore
dd <- women.pclrc$Diff_Conn
zscore <- abs(dd-mean(dd))/sd(dd)
diff.conn$sig <- rep("Not significant", 21)
diff.conn$sig[which(zscore > 1)] <- 'Significant'

# Create a bar plot
ggplot(data = diff.conn, aes(x = Diff_Conn, y = lipids, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("orange", "steelblue")) + 
  # Change the axis labels
  xlab('Differential connectivity') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  xlim(0, 1)
