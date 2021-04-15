###############################################################################

# Yasmijn Balder
# 19-03-2021

# Make differential connectivity plots

# Output
# GGplots for sex, men, and women, with differential connectivity per lipid

###############################################################################

# Load the packages
library(stringr)            # Used to split and shorten the lipid names
library(ggplot2)            # Used to make pretty plots

###############################################################################

# Load the data
sex.pclrc <- readRDS("Results/sex.pclrc.rds")
men.pclrc <- readRDS("Results/men.pclrc.rds")
women.pclrc <- readRDS("Results/women.pclrc.rds")

###############################################################################

# Shorten the names of the lipids
short.names <- NULL
for (long.name in names(sex.pclrc$Diff_Conn)) {
  short.name <- str_split(long.name, ', ')[[1]][2:3]
  combined <- paste(short.name, collapse = ' & ')
  short.names <- c(short.names, combined)
}

###############################################################################

# Create figure for sex

# Create a dataframe with the diff conn values
diff.conn <- as.data.frame(sex.pclrc$Diff_Conn)
# Change the name of the column
colnames(diff.conn)[1] <- 'Diff_Conn'
# Add the lipid names as a column
diff.conn$lipids <- short.names
# Create a column with significance based on the Z-score
dd <- sex.pclrc$Diff_Conn
zscore <- abs(dd-mean(dd))/sd(dd)
diff.conn$sig <- rep("Z-score \U2264 1", 21)
diff.conn$sig[which(zscore > 1)] <- 'Z-score > 1'
diff.conn$sig[which(zscore > 2)] <- 'Z-score > 2'

# Order data by lipids
diff.conn <- diff.conn[order(diff.conn$lipids, decreasing = T),]
diff.conn$lipids <- as.factor(diff.conn$lipids)
diff.conn$lipids <- factor(diff.conn$lipids, levels = rev(levels(diff.conn$lipids)))

# Create a bar plot
ggplot(data = diff.conn, aes(x = Diff_Conn, y = lipids, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("steelblue", "orange", "red")) + 
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
# Create a column with significance based on the Z-score
dd <- men.pclrc$Diff_Conn
zscore <- abs(dd-mean(dd))/sd(dd)
diff.conn$sig <- rep("Z-score < 1", 21)
diff.conn$sig[which(zscore > 1)] <- 'Z-score > 1'
diff.conn$sig[which(zscore > 2)] <- 'Z-score > 2'

# Order data by lipids
diff.conn <- diff.conn[order(diff.conn$lipids, decreasing = T),]
diff.conn$lipids <- as.factor(diff.conn$lipids)
diff.conn$lipids <- factor(diff.conn$lipids, levels = rev(levels(diff.conn$lipids)))

# Create a bar plot
ggplot(data = diff.conn, aes(x = Diff_Conn, y = lipids, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("steelblue", "orange", "red")) + 
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
# Create a column with significance based on the Z-score
dd <- women.pclrc$Diff_Conn
zscore <- abs(dd-mean(dd))/sd(dd)
diff.conn$sig <- rep("Z-score < 1", 21)
diff.conn$sig[which(zscore > 1)] <- 'Z-score > 1'
diff.conn$sig[which(zscore > 2)] <- 'Z-score > 2'

# Order data by lipids
diff.conn <- diff.conn[order(diff.conn$lipids, decreasing = T),]
diff.conn$lipids <- as.factor(diff.conn$lipids)
diff.conn$lipids <- factor(diff.conn$lipids, levels = rev(levels(diff.conn$lipids)))

# Create a bar plot
ggplot(data = diff.conn, aes(x = Diff_Conn, y = lipids, fill = sig)) + 
  geom_bar(stat = "identity") +
  # Use nicer colours
  scale_fill_manual(values = c("steelblue", "red")) + 
  # Change the axis labels
  xlab('Differential connectivity') + ylab('Lipoprotein Main Fractions') + 
  # Make plot black and white
  theme_bw(base_size = 17) +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  xlim(0, 1)

