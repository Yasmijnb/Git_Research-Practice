###############################################################################

# Yasmijn Balder
# 20-01-2021

# Combine the data sets and save separately
# NOTE: The NAs still need to be taken care of in the complete data file

###############################################################################

# Load packages

library(stringr)    # Used to split the row names
library(readxl)     # Used to read excel files

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")

# Load IVDr_lipo (NMR estimation of lipidomes)
lipids <- read.csv("Data/IVDr_lipo.txt", sep = ' ', check.names = FALSE)
# Change the rownames to match (stringr package needed)
for (row in 1:nrow(lipids)) {
  new.name <- str_split(rownames(lipids)[row], '/')[[1]][2]
  rownames(lipids)[row] <- new.name
}
# Load AVIS_2020
clinic <- as.data.frame(read_excel("Data/AVIS_2020.xlsx", sheet = 3))
# Change certain columns to numeric
for (column in c(2:14,18:42)) {
  clinic[,column] <- as.numeric(clinic[,column])
}
# Change a colname
colnames(clinic)[which(colnames(clinic) == '...1')] <- 'ID'

###############################################################################

# Make the first data set: estimated lipids with age and sex

# Merge the data sets into one
data.age.sex <- merge(clinic[,c(1, 2, 13)], lipids, by.x = 'ID', 
                      by.y = 'row.names')
# Save the data set
write.csv(data.age.sex, file = "Data/Lipids_age_sex.csv")

###############################################################################

# Make the second data set: estimated lipids with all the other variables

# NOTE: The NAs still need to be taken care of

# Merge the data sets into one
complete.data <- merge(clinic, lipids, by.x = 'ID', by.y = 'row.names')
# Save the data set
write.csv(complete.data, file = "Data/Lipids_all.csv")

###############################################################################

# Take care of the NAs (which might all be zeros)

# # How many NAs are there in total?
# length(which(is.na(complete.data)))
# 
# # See how many NAs each variable has
# NAs <- NULL
# for (column in 1:ncol(complete.data)) {
#   number <- length(which(is.na(complete.data[,column])))
#   perc <- number / nrow(complete.data) * 100
#   if (perc != 0) {
#     NAs <- rbind(NAs, c(colnames(complete.data)[column], round(perc,2)))
#   }
# }
# 
# # Turn all NAs into a zero
# complete.data[which(is.na(complete.data), arr.ind = T)] <- 0

