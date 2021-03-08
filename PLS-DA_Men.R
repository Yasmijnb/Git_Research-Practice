###############################################################################

# Yasmijn Balder
# 15-02-2021

# PLS-DA analysis between young and old men using lipoprotein main fractions

# Output
# AUC for each number of components

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/Lipids_age_sex.csv", check.names = TRUE)
data[,1] <- NULL

###############################################################################

# Make sure gender is seen as a factor
data$Gender <- as.factor(data$Gender)
# Use only men
men <- data[which(data$Gender == 'man'),]

# Make young and old groups
men.age_groups <- rep(0, nrow(data))
men.age_groups[which(men$Age < 32)] <- 'young'
men.age_groups[which(men$Age > 49)] <- 'old'

# Use only young and old, not the middle age group
men.age.data <- data[which(men.age_groups != 0),]
men.age_groups <- as.factor(men.age_groups[which(men.age_groups != 0)])

###############################################################################

# Make outer loop (CV 2) 
in.CV2 <- sample(seq(along = men.age.data$Gender), nrow(men.age.data)/4)
test <- men.age.data[in.CV2,23:43]
rest <- men.age.data[-in.CV2,23:43]
test.true <- men.age_groups[in.CV2]
rest.true <- men.age_groups[-in.CV2]

# Inner loop for deciding the number of components
n.resamp <- 100
for (comp in 1:21) {
  accuracy.comp <- NULL
  auc.comp <- NULL
  for (i in 1:n.resamp) {
    # Resample
    in.CV1 <- sample(seq(along = in.CV2), nrow(rest)/3)
    validation <- men.age.data[in.CV1,23:43]
    training <- men.age.data[-in.CV1,23:43]
    validation.true <- men.age_groups[in.CV1]
    training.true <- men.age_groups[-in.CV1]
    
    # Make a model
    men.plsda.comp <- plsda(training, training.true, ncomp = comp, 
                            validation = 'CV')
    # Make predictions
    predict.values <- confusionMatrix(predict(men.plsda.comp, validation), 
                                      validation.true)
    # Get accuracy
    accuracy <- predict.values$overall["Accuracy"]
    accuracy.comp <- c(accuracy.comp, accuracy)
    # Get AUC
    auc <- colAUC(as.numeric(predict(men.plsda.comp, validation)), 
           as.numeric(validation.true))
    auc.comp <- c(auc.comp, auc)
  }
  # print(paste('Accuracy for', comp, 'components:', round(mean(accuracy.comp), 3)))
  print(paste('AUC for', comp, 'components:', round(mean(auc.comp), 3)))
}

# # Use 3 components (Edoardo)
# comp <- 3
# 
# # Outer loop 
# max.perm = 1000
# accuracy.perm <- NULL
# auc.perm <- NULL
# for (perm in 1:max.perm) {
#   # Make a model
#   men.plsda.perm <- plsda(training, training.true, ncomp = comp, 
#                           validation = 'CV')
#   # Make predictions
#   predict.values <- confusionMatrix(predict(men.plsda.perm, validation), 
#                                     validation.true)
#   # Get accuracy
#   accuracy <- predict.values$overall["Accuracy"]
#   accuracy.perm <- c(accuracy.perm, accuracy)
#   # Get AUC
#   auc <- colAUC(as.numeric(predict(men.plsda.perm, validation)), 
#                 as.numeric(validation.true))
#   auc.perm <- c(auc.perm, auc)
# }
# 
# mean(accuracy.perm)
# mean(auc.perm)

###############################################################################
