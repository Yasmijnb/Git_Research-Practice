###############################################################################

# Yasmijn Balder
# 15-02-2021

# PLS-DA analysis between men and women using lipoprotein main fractions

# Output
# Accuracy for each number of components

###############################################################################

# Load packages
library(caret)        # Used for PLS-DA
library(caTools)      # Used for AUC

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/Lipids_age_sex.csv")
data[,1] <- NULL

###############################################################################

## Trial

# data$Gender <- as.factor(data$Gender)
# 
# set.seed(1)
# inTrain <- sample(seq(along = data$Gender), 450)
# 
# training <- data[inTrain,23:43]
# test <- data[-inTrain,23:43]
# trainGender <- data$Gender[inTrain]
# testGender <- data$Gender[-inTrain]
# 
# sex.plsda <- plsda(training, as.factor(trainGender))
# 
# confusionMatrix(predict(sex.plsda, test), testGender)
# 
# # AUC
# mean(colAUC(test, testGender))
# # Sensitivity
# confusionMatrix(predict(sex.plsda, test), testGender)$byClass["Sensitivity"]
# # Specificity
# confusionMatrix(predict(sex.plsda, test), testGender)$byClass["Specificity"]
# 
# histogram(~predict(sex.plsda, test, type = "prob")[,"man",]
#           | testGender, xlab = "woman", xlim = c(-.1,1.1),
#           main = "Certainty of predictions")

###############################################################################

# In the outer loop (CV2) the complete data set is split into a test set and a
# rest set: The test set is set aside, and the rest set is used in the loop CV1.
# Within CV1 the rest set is split into a validation set and a training set. The
# training set is used to develop a series of PLSDA models with different number
# of latent variables from which the samples in the validation set are
# predicted: The optimal number of components is chosen that maximizes the AUROC
# (i.e., maximize the prediction power of the model).

# Make sure gender is seen as a factor
data$Gender <- as.factor(data$Gender)

# Set a seed to get the same results
set.seed(1)

# Outer loop (CV 2) for overall accuracy
in.CV2 <- sample(seq(along = data$Gender), nrow(data)/4)
test <- data[in.CV2,23:43]
rest <- data[-in.CV2,23:43]
test.true <- data$Gender[in.CV2]
rest.true <- data$Gender[-in.CV2]

# Inner loop (CV 1) for optimal number of components
in.CV1 <- sample(seq(along = in.CV2), nrow(rest)/3)
validation <- data[in.CV1,23:43]
training <- data[-in.CV1,23:43]
validation.true <- data$Gender[in.CV1]
training.true <- data$Gender[-in.CV1]

accuracy.comp <- NULL
for (comp in 1:21) {
  sex.plsda.comp <- plsda(training, training.true, ncomp = comp)
  predict.values <- confusionMatrix(predict(sex.plsda.comp, validation), 
                                    validation.true)
  accuracy <- predict.values$overall["Accuracy"]
  accuracy.comp <- c(accuracy.comp, accuracy)
}
plot(accuracy.comp, xlab = 'Number of components', ylab = 'Accuracy', pch = 19,
     col = 'white')
text(accuracy.comp, labels = 1:21)

# The choice was made for 5 components
sex.plsda <- plsda(test, test.true, ncomp = 5)


# Now try with resampling
n.resamp <- 100
for (comp in 1:21) {
  accuracy.comp <- NULL
  for (i in 1:n.resamp) {
    # Resample
    in.CV1 <- sample(seq(along = in.CV2), nrow(rest)/3)
    validation <- data[in.CV1,23:43]
    training <- data[-in.CV1,23:43]
    validation.true <- data$Gender[in.CV1]
    training.true <- data$Gender[-in.CV1]
    
    # Make a model
    sex.plsda.comp <- plsda(training, training.true, ncomp = comp, 
                            validation = 'CV')
    predict.values <- confusionMatrix(predict(sex.plsda.comp, validation), 
                                      validation.true)
    accuracy <- predict.values$overall["Accuracy"]
    accuracy.comp <- c(accuracy.comp, accuracy)
  }
  print(paste('Accuracy for', comp, 'components:', round(mean(accuracy.comp), 3)))
}

# Take care of unbalance still
# nsize = round(min(n1,n2) * perc.unbal)

max.perm = 1000
