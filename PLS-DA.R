###############################################################################

# Yasmijn Balder
# 15-02-2021

# PLS-DA analysis between sex and age using main fractions

# Output
# 

###############################################################################

# Load packages
library(caret)

###############################################################################

# Load data

# Set working directory
setwd("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/")
# Load data
data <- read.csv("Data/Lipids_age_sex.csv", check.names = FALSE)
data[,1] <- NULL

###############################################################################

sex.plsda <- plsda(data[,23:43], as.factor(data$Gender))
predict(sex.plsda)

# In the outer loop (CV2) the complete data set is split into a test set and a
# rest set: The test set is set aside, and the rest set is used in the loop CV1.
# Within CV1 the rest set is split into a validation set and a training set. The
# training set is used to develop a series of PLSDA models with different number
# of latent variables from which the samples in the validation set are
# predicted: The optimal number of components is chosen that maximizes the AUROC
# (i.e., maximize the prediction power of the model).

# Set a seed to get the same results
set.seed(1)

# Sample 
in.CV1 <- sample(seq(along = data$Gender), nrow(data)/4*3)
test <- data[in.CV1,]
rest <- data[-in.CV1,]
test.true <- data$Gender[in.CV1]
rest.true <- data$Gender[-in.CV1]

preprocess.test <- preProcess(test)
test.pred <- predict(preprocess.test, test)
rest.pred <- predict(preprocess.test, rest)

sex.plsda <- plsda(test.pred, test.true, ncomp = 5)

in.CV2
validation
training 

###############################################################################

data(mdrr)
set.seed(1)
inTrain <- sample(seq(along = mdrrClass), 450)

nzv <- nearZeroVar(mdrrDescr)
filteredDescr <- mdrrDescr[, -nzv]

training <- filteredDescr[inTrain,]
test <- filteredDescr[-inTrain,]
trainMDRR <- mdrrClass[inTrain]
testMDRR <- mdrrClass[-inTrain]

preProcValues <- preProcess(training)

trainDescr <- predict(preProcValues, training)
testDescr <- predict(preProcValues, test)

useSoftmax <- plsda(trainDescr, trainMDRR, ncomp = 5)

confusionMatrix(predict(useBayes, testDescr),
                testMDRR)

confusionMatrix(predict(useSoftmax, testDescr),
                testMDRR)

histogram(~predict(useSoftmax, testDescr, type = "prob")[,"Active",]
          | testMDRR, xlab = "Active Prob", xlim = c(-.1,1.1))
