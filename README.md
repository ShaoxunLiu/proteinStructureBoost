# proteinStructureBoost
Package for boosting machine learning models for function/family recognition using protein structures. 

# Usage
install.packages('devtools')
devtools::install_github('ShaoxunLiu/proteinStructureBoost')
library(proteinStructureBoost)

# Inputs to prepare
List of pdb ids for positive and nagative labels. 

# Procudure
fM <- featureMatrix(postive, negative)
model <- adaboostTrain(fM)
