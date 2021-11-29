#test for proteinStructure boost
install.packages('devtools')
devtools::install_github('ShaoxunLiu/proteinStructureBoost')
library(proteinStructureBoost)

#generate adaboost model
prion <- as.character(read.csv("E:/Lin Lab/prion_r/rcsb_pdb_ids_20200623212016.txt", header=FALSE, stringsAsFactors = FALSE)[1,])
non_prion <- as.character(read.csv("E:/Lin Lab/prion_r/rcsb_pdb_ids_20200624060710.txt", header=FALSE, stringsAsFactors = FALSE)[1,])
p_pdbs <- prion[]
n_pdbs <- non_prion[c(1:500)]

fM <- featureMatrix(positive = p_pdbs, negative = n_pdbs, 
                     projectFile = 'E:/Desktop/proteinStructureBoost/test/PDBs', size = 32, as.file = FALSE, balance = TRUE, rotations = 5)
model <- adaboostTrain(fM, mfinal.sequence = c(1:10)*10)
adaboostModel <- model$model

pred <- predict(adaboostModel, attr(model, 'Test'))
pred$prob
####test LIME ####
p53Test <- featureMatrix(positive = '1TSR',  projectFile = 'E:/Desktop/proteinStructureBoost/test/PDBs', size = 32, 
                         as.file = FALSE, balance =  FALSE, rotation = 5, combine.chains = FALSE )
dna <- c(grep('_E', row.names(p53Test)),grep('_F', row.names(p53Test)))
p53Test <- p53Test[-dna,]
pred <- predict(model, p53Test)
mean(pred$prob[,2])
mean(as.numeric(pred$class))
