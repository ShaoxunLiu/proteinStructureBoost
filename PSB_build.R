#package creation
require("devtools")
require("roxygen2")
setwd('E:/Desktop/proteinStructureBoost')
devtools::create("proteinStructureBoost")

#update documentation
setwd('E:/Desktop/proteinStructureBoost/proteinStructureBoost')
devtools::document()

  
#update and test package
detach("package:proteinStructureBoost",unload = TRUE)
devtools::install('E:/Desktop/proteinStructureBoost/proteinStructureBoost')
 
library(proteinStructureBoost)


