#' Download vector of PDB files
#'
#' This function downloads a vector of PDB files into the given folder, and
#' returns a vector of the directories of the PDB files.
#'
#' @param pdbids Vector of pdb ids to be downloaded.
#' @param projectFile The directory in which the PDB files are downloaded.
#' @param as.file If true, creates a /PDBs folder to store downloaded files.
#' @return A vector of the directories of the downloaded PDB files.
#' @export
download.pdb <- function(pdbids,
                         projectFile = '~/',
                         as.file = TRUE){
  if(as.file){
    pdb_file <- paste(projectFile,'/PDBs/', sep = '')
    if(!dir.exists(pdb_file)){
      dir.create(pdb_file)
    }
  }else{
    pdb_file <- paste(projectFile, '/',sep = '')
  }

  options(scipen = 999)
  out <- c()
  for(i in 1:length(pdbids)){
    pdbid <- pdbids[i]
    new <- length(grep(pdbid, list.files(pdb_file))) == 0
    dir <- paste(pdb_file, pdbid,".pdb", sep = "")
    if(new){
      tryCatch({
        url <- paste("https://files.rcsb.org/download/", pdbid ,".pdb", sep = "")
        download.file(url, destfile = dir)
      },error = function(e) e)
    }
    if(file.exists(dir)){
      out <- c(out, dir)
    }else{
      cat(pdbid, 'not successfully downloaded \n')
    }

  }
  options(scipen = 0)
  return(out)
}


#'Default atom.value
#'
#'@description data.frame for assigning numeric values for atoms in form.grid
#'function. row.names are fixed. col.name = 'value', row.names = atom atom (H, N
#', C, O, S, CA and CB for alpha and beta carbon), value = assigned value for
#'each atom (0~1).
#'@export
atom.value <- data.frame(value = c(0.1, 0.3, 0.5, 0.7, 0.75, 0.9, 1),
                         row.names = c('H', 'N', 'C', 'O', 'S', 'CB', 'CA'))


#'Default atom.color
#'
#'@description data.frame for assigning color codes for atoms in show.grid
#'function. row.names are fixed. col.name = 'value', row.names = atom (H, N, C,
#'O, S, CA and CB for alpha and beta carbon), value = string for valid colors
#'(i.e. 'red').
#'@export
atom.color <- data.frame(color = c('yellow', 'red', 'green', 'blue', 'orange', 'green', 'green'),
                         row.names = c('H', 'N', 'C', 'O', 'S', 'CB', 'CA'))



#'Form 3D grid for pdb file
#'
#'This function generates a 3D grid for a given pdb structure or pdb file
#'directory, with each atom shown as the designated atom.value.
#'
#'@param pdb pdb structure from bio3d::read.pdb. (if entered, should not enter
#'pdb.file).
#'@param pdb.name Name of the pdb structure entered. Only enter when pdb
#'argument is used, and pdb.file is not.
#'@param pdb.file Directory to the pdb file. (if entered, should not enter pdb).
#'@param combine.chains If true, assign all chains in the pdb file in one grid.
#'If false, generates a separate grid for each chain.
#'@param size The length of the recognition cube in Angstroms
#'@param units The number of blocks in each axis.
#'@param center.median If true, centers each chain (or all chains depending on
#'combine.chains) to the median of their xyz.
#'@param center Designated center for all chains in pdb. Only enter value when
#'center.median = FALSE.
#'@param show.grid If true, outputs an .html file that contains the 3D
#'representation of the grid. Directory to output 3D representation .html file
#'is the same as the .pdb file.
#'@return A list of 3d arrays containing the 3D grids of each chain (or all
#'chains as a whole when combine.chains = TRUE) of the pdb file with attributes
#'of input parameters.
#'@export
form.grid <- function(pdb = NA, pdb.name = 'pdb0', pdb.file = NA,
                      combine.chains = TRUE,
                      size = 32,
                      units = 16,
                      center.median = TRUE,
                      center = NA,
                      show.grid = TRUE
                      ){
  if(is.na(pdb[1]) && is.na(pdb.file)){
    cat('Provide one of pdb for pdb.file. \n')
    return()
  }else if(is.na(pdb[1])){
    pdb <- bio3d::read.pdb(pdb.file)
    placehold <- FALSE
  }else if (is.na(pdb.file)){
    placehold <- TRUE
    pdb <- pdb
  }else{
    cat('pdb arugument is used over pdb.file \n')
    pdb <- pdb
    placehold <- TRUE
  }
  pdb_info <- pdb$atom
  if(combine.chains){
    pdb_info$chain <- 'A'
  }
  chains <- unique(pdb_info$chain)
  all_chains <- pdb_info
  out <- list()
  for(i in 1:length(chains)){
    pdb_info <- all_chains[all_chains$chain == chains[i],]
    if(center.median){
      x.center <- median(pdb_info$x)
      y.center <- median(pdb_info$y)
      z.center <- median(pdb_info$z)
      pdb_info$x <- pdb_info$x - x.center
      pdb_info$y <- pdb_info$y - y.center
      pdb_info$z <- pdb_info$z - z.center
    }else if(length(center) == 3){
      x.center <- center[1]
      y.center <- center[2]
      z.center <- center[3]
      pdb_info$x <- pdb_info$x - x.center
      pdb_info$y <- pdb_info$y - y.center
      pdb_info$z <- pdb_info$z - z.center
    }else{
      cat('Center argumet not correctly entered. Should be list of 3. \n')
    }
    #trim pdb according to size
    pdb_info$x <- pdb_info$x + size/2
    pdb_info$y <- pdb_info$y + size/2
    pdb_info$z <- pdb_info$z + size/2
    pdb_info <- pdb_info[pdb_info$x > 0,]
    pdb_info <- pdb_info[pdb_info$x < size,]
    pdb_info <- pdb_info[pdb_info$y > 0,]
    pdb_info <- pdb_info[pdb_info$y < size,]
    pdb_info <- pdb_info[pdb_info$z > 0,]
    pdb_info <- pdb_info[pdb_info$z < size,]
    #change xyz to coordinates
    pdb_info$x <- floor(pdb_info$x*units/size)
    pdb_info$y <- floor(pdb_info$y*units/size)
    pdb_info$z <- floor(pdb_info$z*units/size)
    #extract data.frame and add value
    pdb_info <- pdb_info[,c('elety', 'chain', 'x','y','z')]
    if(nrow(pdb_info) == 0){
      cat('chain', chains[i], 'is empty \n')
      next()
    }
    value <- c()
    for(j in 1:nrow(pdb_info)){
      add <- 0
      atom <- pdb_info$elety[j]
      if(atom != 'CA' && atom != 'CB'){
        atom <- substr(atom,1,1)
      }
      value <- c(value, atom.value[atom, 'value'])
    }
    pdb_info <- cbind.data.frame(pdb_info, value = value)
    #trim pdb_info
    pdb_info <- dplyr::arrange(pdb_info, desc(pdb_info$value))
    pdb_info <- pdb_info[!duplicated(pdb_info[,c('x','y','z')]),]


    #fit pdb into grid
    grid <- array(dim = c(units, units, units), data = 0)
    for(j in 1:nrow(pdb_info)){
      pdb_line <- pdb_info[j,]
      grid[pdb_line$x + 1, pdb_line$y + 1, pdb_line$z + 1] <- pdb_line$value
    }
    if(placehold){
      attr(grid, 'file') <- paste('placeHolder/', pdb.name, '.pdb', sep = '')
    }else{
      attr(grid, 'file') <- pdb.file
    }
    if(combine.chains){
      ch <- 'all'
      attr(grid, 'chain') <- ch
    }else{
      ch <- chains[i]
      attr(grid, 'chain') <- ch
    }
    attr(grid, 'size') <- size
    attr(grid, 'units') <- units
    attr(grid, 'center') <- c(x.center, y.center, z.center)
    grid[is.na(grid)] <- 0
    out[[i]] <- grid
    #call show.grid
    if(show.grid){
      #show the grid by calling grid3D
      plot_file <- paste(pdb.file,'.', ch,'.html', sep = '')
      plot3D <- grid3D(grid)
      htmlwidgets::saveWidget(plot3D, plot_file)
      cat('3D plot for ', ch, ' chain(s) saved in ', plot_file, '\n',sep = '')
    }
  }
  return(out)
}


#'3D visualization for pdb grid
#'
#'This function generates a 3D visualization of a 3D grid generated from a
#'form.grid function.
#'
#'@param grid Out put of form.grid function. A 3D grid of pdb information.
#'@param show.cude If true, shows the enclosing cube of recognition area.
#'@param show.color If true, uses different colors to represent different atoms.
#'Uses the data set atom.color as index. More information in ?atom.color.
#'@return A plotly 3D scatter plot to show the information in a 3D grid.
#'@export
grid3D <- function(grid,
                   show.cube = TRUE,
                   show.color = TRUE){
  table_grid <- as.data.frame(data.table::as.data.table(grid))
  colnames(table_grid) <- c('x', 'y', 'z', 'value')
  units <- dim(grid)[1]
  l <- units + 1
  cube <- matrix(ncol = 3, nrow = 16, byrow = TRUE,
                 data = c(0,0,0,
                          l,0,0,
                          l,l,0,
                          l,l,l,
                          l,0,l,
                          0,0,l,
                          0,l,l,
                          0,l,0,
                          0,0,0,
                          0,0,l,
                          0,l,l,
                          l,l,l,
                          l,0,l,
                          l,0,0,
                          l,l,0,
                          0,l,0))

  axx <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  color.value <- cbind.data.frame(atom.value, atom.color)
  color.value <- color.value[!duplicated(color.value$value),]

  #form grid
  grid3D <- plotly::plot_ly()
  grid3D <- plotly::layout(p = grid3D, scene = list(xaxis = axx, yaxis = axx, zaxis = axx))
  if(show.cube){
    grid3D <- plotly::add_trace(p = grid3D, data = cube, x = cube[,1],
                                y=cube[,2],z=cube[,3],color = I('black'),
                                type = 'scatter3d',
                                mode = 'lines',
                                line = list(width = 1),
                                opacity = 0.5,
                                name = 'cube')
  }
  if(!show.color){
    color.value$color <- 'green'
  }
  for(i in 1:nrow(color.value)){
    table_grid_v <- table_grid[table_grid$value == color.value$value[i],]
    if(nrow(table_grid_v) < 1){
      next()
    }
    grid3D <- plotly::add_trace(p = grid3D, table_grid_v, x = as.numeric(table_grid_v$x), y = as.numeric(table_grid_v$y),
                                z = as.numeric(table_grid_v$z), color = I(color.value$color[i]),
                                type = 'scatter3d',
                                mode = 'markers',
                                marker = list(size = 10, opacity = 0.8*as.numeric(color.value$value[i])),
                                name = row.names(color.value)[i])
  }
  return(grid3D)
}



#'3D comparison coordinates for 3D-LBP
#'
#'@description Functional matrix for 3D-LBP algorithm. Shall not alter.
#'@export
coord <- matrix(ncol = 3, nrow = 27, data = 0)
for(j in 1:27){
  i <- j - 1
  x <- floor(i/9)
  coord[j,1] <- x
  y <- floor((i-9*x)/3 )
  coord[j,2] <- y
  coord[j,3] <- i - 9*x - 3*y
}
coord <- coord[-14,]


#'Form LBP matrix from 3D grid
#'
#'Generate the LBP matrix from a 3D grid generated from from grid.
#'
#'@param grid Out put of form.grid function. A 3D grid of pdb information.
#'@return LBP matrix 2^26 possible levels
#'@export
lbp3d <- function(grid){
  dim <- dim(grid)[1] - 2
  lbp <- array(0, dim = c(dim,dim,dim))
  for(i in 1:dim){
    for(j in 1: dim){
      for(k in 1:dim){
        for(co in 1:nrow(coord)){
          if(grid[i+1,j+1,k+1] >= grid[i+coord[co,1],j+coord[co,2],k+coord[co,3]]){
            lbp[i,j,k] <- lbp[i,j,k] + 2^(co-1)
          }
        }
      }
    }
  }
  attr(lbp,'name') <- gsub('.pdb', '',strsplit(attr(grid,'file'), '/')[[1]][length(strsplit(attr(grid,'file'), '/')[[1]])])
  attr(lbp, 'chain') <- attr(grid, 'chain')
  return(lbp)
}


#' Form a row in feature matrix from lbp
#'
#' This function transforms a single lbp matrix into 1 line in the feature
#' matrix
#'
#' @param lbp lbp matrix generated by lbp3D function
#' @param n.feature Suggested number of features in the feature matrix. Will be
#' fixed to multiples of 68 (a result of lbp3D algorithm).
#' @param label Label of the sample features. Should be entered as 0 or 1 for
#' non-factorized models.
#' @param name Name of the feature line. Will be used as rowname in feature
#' matrix. Default is file name.
#' @return A line of lbp features of the given lbp matrix, with the label as the
#'  last element. Also shows the lbp histogram.
#' @export
lbp2Feature <- function(lbp, n.feature = 68, label = 1, name = attr(lbp, 'name')){
  pattern <- c(1,lbp,67108864)
  h <- graphics::hist(pattern, breaks=n.feature, ylim = c(0,100))
  feature <- h$counts
  label <- label
  feature.line <- matrix(nrow = 1, ncol = length(feature)+1, data = c(feature, label))
  feature.line <- as.data.frame(feature.line)
  colnames(feature.line) <- c(paste('F', c(1:length(feature)), sep = ''), 'Label')
  rownames(feature.line) <- paste(name, '_', attr(lbp, 'chain'), sep = '')
  return(feature.line)
}


#' Rotate a pdb structure
#'
#' This function takes a pdb structure read by bio3d::read.pdb function and
#' rotates it.
#'
#' @param pdb pdb structure output from bio3d::read.pdb
#' @param degree A vector of 3 numeric elements indicating the rotation degree
#' along x, y, and z axis.
#' @return A pdb structure with $Atom$x(y,z) values changed
#' @export
rotatePDB <- function(pdb, degree = c(0,0,0)){
  for(i in 1:nrow(pdb$atom)){
    #change x
    a <- as.numeric(pdb$atom$y[i])
    b <- as.numeric(pdb$atom$z[i])
    len <- sqrt(a^2 + b^2)
    deg <- degree[1] * pi / 180
    ang <- atan2(a,b)
    newang <- ang + deg
    newa <- sin(newang) * len
    newb <- cos(newang) * len
    pdb$atom$y[i] <- newa
    pdb$atom$z[i] <- newb
    #change y
    a <- as.numeric(pdb$atom$x[i])
    b <- as.numeric(pdb$atom$z[i])
    len <- sqrt(a^2 + b^2)
    deg <- degree[2] * pi / 180
    ang <- atan2(a,b)
    newang <- ang + deg
    newa <- sin(newang) * len
    newb <- cos(newang) * len
    pdb$atom$x[i] <- newa
    pdb$atom$z[i] <- newb
    #change z
    a <- as.numeric(pdb$atom$x[i])
    b <- as.numeric(pdb$atom$y[i])
    len <- sqrt(a^2 + b^2)
    deg <- degree[3] * pi / 180
    ang <- atan2(a,b)
    newang <- ang + deg
    newa <- sin(newang) * len
    newb <- cos(newang) * len
    pdb$atom$x[i] <- newa
    pdb$atom$y[i] <- newb
  }
  return(pdb)
}



#' Form a feature matrix using lists of positive and negative pdb samples
#'
#' This function takes a vector of positively labeled pdb file names and a
#' negatively labeled pdb file names to generate a feature matrix for machine
#' learning model formation. pbd files can be randomly rotated to increase
#' sample size.
#'
#' @param positive Character vector of positively labeled pdb files.
#' @param negative Character vector of negatively labeled pdb files.
#' @param projectFile The directory in which the PDB files are downloaded.
#' @param as.file If true, creates a /PDBs folder to store downloaded files.
#' @param balance If true balances the number of positive and negative samples.
#' @param positive.label Label assigned for positive samples.
#' @param negative.label Label assigned for negative samples.
#' @param rotations Number of rotational augmentation samples added to the
#' matrix.
#' @param size The length of the recognition cube in Angstroms
#' @param units The number of blocks in each axis.
#' @param combine.chains If true, assign all chains in the pdb file in one grid.
#'If false, generates a separate grid for each chain.
#' @param center.median If true, centers each chain (or all chains depending on
#'combine.chains) to the median of their xyz.
#' @param show.grid If true, outputs an .html file that contains the 3D
#'representation of the grid. Directory to output 3D representation .html file
#'is the same as the .pdb file.
#' @param n.feature Suggested number of features in the feature matrix. Will be
#' fixed to multiples of 68 (a result of lbp3D algorithm).
#' @return A data.frame containing the feature matrix of the input samples
#' @export
featureMatrix <- function(
  positive,
  negative,
  projectFile = '~/',
  as.file = TRUE,
  balance = FALSE,
  positive.label = 1,
  negative.label = 0,
  rotations = 0,
  size = 32,
  units = 16,
  combine.chains = TRUE,
  show.grid = FALSE,
  n.feature = 68){
  #function
  if(balance){
    np <- length(positive)
    nn <- length(negative)
    if(np > nn){
      r <- sample(c(1:np), nn)
      positive <- positive[r]
      np <- nn
    }else if(nn > np){
      r <- sample(c(1:nn), np)
      negative <- negative[r]
      nn <- np
    }
  }
  p_pdbs <- download.pdb(positive, projectFile, as.file)
  n_pdbs <- download.pdb(negative, projectFile, as.file)

  featureMatrix <- data.frame(NULL)
  for(i in 1:length(p_pdbs)){
    grid <- form.grid(pdb.file = p_pdbs[i], size = size, units = units, combine.chains = combine.chains,  show.grid = show.grid)
    for(j in 1:length(grid)){
      lbp <- lbp3d(grid[[j]])
      addLine <- lbp2Feature(lbp, n.feature, label = positive.label)
      featureMatrix <- rbind.data.frame(featureMatrix, addLine)
      pdbname <- attr(lbp, 'name')
    }
    if(rotations > 0){
      ppdb <- bio3d::read.pdb(p_pdbs[i])
      for(k in 1:rotations){
        rppdb <- rotatePDB(ppdb, degree = sample(c(1:360), 3, replace = TRUE))
        grid <- form.grid(pdb = rppdb, pdb.name = pdbname, size = size, units = units, combine.chains = combine.chains,  show.grid = show.grid)
        for(j in 1:length(grid)){
          lbp <- lbp3d(grid[[j]])
          addLine <- lbp2Feature(lbp, n.feature, label = positive.label)
          featureMatrix <- rbind.data.frame(featureMatrix, addLine)
        }
      }
    }
  }
  for(i in 1:length(n_pdbs)){
    grid <- form.grid(pdb.file = n_pdbs[i], size = size, units = units, combine.chains = combine.chains,  show.grid = show.grid)
    for(j in 1:length(grid)){
      lbp <- lbp3d(grid[[j]])
      addLine <- lbp2Feature(lbp, n.feature, label = negative.label)
      featureMatrix <- rbind.data.frame(featureMatrix, addLine)
      pdbname <- attr(lbp, 'name')
    }
    if(rotations > 0){
      ppdb <- bio3d::read.pdb(n_pdbs[i])
      for(k in 1:rotations){
        rppdb <- rotatePDB(ppdb, degree = sample(c(1:360), 3, replace = TRUE))
        grid <- form.grid(pdb = rppdb, pdb.name = pdbname, size = size, units = units, combine.chains = combine.chains,  show.grid = show.grid)
        for(j in 1:length(grid)){
          lbp <- lbp3d(grid[[j]])
          addLine <- lbp2Feature(lbp, n.feature, label = negative.label)
          featureMatrix <- rbind.data.frame(featureMatrix, addLine)
        }
      }
    }
  }
  cat('Feature Matrix with ', sum(featureMatrix$Label == positive.label), ' positive entries, ',
      sum(featureMatrix$Label == negative.label), ' negative entries, and ', (ncol(featureMatrix) - 1), ' features is gerneated. \n',
      sep = '')
  return(featureMatrix)

}



#' Separate featureMatrix into Train and Test sets
#'
#' This function takes the data.frame generated from featureMatrix and separates
#'  the data.frame in two Train and Test sets for cross-validation.
#'
#' @param featureMatrix The data.frame that featureMatrix function outputs.
#' @param positive.label The positive label used in featureMatrix function.
#' @param negative.label The negative label used in featureMatrix function.
#' @param testFold Portion of data.frame to be used as Test set.
#' @param labelAsFactor Set the column of Label into factors.
#' @return A list containing the Train and Test sets. no.Aug attribute contains
#' the number of base and augmented sample in each group.
#' @export
trainTest <- function(featureMatrix,
                      positive.label = 1,
                      negative.label = 0,
                      testFold = 5,
                      labelAsFactor = TRUE){
  pSet <- featureMatrix[featureMatrix$Label == positive.label,]
  pNames <- unique(sapply(strsplit(rownames(pSet),"_",),'[',1))
  nSet <- featureMatrix[featureMatrix$Label == negative.label,]
  nNames <- unique(sapply(strsplit(rownames(nSet),"_",),'[',1))
  n.pTest <- floor(length(pNames)/testFold)
  n.nTest <- floor(length(nNames)/testFold)
  #save parameters
  no.aug <- length(grep(pNames[1], rownames(pSet)))
  #separate positive train and test
  ran.pSample <- sample(c(1:length(pNames)),n.pTest)
  name.pTrain <- pNames[-c(ran.pSample)]
  pSet.Train <- data.frame()
  for(i in 1:length(name.pTrain)){
    pSet.Train <- rbind.data.frame(pSet.Train, pSet[grep(name.pTrain[i], rownames(pSet)),])
  }
  name.pTest <- pNames[c(ran.pSample)]
  pSet.Test <- data.frame()
  for(i in 1:length(name.pTest)){
    pSet.Test <- rbind.data.frame(pSet.Test, pSet[grep(name.pTest[i], rownames(pSet)),])
  }
  #separate negative train and test
  ran.nSample <- sample(c(1:length(nNames)),n.nTest)
  name.nTrain <- nNames[-c(ran.nSample)]
  nSet.Train <- data.frame()
  for(i in 1:length(name.nTrain)){
    nSet.Train <- rbind.data.frame(nSet.Train, nSet[grep(name.nTrain[i], rownames(nSet)),])
  }
  name.nTest <- nNames[c(ran.nSample)]
  nSet.Test <- data.frame()
  for(i in 1:length(name.nTest)){
    nSet.Test <- rbind.data.frame(nSet.Test, nSet[grep(name.nTest[i], rownames(nSet)),])
  }

  Train <- rbind.data.frame(pSet.Train, nSet.Train)
  Test <- rbind.data.frame(pSet.Test, nSet.Test)
  if(labelAsFactor){
    Train$Label <- as.factor(Train$Label)
    Test$Label <- as.factor(Test$Label)
  }
  out <- list(Train = Train, Test = Test)
  attr(out,'no.Aug') <- no.aug
  return(out)
}

#' Get votes from augmented data
#'
#' This function takes the augmented prediction results and votes for the
#' consensus prediction result
#'
#' @param x The vector of prediction results.
#' @param no.Aug Number of augmented samples, includes the original sample.
#' @param threshold The proportion of positive votes to be seen as consensus
#' @param positive.label The positive label of the data set (used in
#' featureMatrix function).
#' @param negative.label The negative label of the data set (used in
#' featureMatrix function).
#' @return A vector of voted results.
#' @export
augVote <- function(x, no.Aug, threshold = 0.5, positive.label = 1, negative.label = 0){
len <- length(x)
augM <- matrix(ncol = no.Aug, nrow = len/no.Aug, data = x, byrow =TRUE)
result <- c()
for(i in 1:nrow(augM)){
  tn <- length(grep(positive.label, augM[i,]))
  if( tn > (ncol(augM)*threshold)){
    result <- c(result, positive.label)
  }else{
    result <- c(result, negative.label)
  }
}
return(result)
}



#' Train a adaboost model using feature Matrix
#'
#' This function takes the data set generated by featureMatrix and processed by
#' trainTest and trains for a optimized adaboost model.
#'
#' @param set Set of train and test data.frames, output of trainTest.
#' @param Train The train data.frame. Only enter when set is not provided.
#' @param Test The test data.frame. Only enter when set is not provided.
#' @param no.Aug Numebr of augmentations plus the original. Only enter when set
#' is not provided.
#' @param positive.label The positive label of the data set (used in
#' featureMatrix function).
#' @param negative.label The negative label of the data set (used in
#' featureMatrix function).
#' @param iteration Number of iterations in training round.
#' @param mfinal The mfinal value in boosting function. See more in
#' ?adabog::boosting
#' @param mfinal.sequence The list of mfinal values used. Length must be the
#' same as iteration. Only enter when mfinal is not provided.
#' @param threshold The proportion of positive votes to be seen as consensus.
#' @param verbose If true, prints the error of each iteration.
#' @return A vector of voted results.
#' @export
adaboostModel <- function(set = NA,
                          Train = NA,
                          Test = NA,
                          no.Aug = NA,
                          positive.label = 1,
                          negative.label = 0,
                          iteration = 10,
                          mfinal = 100,
                          mfinal.sequence = NA,
                          threshold = 0.5,
                          verbose = FALSE){
  if(!is.na(set)[1]){
    Train <- set$Train
    Test <- set$Test
    no.Aug <- attr(set,'no.Aug')
  }else if(is.na(Train) || is.na(Test) || is.na(no.Aug)){
    cat('Missing aruguments Train/Test/no.Aug \n')
    return()
  }
  loss <- Inf
  bestModel <- NA
  for(i in 1:iteration){
    if(!is.na(mfinal.sequence[1])){
      if(length(mfinal.sequence) == iteration){
        mfinal <- mfinal.sequence[i]
      }else{
        cat('mfinal.sequence length is different from iterations')
        return()
      }
    }
    model <- adabag::boosting(Label~.,data = Train, boos = TRUE, mfinal = mfinal)
    pred <- stats::predict(model, Test)
    error <- pred$error
    if(error < loss){
      loss <- error
      bestModel <- model
    }
    if(verbose){
      cat('Iteration ', i, ' completed with error: ', error,'. \n', sep = '')
    }
  }
  bestPred <- stats::predict(bestModel, Test)
  if(length(unique(bestPred$class)) == 1 || length(unique(Test$Label)) == 1){
    preAug_F1 <- NA
  }else{
    preAug_F1 <- MLmetrics::F1_Score(bestPred$class, Test$Label)
  }
  result <- augVote(bestPred$class, no.Aug, threshold = threshold, positive.label = positive.label, negative.label = negative.label)
  trueLabel <- augVote(Test$Label, no.Aug, positive.label = positive.label, negative.label = negative.label)
  if(length(unique(result)) == 1 || length(unique(trueLabel)) == 1){
    F1 <- NA
  }else{
    F1 <- MLmetrics::F1_Score(result, trueLabel)
  }
  out <- list(model = bestModel, pred = bestPred, result = result, F1 = F1, preAug_F1 = preAug_F1)
  return(out)
}


#' Train a set of Adaboost models
#'
#' Train Adaboost models for multiple rounds and get model F1 for each round.
#'
#' @param featureMatrix The data.frame that featureMatrix function outputs.
#' @param testFold Portion of data.frame to be used as Test set.
#' @param round Number of training rounds.
#' @param threshold The proportion of positive votes to be seen as consensus.
#' @param iteration Number of iterations in training round.
#' @param mfinal The mfinal value in boosting function. See more in
#' ?adabog::boosting
#' @param mfinal.sequence The list of mfinal values used. Length must be the
#' same as iteration. Only enter when mfinal is not provided.
#' @return The adaboost model with the highest augmented F1 value and vecters of
#'  all pre-augmented and augmented rounds.
#' @export
adaboostTrain <- function(featureMatrix,
                          testFold = 5,
                          rounds = 10,
                          threshold = 0.5,
                          iteration = 10,
                          mfinal = 100,
                          mfinal.sequence = NA){
  preAugs <- c()
  Augs <- c()
  for(i in 1:rounds){
    set <- trainTest(fM, testFold = testFold)
    ABM <- adaboostModel(set = set, threshold = threshold, iteration = iteration, verbose = F, mfinal = mfinal, mfinal.sequence = mfinal.sequence)
    preAug <- ABM$preAug_F1
    Aug <- ABM$F1
    preAugs <- c(preAugs, preAug)
    Augs <- c(Augs, Aug)
    if(i == 1){
      BABM <- ABM
    }
    if(Aug == max(Augs)){
      BABM <- ABM
    }
    cat('round ', i, ' preAug: ', preAug, ', Aug: ', Aug, '\n', sep = '')
  }
  attr(BABM,'preAugs') <- preAugs
  attr(BABM,'Augs') <- Augs
  return(BABM)
}
