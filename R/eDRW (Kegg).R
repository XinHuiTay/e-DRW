library(lattice)
library(Matrix)
library(igraph)

# Load eDRWData
data("LungExpData_ID")
load(file="eDRWData.RData")
gMetabolic <- get("graphMetabolic", envir=eDRWData)
gNonMetabolic <- get("graphNonMetabolic", envir=eDRWData)
DirectGraph <- get("DirectGraph", envir=eDRWData)

getnormalizedMatrix <- function(LungExpData_ID, d=1)
{
  # Normalization of each gene in the expression profile
  # input:
  # LungExpData_ID: the gene expression profile (labelled with gene ID)
  # d: the dimension of the matrix, d = 1: normalization according to row, d = 2: normalization according to column
  # output:
  # the normalized expression profile

  norm_mat <- matrix(nrow=dim(LungExpData_ID)[1], ncol=dim(LungExpData_ID)[2], data=NA)
  rownames(norm_mat) <- rownames(LungExpData_ID)
  colnames(norm_mat) <- colnames(LungExpData_ID)
  mu <- apply(LungExpData_ID, d, mean)
  sigma <- apply(LungExpData_ID, d, sd)
  if(d == 1)
  {
    for (i in 1 : dim(LungExpData_ID)[d]){
      norm_mat[i,] <- (LungExpData_ID[i,] - mu[i])/sigma[i]
    }
  }
  if(d == 2)
  {
    for (i in 1 : dim(LungExpData_ID)[d]){
      norm_mat[,i] <- (LungExpData_ID[,i] - mu[i])/sigma[i]
    }
  }
  return(norm_mat)
}
normalizedMatrix <- getnormalizedMatrix(LungExpData_ID, d=1)
#normalizedMatrix

calculatePBCScore <- function(normalizedMatrix, LungDiseaseLabel)
{
  # Point Biserial Correlation (PBC) of each gene in the expression profile
  # input:
  # normalizedMatrix: the normalized expression profile
  # LungDiseaseLabel: the disease sample index (Correlation)
  # output:
  # the PBC score of each gene in the expression profile

  PCscore <- matrix(NA, nrow=nrow(normalizedMatrix), ncol=2)
  rownames(PCscore) <- rownames(normalizedMatrix)
  for (i in 1:nrow(normalizedMatrix))
  {
    PCscore_tmp <- cor.test(normalizedMatrix[i,], LungDiseaseLabel)
    PCscore[i, 1] <- PCscore_tmp$estimate
    PCscore[i, 2] <- PCscore_tmp$p.value
  }
  return(abs(PCscore))
}
PBCScore <- calculatePBCScore(normalizedMatrix, LungDiseaseLabel)
#PBCScore

calculatetScore <- function(normalizedMatrix, LungNormSample, LungDiseaseSample)
{
  # T-test statistic of each gene in the expression profile
  # input:
  # normalizedMatrix: the normalized expression profile
  # LungNormSample: the index of normal samples
  # LungDiseaseSample: the index of disease samples
  # output:
  # the t-test statistic of each gene in the expression profile

  tscore <- matrix(NA, nrow=nrow(normalizedMatrix), ncol=2)
  rownames(tscore) <- rownames(normalizedMatrix)
  for (i in 1:nrow(normalizedMatrix))
  {
    tscore_tmp <- t.test(normalizedMatrix[i, LungDiseaseSample], normalizedMatrix[i, LungNormSample], var.equal=TRUE)
    tscore[i, 1] <- tscore_tmp$statistic
    tscore[i, 2] <- tscore_tmp$p.value
  }
  return(tscore)
}
Tscore <- calculatetScore(normalizedMatrix, LungNormSample, LungDiseaseSample)
#Tscore

calculatePCTScore <- function(Tscore, PBCScore)
{
  # Point Biserial Correlation Coefficient (PBC) and T-test statistic of each gene in the expression profile
  # input:
  # Tscore: the T-test statistic of each gene in expression profile
  # PBCScore: the PBC score of each gene in expression profiles
  # output:
  # the Point Biserial Correlation Coefficient and T-test statistic (PCT) of each gene in the expression profile

  PCTscore <- matrix(nrow=dim(Tscore)[1], ncol=2, data=Tscore)
  rownames(PCTscore) <- rownames(Tscore)
  for (i in 1:nrow(PCTscore))
  {
    PCTscore[i, 1] <- (Tscore[i, 1]^2) + PBCScore[i, 1]
  }

  return(PCTscore)
}
PCTScore <- calculatePCTScore(Tscore, PBCScore)
#PCTScore

getVertexTScore <- function(PCTScore, DirectGraph, mRNANames_ID)
{
  # Get the PCT and t-test statistic of each gene in Graph
  # input:
  # PCTScore: the Point Biserial Correlation Coefficient and T-test statistic (PCT) of each gene in the expression profile
  # DirectGraph: the adjacency matrix of global pathway network
  # mRNANames_ID: mRNA names (gene ID) in the expression profile
  # output:
  # the PCT and t-test statistic of each gene in the global pathway network

  #get vertex sequences
  Ve <- V(DirectGraph)

  VertexTScore <- matrix(NA, nrow=length(Ve), ncol=2)
  rownames(VertexTScore) <- get.vertex.attribute(DirectGraph, "name", index=V(DirectGraph))
  for (i in 1 : length(mRNANames_ID))
  {
    tmpgeneID <- mRNANames_ID[i]
    if (tmpgeneID %in% rownames(VertexTScore))
    {
      VertexTScore[tmpgeneID, 1] <- PCTScore[i, 1]
      VertexTScore[tmpgeneID, 2] <- PCTScore[i, 2]
    }
  }
  Idx <- which(is.na(VertexTScore[,1]))
  VertexTScore[Idx, 1] = 0
  VertexTScore[Idx, 2] = 1

  return(VertexTScore)
}
vertexTScore <- getVertexTScore(PCTScore, DirectGraph, mRNANames_ID)
#vertexTScore

getInitialVertexWeightTScore <- function(vertexTScore)
{
  # Calculate initial weight of genes in global pathway network
  # input:
  # vertexTscore: the PCT of each gene in the global pathway network
  # output:
  # the initial weight of genes in global pathway network

  vertexTScore[,1] <- abs(vertexTScore[,1])
  vertexTScore[,1] <- (vertexTScore[,1] - min(vertexTScore[,1]))/(max(vertexTScore[,1]) - min(vertexTScore[,1]))

  return(vertexTScore)
}
initialWeight <- getInitialVertexWeightTScore(vertexTScore)
#initialWeight

getp0 <- function(initialWeight)
{
  # Calculate initial weight of genes in global pathway network
  # input:
  # initialWeight: the PCT of each gene in the global pathway network
  # output:
  # the initial weight of genes in global pathway network

  GeneWeight <- initialWeight[,1]

  return(GeneWeight)
}
p0 <- getp0(initialWeight)
#p0

getEntropyProb <- function(p0)
{
  # Calculate entropy probability vector for each gene in global pathway network
  # Entropy formula: -(pij*(log2(pij))), pij: probability of each gene (edge weight) across the graph
  # input:
  # p0: initial weight of genes in global pathway network
  # output:
  # the entropy of genes (probability) in global pathway network

  prob <- p0/sum(p0)
  entropy <- -(prob*(log2(prob)))
  entropy[is.nan(entropy)] = 0

  return(entropy)
}
entropyProb <- getEntropyProb(p0)
#entropyProb

igraph2Adjmatrix<-function(DirectGraph)
{
  # Convert igraph to row-normalized adjacency matrix
  # input:
  # DirectGraph: the global pathway network (328 KEGG pathways)
  # output:
  # the row-normalized adjacency matrix of global pathway network

  AdjM <- as_adjacency_matrix(DirectGraph)

  # row-normalization
  for (i in 1:dim(AdjM)[1]){
    sumr <- sum(AdjM[i,])
    if(sumr == 0)
    {
      AdjM[i,] <- numeric(length=length(AdjM[i,]))
    }
    if(sumr > 0)
    {
      AdjM[i,] <- AdjM[i,]/sumr
    }
  }
  return(AdjM)
}
W <- igraph2Adjmatrix(DirectGraph)
#W

EntAdjmatrix<-function(W)
{
  # Convert igraph to entropy of row-normalized adjacency matrix
  # input:
  # W: the row-normalized adjacency matrix
  # output:
  # the entropy row-normalized adjacency matrix of global pathway network
  # the entropy of edge-weighted graph

  # Convert W sparse matrix to matrix
  NetEntropy <- as.matrix(W)

  # get adjacency information entropy
  for (i in 1:dim(W)[1]){
    NetEntropy[i,] <- -(W[i,]*(log2(W[i,])))
  }
  # Replace Nan to 0
  NetEntropy[is.nan(NetEntropy)] <- 0
  # Convert matrix to sparse matrix
  NetEntropy <- as(NetEntropy, "sparseMatrix")

  return(NetEntropy)
}
NetEntropy <- EntAdjmatrix(W)
#NetEntropy

rw_direct <- function(NetEntropy, entropyProb, gamma=0.5)
{
  # Calculate the weight of genes based on random walk
  # input:
  # NetEntropy: the entropy row-normalized row-normalized adjacency matrix
  # entropyProb: the entropy of genes (probability) in global pathway network
  # gamma: the restart probability (0.1-0.9)
  # output:
  # the weight of nodes after random walk

  # Add Ground Node, construct new adjacent matrix
  newrow <- matrix(1,1,dim(NetEntropy)[2])
  rownames(newrow) <- c("GN")
  W1 <- rbind(NetEntropy, newrow)
  newcol <- matrix(1,dim(W1)[1],1)
  colnames(newcol) <- c("GN")
  WGN <- cbind(W1,newcol) # adjacency matrix after adding ground node

  entropyProb <- t(as.matrix(entropyProb))
  entropyProb <- cbind(entropyProb,0) # The initial probability of the ground node is 0
  colnames(entropyProb)[dim(entropyProb)[2]] = "GN"

  PT <- entropyProb

  k <- 0
  delta <- 1

  # iteration
  while(is.na(delta > 1e-10))
  {
    PT1 <- (1-gamma)*WGN
    PT2 <- PT1 %*% t(PT)
    PT3 <- (gamma * entropyProb)
    PT4 <- t(PT2) + PT3
    delta <- sum(abs(PT4 - PT))
    PT <- PT4
    k <- k + 1
    delta <- FALSE
  }
  cat('converged\n')

  PT <- t(PT)
  rownames(PT) <- NULL

  # distribute the probability of the ground node back to all other nodes
  PT[1:(dim(PT)[1]-1)] <- PT[1:(dim(PT)[1]-1)] + PT[dim(PT)[1]]/(dim(PT)[1]-1)
  res <- drop(PT[1:(dim(PT)[1]-1)])

  return(res)
}
eDRW <- rw_direct(NetEntropy, entropyProb, gamma=0.5)
#eDRW

getEntropyWeight <- function(eDRW)
{
  # Calculate the entropy weight of each gene in global pathway network
  # input:
  # eDRW: the weight of each gene after random walk in global pathway network
  # output:
  # the entropy weight of genes in global pathway network

  for (i in 1:length(eDRW))
  {
    weight <- 1-(eDRW)
  }
  sumt <- sum(weight)
  for (i in 1:length(weight))
  {
    EntWeight <- weight / sumt
  }
  return(EntWeight)
}
EntropyWeight <- getEntropyWeight(eDRW)
#EntropyWeight

getVertexWeightFromRandWalk <- function(eDRW, DirectGraph)
{
  # Get the weight of genes based on random walk
  # input:
  # eDRW: the weight of nodes after random walk
  # DirectGraph: the global pathway network
  # output:
  # the weight of genes based on random walk

  VertexWeight <- eDRW
  names(VertexWeight) <- get.vertex.attribute(DirectGraph, "name", index=V(DirectGraph))

  return(VertexWeight)
}
vertexWeight <- getVertexWeightFromRandWalk(eDRW, DirectGraph)
vertexWeight

getvertexEntWeightFromRandWalk <- function(EntropyWeight, DirectGraph)
{
  # Get the entropy weight of genes based on random walk
  # input:
  # EntropyWeight: the entropy weight of nodes after random walk
  # DirectGraph: the global pathway network
  # output:
  # the entropy weight of genes based on random walk

  vertexEntWeightKegg <- EntropyWeight
  names(vertexEntWeightKegg) <- get.vertex.attribute(DirectGraph, "name", index=V(DirectGraph))

  return(vertexEntWeightKegg)
}
vertexEntWeight <- getvertexEntWeightFromRandWalk(EntropyWeight, DirectGraph)
vertexEntWeight

getDataFrames <- function(normalizedMatrix)
{
  # Split data frames for classification
  # input:
  # normalized_Matrix: the normalised gene expression profile
  # output:
  # Data frames for training, test set and cross test set

  # Create random training, validation, and test sets
  df <- normalizedMatrix

  # Set the fractions of the dataframe you want to split into training, validation, and test
  fractionTraining   <- 0.60
  fractionValidation <- 0.20
  fractionTest       <- 0.20

  # Compute sample sizes
  sampleSizeTraining   <- floor(fractionTraining   * nrow(df))
  sampleSizeValidation <- floor(fractionValidation * nrow(df))
  sampleSizeTest       <- floor(fractionTest       * nrow(df))

  # Create the randomly-sampled indices for the dataframe.
  # Use setdiff() to avoid overlapping subsets of indices
  indicesTraining    <- sort(sample(seq_len(nrow(df)), size=sampleSizeTraining))
  indicesNotTraining <- setdiff(seq_len(nrow(df)), indicesTraining)
  indicesValidation  <- sort(sample(indicesNotTraining, size=sampleSizeValidation))
  indicesTest        <- setdiff(indicesNotTraining, indicesValidation)

  # Output the three dataframes for training, validation and test
  matrix_training <- df[indicesTraining, ]
  matrix_test <- df[indicesValidation, ]
  matrix_crosstest <- df[indicesTest, ]

  return(list(matrix_training, matrix_test, matrix_crosstest))
}
dataframes <- getDataFrames(normalizedMatrix)
dataframes

getTrainSet <- function(dataframes)
{
  # Extract training set from data frames
  # input:
  # dataframes: data frames for training, test set and cross test set
  # output:
  # Data frames for training set

  TrainExpr <- dataframes[[1]]

  return(TrainExpr)
}
mRNA_matrix_training <- getTrainSet(dataframes)
mRNA_matrix_training

getTestSet <- function(dataframes)
{
  # Extract test set from data frames
  # input:
  # dataframes: data frames for training, test set and cross test set
  # output:
  # Data frames for test set

  TestExpr <- dataframes[[2]]

  return(TestExpr)
}
mRNA_matrix_test <- getTestSet(dataframes)
mRNA_matrix_test

getCrossTestSet <- function(dataframes)
{
  # Extract cross test set from data frames
  # input:
  # dataframes: data frames for training, test set and cross test set
  # output:
  # Data frames for cross test set

  CrossTestExpr <- dataframes[[3]]

  return(CrossTestExpr)
}
mRNA_matrix_crosstest <- getCrossTestSet(dataframes)
mRNA_matrix_crosstest

getPathwayRWTrain <- function(gNonMetabolic, gMetabolic, mRNA_matrix_training, vertexWeight, vertexTScore, vertexEntWeight)
{
  # Infer pathway expression profile for training set
  # input:
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # mRNA_matrix_training: the expression profile of training set
  # vertexWeight: the weight of genes from random walk
  # vertexTScore: the t-test statistic of each gene in the global pathway network
  # vertexEntWeight: the entropy weight of genes from random walk
  # output:
  # the pathway expression profile for training set
  # the t-test statistics of each pathways for training set

  # infer non metabolic pathway expression profile
  pathwayRW_training_gNonMetabolic <- c()
  TValueRW_training_gNonMetabolic <- matrix(NA, nrow=length(gNonMetabolic), ncol=1)
  rownames(TValueRW_training_gNonMetabolic) <- names(gNonMetabolic)
  for (i in 1 : length(gNonMetabolic))
  {
    Vpathwayi <- get.vertex.attribute(gNonMetabolic[[i]], "name", index=V(gNonMetabolic[[i]]))
    if (length(Vpathwayi) > 0)
    {
      n <- 0    # the number of differential genes in ith pathway
      pathway_training_tmp <- matrix(nrow=1, ncol=dim(mRNA_matrix_training)[2], data=0)
      TValueRW_training_tmp <- 0

      Idx_pathwayi <- c()
      for (j in 1 : length(Vpathwayi))
      {
        Idx <- which(rownames(mRNA_matrix_training)==Vpathwayi[j])
        if (length(Idx) > 0)
        {
          if (rownames(mRNA_matrix_training)[Idx] %in% names(vertexWeight))
          {
            if(vertexTScore[rownames(mRNA_matrix_training)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_training_tmp <- pathway_training_tmp + vertexTScore[rownames(mRNA_matrix_training)[Idx],1] * sum(vertexWeight[rownames(mRNA_matrix_training)[Idx]]) * mRNA_matrix_training[Idx,]
              TValueRW_training_tmp <- TValueRW_training_tmp + vertexWeight[rownames(mRNA_matrix_training)[Idx]] * abs(vertexTScore[rownames(mRNA_matrix_training)[Idx],1])
              n <- n + 1
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx)
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_training_tmp <- pathway_training_tmp / sqrt(sum(vertexEntWeight[rownames(mRNA_matrix_training)[Idx_pathwayi]]^2))
        TValueRW_training_gNonMetabolic[i, 1] <- TValueRW_training_tmp / sum(vertexWeight[rownames(mRNA_matrix_training)[Idx_pathwayi]])
        rownames(pathway_training_tmp) <- names(gNonMetabolic)[i]
        pathwayRW_training_gNonMetabolic <- rbind(pathwayRW_training_gNonMetabolic, pathway_training_tmp)
      }#else
    }#else
  }


  # infer metabolic pathway expression profile
  pathwayRW_training_gMetabolic <- c()
  TValueRW_training_gMetabolic <- matrix(NA, nrow=length(gMetabolic), ncol=1)
  rownames(TValueRW_training_gMetabolic) <- names(gMetabolic)
  for (i in 1 : length(gMetabolic))
  {
    Vpathwayi <- get.vertex.attribute(gMetabolic[[i]], "name", index=V(gMetabolic[[i]]))
    if (length(Vpathwayi) > 0)
    {
      n <- 0     # the number of differential genes in ith pathway
      pathway_training_tmp <- matrix(nrow=1, ncol=dim(mRNA_matrix_training)[2], data=0)
      TValueRW_training_tmp <- 0
      Idx_pathwayi <- c()   # record the index of nodes
      for (j in 1 : length(Vpathwayi))
      {
        Idx <- which(rownames(mRNA_matrix_training)==Vpathwayi[j])
        if (length(Idx) > 0)
        {
          if (rownames(mRNA_matrix_training)[Idx] %in% names(vertexWeight))
          {
            if(vertexTScore[rownames(mRNA_matrix_training)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_training_tmp <- pathway_training_tmp + vertexTScore[rownames(mRNA_matrix_training)[Idx],1] * sum(vertexWeight[rownames(mRNA_matrix_training)[Idx]]) * mRNA_matrix_training[Idx,]
              TValueRW_training_tmp <- TValueRW_training_tmp + vertexWeight[rownames(mRNA_matrix_training)[Idx]] * abs(vertexTScore[rownames(mRNA_matrix_training)[Idx],1])
              n <- n + 1
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx)
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_training_tmp <- pathway_training_tmp / sqrt(sum(vertexEntWeight[rownames(mRNA_matrix_training)[Idx_pathwayi]]^2))
        TValueRW_training_gMetabolic[i, 1] <- TValueRW_training_tmp / sum(vertexWeight[rownames(mRNA_matrix_training)[Idx_pathwayi]])
        rownames(pathway_training_tmp) <- names(gMetabolic)[i]
        pathwayRW_training_gMetabolic <- rbind(pathwayRW_training_gMetabolic, pathway_training_tmp)
      }#else
    }#else
  }

  # combine the nonMetabolic pathway and metabolic pathway expression profiles
  pathwayRW_training <- rbind(pathwayRW_training_gNonMetabolic, pathwayRW_training_gMetabolic)
  TValueRW_training <- rbind(TValueRW_training_gNonMetabolic, TValueRW_training_gMetabolic)
  TValueRW_training <- TValueRW_training[!is.na(TValueRW_training), ]

  Idx <- sort(TValueRW_training, decreasing = TRUE, index.return=TRUE)$ix
  TValueRW_training <- TValueRW_training[Idx]

  return(list(pathwayRW_training, TValueRW_training))
}
PathwayRWTrain <- getPathwayRWTrain(gNonMetabolic, gMetabolic, mRNA_matrix_training, vertexWeight, vertexTScore, vertexEntWeight)
#PathwayRWTrain

getPathwayRWTest <- function(gNonMetabolic, gMetabolic, mRNA_matrix_test, vertexWeight, vertexTScore, vertexEntWeight)
{
  # Infer pathway expression profile for test set
  # input:
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # mRNA_matrix_test: the expression profile of test set
  # vertexWeight: the weight of genes from random walk
  # vertexTScore: the t-test statistic of each gene in the global pathway network
  # vertexEntWeight: the entropy weight of genes from random walk
  # output:
  # the pathway expression profile of test set
  # the t-test statistics of each pathways for test set

  # infer non metabolic pathway expression profile
  pathwayRW_test_gNonMetabolic <- c()
  TValueRW_test_gNonMetabolic <- matrix(NA, nrow=length(gNonMetabolic), ncol=1)
  rownames(TValueRW_test_gNonMetabolic) <- names(gNonMetabolic);
  for (i in 1 : length(gNonMetabolic))
  {
    Vpathwayi <- get.vertex.attribute(gNonMetabolic[[i]], "name", index=V(gNonMetabolic[[i]]))
    if (length(Vpathwayi) > 0)
    {
      n <- 0    # the number of differential genes in ith pathway
      pathway_test_tmp <- matrix(nrow=1, ncol=dim(mRNA_matrix_test)[2], data=0)
      TValueRW_test_tmp <- 0

      Idx_pathwayi <- c()
      for (j in 1 : length(Vpathwayi))
      {
        Idx <- which(rownames(mRNA_matrix_test)==Vpathwayi[j])
        if (length(Idx) > 0)
        {
          if (rownames(mRNA_matrix_test)[Idx] %in% names(vertexWeight))
          {
            if(vertexTScore[rownames(mRNA_matrix_test)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_test_tmp <- pathway_test_tmp + vertexTScore[rownames(mRNA_matrix_test)[Idx],1] * sum(vertexWeight[rownames(mRNA_matrix_test)[Idx]]) * mRNA_matrix_test[Idx,]
              TValueRW_test_tmp <- TValueRW_test_tmp + vertexWeight[rownames(mRNA_matrix_test)[Idx]] * abs(vertexTScore[rownames(mRNA_matrix_test)[Idx],1])
              n <- n + 1;
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx);
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_test_tmp <- pathway_test_tmp / sqrt(sum(vertexEntWeight[rownames(mRNA_matrix_test)[Idx_pathwayi]]^2))
        TValueRW_test_gNonMetabolic[i, 1] <- TValueRW_test_tmp / sum(vertexWeight[rownames(mRNA_matrix_test)[Idx_pathwayi]])
        rownames(pathway_test_tmp) <- names(gNonMetabolic)[i]
        pathwayRW_test_gNonMetabolic <- rbind(pathwayRW_test_gNonMetabolic, pathway_test_tmp)
      }#else
    }#else
  }


  # infer metabolic pathway expression profile
  pathwayRW_test_gMetabolic <- c()
  TValueRW_test_gMetabolic <- matrix(NA, nrow=length(gMetabolic), ncol=1);
  rownames(TValueRW_test_gMetabolic) <- names(gMetabolic)
  for (i in 1 : length(gMetabolic))
  {
    Vpathwayi <- get.vertex.attribute(gMetabolic[[i]], "name", index=V(gMetabolic[[i]]))
    if (length(Vpathwayi) > 0)
    {
      n <- 0     # the number of differential genes in ith pathway
      pathway_test_tmp <- matrix(nrow=1, ncol=dim(mRNA_matrix_test)[2], data=0)
      TValueRW_test_tmp <- 0
      Idx_pathwayi <- c()   # record the index of nodes
      for (j in 1 : length(Vpathwayi))
      {
        Idx <- which(rownames(mRNA_matrix_test)==Vpathwayi[j])
        if (length(Idx) > 0)
        {
          if (rownames(mRNA_matrix_test)[Idx] %in% names(vertexWeight))
          {
            if(vertexTScore[rownames(mRNA_matrix_test)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_test_tmp <- pathway_test_tmp + vertexTScore[rownames(mRNA_matrix_test)[Idx],1] * sum(vertexWeight[rownames(mRNA_matrix_test)[Idx]]) * mRNA_matrix_test[Idx,]
              TValueRW_test_tmp <- TValueRW_test_tmp + vertexWeight[rownames(mRNA_matrix_test)[Idx]] * abs(vertexTScore[rownames(mRNA_matrix_test)[Idx],1])
              n <- n + 1;
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx);
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_test_tmp <- pathway_test_tmp / sqrt(sum(vertexEntWeight[rownames(mRNA_matrix_test)[Idx_pathwayi]]^2))
        TValueRW_test_gMetabolic[i, 1] <- TValueRW_test_tmp / sum(vertexWeight[rownames(mRNA_matrix_test)[Idx_pathwayi]])
        rownames(pathway_test_tmp) <- names(gMetabolic)[i]
        pathwayRW_test_gMetabolic <- rbind(pathwayRW_test_gMetabolic, pathway_test_tmp)
      }#else
    }#else
  }

  # combine the nonMetabolic pathway and metabolic pathway expression profiles
  pathwayRW_test <- rbind(pathwayRW_test_gNonMetabolic, pathwayRW_test_gMetabolic)
  TValueRW_test <- rbind(TValueRW_test_gNonMetabolic, TValueRW_test_gMetabolic)
  TValueRW_test <- TValueRW_test[!is.na(TValueRW_test), ]

  Idx <- sort(TValueRW_test, decreasing = TRUE, index.return=TRUE)$ix
  TValueRW_test <- TValueRW_test[Idx]

  return(list(pathwayRW_test, TValueRW_test))
}
PathwayRWTest <- getPathwayRWTest(gNonMetabolic, gMetabolic, mRNA_matrix_test, vertexWeight, vertexTScore, vertexEntWeight)
#PathwayRWTest

getPathwayRWCross <- function(gNonMetabolic, gMetabolic, mRNA_matrix_crosstest, vertexWeight, vertexTScore, vertexEntWeight)
{
  # Infer pathway expression profile for test set
  # input:
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # mRNA_matrix_test: the expression profile of test set
  # vertexWeight: the weight of genes from random walk
  # vertexTScore: the t-test statistic of each gene in the global pathway network
  # vertexEntWeight: the entropy weight of genes from random walk
  # output:
  # the pathway expression profile of test set
  # the t-test statistics of each pathways for test set

  # infer non metabolic pathway expression profile
  pathwayRW_crosstest_gNonMetabolic <- c()
  TValueRW_crosstest_gNonMetabolic <- matrix(NA, nrow=length(gNonMetabolic), ncol=1)
  rownames(TValueRW_crosstest_gNonMetabolic) <- names(gNonMetabolic);
  for (i in 1 : length(gNonMetabolic))
  {
    Vpathwayi <- get.vertex.attribute(gNonMetabolic[[i]], "name", index=V(gNonMetabolic[[i]]))
    if (length(Vpathwayi) > 0)
    {
      n <- 0    # the number of differential genes in ith pathway
      pathway_crosstest_tmp <- matrix(nrow=1, ncol=dim(mRNA_matrix_crosstest)[2], data=0)
      TValueRW_crosstest_tmp <- 0

      Idx_pathwayi <- c()
      for (j in 1 : length(Vpathwayi))
      {
        Idx <- which(rownames(mRNA_matrix_crosstest)==Vpathwayi[j])
        if (length(Idx) > 0)
        {
          if (rownames(mRNA_matrix_crosstest)[Idx] %in% names(vertexWeight))
          {
            if(vertexTScore[rownames(mRNA_matrix_crosstest)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_crosstest_tmp <- pathway_crosstest_tmp + vertexTScore[rownames(mRNA_matrix_crosstest)[Idx],1] * sum(vertexWeight[rownames(mRNA_matrix_crosstest)[Idx]]) * mRNA_matrix_crosstest[Idx,]
              TValueRW_crosstest_tmp <- TValueRW_crosstest_tmp + vertexWeight[rownames(mRNA_matrix_crosstest)[Idx]] * abs(vertexTScore[rownames(mRNA_matrix_crosstest)[Idx],1])
              n <- n + 1;
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx);
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_crosstest_tmp <- pathway_crosstest_tmp / sqrt(sum(vertexEntWeight[rownames(mRNA_matrix_crosstest)[Idx_pathwayi]]^2))
        TValueRW_crosstest_gNonMetabolic[i, 1] <- TValueRW_crosstest_tmp / sum(vertexWeight[rownames(mRNA_matrix_crosstest)[Idx_pathwayi]])
        rownames(pathway_crosstest_tmp) <- names(gNonMetabolic)[i]
        pathwayRW_crosstest_gNonMetabolic <- rbind(pathwayRW_crosstest_gNonMetabolic, pathway_crosstest_tmp)
      }#else
    }#else
  }


  # infer metabolic pathway expression profile
  pathwayRW_crosstest_gMetabolic <- c()
  TValueRW_crosstest_gMetabolic <- matrix(NA, nrow=length(gMetabolic), ncol=1);
  rownames(TValueRW_crosstest_gMetabolic) <- names(gMetabolic)
  for (i in 1 : length(gMetabolic))
  {
    Vpathwayi <- get.vertex.attribute(gMetabolic[[i]], "name", index=V(gMetabolic[[i]]))
    if (length(Vpathwayi) > 0)
    {
      n <- 0     # the number of differential genes in ith pathway
      pathway_crosstest_tmp <- matrix(nrow=1, ncol=dim(mRNA_matrix_crosstest)[2], data=0)
      TValueRW_crosstest_tmp <- 0
      Idx_pathwayi <- c()   # record the index of nodes
      for (j in 1 : length(Vpathwayi))
      {
        Idx <- which(rownames(mRNA_matrix_crosstest)==Vpathwayi[j])
        if (length(Idx) > 0)
        {
          if (rownames(mRNA_matrix_crosstest)[Idx] %in% names(vertexWeight))
          {
            if(vertexTScore[rownames(mRNA_matrix_crosstest)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_crosstest_tmp <- pathway_crosstest_tmp + vertexTScore[rownames(mRNA_matrix_crosstest)[Idx],1] * sum(vertexWeight[rownames(mRNA_matrix_crosstest)[Idx]]) * mRNA_matrix_crosstest[Idx,]
              TValueRW_crosstest_tmp <- TValueRW_crosstest_tmp + vertexWeight[rownames(mRNA_matrix_crosstest)[Idx]] * abs(vertexTScore[rownames(mRNA_matrix_crosstest)[Idx],1])
              n <- n + 1;
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx);
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_crosstest_tmp <- pathway_crosstest_tmp / sqrt(sum(vertexEntWeight[rownames(mRNA_matrix_crosstest)[Idx_pathwayi]]^2))
        TValueRW_crosstest_gMetabolic[i, 1] <- TValueRW_crosstest_tmp / sum(vertexWeight[rownames(mRNA_matrix_crosstest)[Idx_pathwayi]])
        rownames(pathway_crosstest_tmp) <- names(gMetabolic)[i]
        pathwayRW_crosstest_gMetabolic <- rbind(pathwayRW_crosstest_gMetabolic, pathway_crosstest_tmp)
      }#else
    }#else
  }

  # combine the nonMetabolic pathway and metabolic pathway expression profiles
  pathwayRW_crosstest <- rbind(pathwayRW_crosstest_gNonMetabolic, pathwayRW_crosstest_gMetabolic)
  TValueRW_crosstest <- rbind(TValueRW_crosstest_gNonMetabolic, TValueRW_crosstest_gMetabolic)
  TValueRW_crosstest <- TValueRW_crosstest[!is.na(TValueRW_crosstest), ]

  Idx <- sort(TValueRW_crosstest, decreasing = TRUE, index.return=TRUE)$ix
  TValueRW_crosstest <- TValueRW_crosstest[Idx]

  return(list(pathwayRW_crosstest, TValueRW_crosstest))
}
PathwayRWCross <- getPathwayRWCross(gNonMetabolic, gMetabolic, mRNA_matrix_crosstest, vertexWeight, vertexTScore, vertexEntWeight)
#PathwayRWCross

predictExprCross <- function(mRNA_matrix_training, mRNA_matrix_test,  mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
{
  # Label normal and disease samples index in pathway expression profile
  # mRNA_matrix_training: the expression profile of training set
  # mRNA_matrix_test: the expression profile of feature evaluation set
  # mRNA_matrix_crosstest: the expression profile of test set
  # LungNormSample: the index of normal samples in training set
  # LungDiseaseSample: the index of disease samples in training set
  # LungNormSample_test: the index of normal samples in test set
  # LungDiseaseSample_test: the index of disease samples in test set
  # LungNormSample_crosstest: the index of normal samples in cross test set
  # LungDiseaseSample_crosstest: the index of disease samples in cross test set
  # output:
  # the labelled normal and disease samples index of pathway expression profile

  tScoreTrain <- calculatetScore(mRNA_matrix_training, LungNormSample, LungDiseaseSample)
  tScoreTest <- calculatetScore(mRNA_matrix_test, LungNormSample_test, LungDiseaseSample_test)
  tScoreCrossTest <- calculatetScore(mRNA_matrix_crosstest, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  vertexTScoreTrain <- getVertexTScore(tScoreTrain, DirectGraph, rownames(mRNA_matrix_training))
  vertexTScoreTest <- getVertexTScore(tScoreTest, DirectGraph, rownames(mRNA_matrix_test))
  vertexTScoreCrossTest <- getVertexTScore(tScoreCrossTest, DirectGraph, rownames(mRNA_matrix_crosstest))

  # extract pathway expression profiles
  PathwayRWTrain <- getPathwayRWTrain(gNonMetabolic, gMetabolic, mRNA_matrix_training, vertexWeight, vertexTScoreTrain, vertexEntWeight)
  PathwayRWTest <- getPathwayRWTest(gNonMetabolic, gMetabolic, mRNA_matrix_test, vertexWeight, vertexTScoreTest, vertexEntWeight)
  PathwayRWCross <- getPathwayRWCross(gNonMetabolic, gMetabolic, mRNA_matrix_crosstest, vertexWeight, vertexTScoreCrossTest, vertexEntWeight)

  # calculate the t-test statistic of each pathway in the pathway expression profiles for training set
  tScoreP_pathway_train <- calculatetScore(PathwayRWTrain[[1]], LungNormSample, LungDiseaseSample)
  tScore_pathway_train <- tScoreP_pathway_train[ ,1]
  Idx <- sort(tScore_pathway_train, decreasing=TRUE, index.return=TRUE)$ix
  tScore_pathway_train <- tScore_pathway_train[Idx]

  # calculate the t-test statistic of each pathway in the pathway expression profiles for feature set
  tScoreP_pathway_test <- calculatetScore(PathwayRWTest[[1]], LungNormSample_test, LungDiseaseSample_test)
  tScore_pathway_test <- tScoreP_pathway_test[ ,1]
  Idx <- sort(tScore_pathway_test, decreasing=TRUE, index.return=TRUE)$ix
  tScore_pathway_test <- tScore_pathway_test[Idx]

  # calculate the t-test statistic of each pathway in the pathway expression profiles for test set
  tScoreP_pathway_cross <- calculatetScore(PathwayRWCross[[1]], LungNormSample_crosstest, LungDiseaseSample_crosstest)
  tScore_pathway_cross <- tScoreP_pathway_cross[ ,1]
  Idx <- sort(tScore_pathway_cross, decreasing=TRUE, index.return=TRUE)$ix
  tScore_pathway_cross <- tScore_pathway_cross[Idx]

  tScore_pathway <- c(tScore_pathway_train, tScore_pathway_test, tScore_pathway_cross)

  # make data for RWeka
  classType_training <- rep(NA, (length(LungNormSample) + length(LungDiseaseSample)))
  classType_training[LungNormSample] <- "normal"
  classType_training[LungDiseaseSample] <- "disease"
  arffRW_training <- data.frame(t(PathwayRWTrain[[1]]), "class"=classType_training, check.names=F)

  classType_test <- rep(NA, (length(LungNormSample_test) + length(LungDiseaseSample_test)))
  classType_test[LungNormSample_test] <- "normal"
  classType_test[LungDiseaseSample_test] <- "disease"
  arffRW_test <- data.frame(t(PathwayRWTest[[1]]), "class"=classType_test, check.names=F)

  classType_crosstest <- rep(NA, (length(LungNormSample_crosstest) + length(LungDiseaseSample_crosstest)))
  classType_crosstest[LungNormSample_crosstest] <- "normal"
  classType_crosstest[LungDiseaseSample_crosstest] <- "disease"
  arffRW_crosstest <- data.frame(t(PathwayRWCross[[1]]), "class" = classType_crosstest, check.names=F)

  resPredict <- list(arffRW_training, arffRW_test, arffRW_crosstest, tScore_pathway_train, tScore_pathway_test, tScore_pathway_cross)
  return(resPredict)
}
predictExprKegg <- predictExprCross(mRNA_matrix_training, mRNA_matrix_test,  mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
#predictExprKegg
