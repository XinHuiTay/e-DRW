library(lattice)
library(Matrix)
library(igraph)

# Load eDRWData
data("LungExpData_GS")
load(file="eDRWData.RData")
gPid1 <- get("graphPid1", envir=eDRWData)
gPid2 <- get("graphPid2", envir=eDRWData)
PidGraph <- get("PidGraph", envir=eDRWData)

normalisation <- function(LungExpData_GS, d=1)
{
  # Normalization of each gene in the expression profile
  # input:
  # LungExpData_GS: the gene expression profile (labelled with gene symbol)
  # d: the dimension of the matrix, d = 1: normalization according to row, d = 2: normalization according to column
  # output:
  # the normalized expression profile

  norm_mat <- matrix(nrow=dim(LungExpData_GS)[1], ncol=dim(LungExpData_GS)[2], data=NA)
  rownames(norm_mat) <- rownames(LungExpData_GS)
  colnames(norm_mat) <- colnames(LungExpData_GS)
  mu <- apply(LungExpData_GS, d, mean)
  sigma <- apply(LungExpData_GS, d, sd)
  if(d == 1)
  {
    for (i in 1 : dim(LungExpData_GS)[d]){
      norm_mat[i,] <- (LungExpData_GS[i,] - mu[i])/sigma[i]
    }
  }
  if(d == 2)
  {
    for (i in 1 : dim(LungExpData_GS)[d]){
      norm_mat[,i] <- (LungExpData_GS[,i] - mu[i])/sigma[i]
    }
  }
  return(norm_mat)
}
normalized_Matrix <- normalisation(LungExpData_GS, d=1)
#normalized_Matrix

PBCScoreCalculation <- function(normalized_Matrix, LungDiseaseLabel)
{
  # Point Biserial Correlation (PBC) of each gene in the expression profile
  # input:
  # normalized_Matrix: the normalized expression profile
  # LungDiseaseLabel: the disease sample index (Correlation)
  # output:
  # the PBC score of each gene in the expression profile

  PCscore <- matrix(NA, nrow=nrow(normalized_Matrix), ncol=2)
  rownames(PCscore) <- rownames(normalized_Matrix)
  for (i in 1:nrow(normalized_Matrix))
  {
    PCscore_tmp <- cor.test(normalized_Matrix[i,], LungDiseaseLabel)
    PCscore[i, 1] <- PCscore_tmp$estimate
    PCscore[i, 2] <- PCscore_tmp$p.value
  }
  return(abs(PCscore))
}
PBC_Score <- PBCScoreCalculation(normalized_Matrix, LungDiseaseLabel)
#PBC_Score

TScoreCalculation <- function(normalized_Matrix, LungNormSample, LungDiseaseSample)
{
  # T-test statistic of each gene in the expression profile
  # input:
  # normalized_Matrix: the normalized expression profile
  # LungNormSample: the index of normal samples
  # LungDiseaseSample: the index of disease samples
  # output:
  # the t-test statistic of each gene in the expression profile

  tscore <- matrix(NA, nrow=nrow(normalized_Matrix), ncol=2)
  rownames(tscore) <- rownames(normalized_Matrix)
  for (i in 1:nrow(normalized_Matrix))
  {
    tscore_tmp <- t.test(normalized_Matrix[i, LungDiseaseSample], normalized_Matrix[i, LungNormSample], var.equal=TRUE)
    tscore[i, 1] <- tscore_tmp$statistic
    tscore[i, 2] <- tscore_tmp$p.value
  }
  return(tscore)
}
T_score <- TScoreCalculation(normalized_Matrix, LungNormSample, LungDiseaseSample)
#T_score

PCTScoreCalculation <- function(T_score, PBC_Score)
{
  # Point Biserial Correlation Coefficient (PBC) and T-test statistic of each gene in the expression profile
  # input:
  # T_score: the T-test statistic of each gene in expression profile
  # PBC_Score: the PBC score of each gene in expression profiles
  # output:
  # the Point Biserial Correlation Coefficient and T-test statistic (PCT) of each gene in the expression profile

  PCTscore <- matrix(nrow=dim(T_score)[1], ncol=2, data=T_score)
  rownames(PCTscore) <- rownames(T_score)
  for (i in 1:nrow(PCTscore))
  {
    PCTscore[i, 1] <- (T_score[i, 1]^2) + PBC_Score[i, 1]
  }

  return(PCTscore)
}
PCT_Score <- PCTScoreCalculation(T_score, PBC_Score)
#PCT_Score

vertexTScoreCalculation <- function(PCT_Score, PidGraph, mRNANames_GS)
{
  # Get the PCT and t-test statistic of each gene in Graph
  # input:
  # PCT_Score: the Point Biserial Correlation Coefficient and T-test statistic (PCT) of each gene in the expression profile
  # PidGraph: the adjacency matrix of global pathway network
  # mRNANames_GS: mRNA names (gene symbol) in the expression profile
  # output:
  # the PCT and t-test statistic of each gene in the global pathway network

  #get vertex sequences
  Ve <- V(PidGraph)

  vertexTScore <- matrix(NA, nrow=length(Ve), ncol=2)
  rownames(vertexTScore) <- get.vertex.attribute(PidGraph, "name", index=V(PidGraph))
  for (i in 1 : length(mRNANames_GS))
  {
    tmpgeneID <- mRNANames_GS[i]
    if (tmpgeneID %in% rownames(vertexTScore))
    {
      vertexTScore[tmpgeneID, 1] <- PCT_Score[i, 1]
      vertexTScore[tmpgeneID, 2] <- PCT_Score[i, 2]
    }
  }
  Idx <- which(is.na(vertexTScore[,1]))
  vertexTScore[Idx, 1] = 0
  vertexTScore[Idx, 2] = 1

  return(vertexTScore)
}
vertex_TScore <- vertexTScoreCalculation(PCT_Score, PidGraph, mRNANames_GS)
#vertex_TScore

InitialvertexWeightTScoreCalculation <- function(vertex_TScore)
{
  # Calculate initial weight of genes in global pathway network
  # input:
  # vertex_TScore: the PCT of each gene in the global pathway network
  # output:
  # the initial weight of genes in global pathway network

  vertex_TScore[,1] <- abs(vertex_TScore[,1])
  vertex_TScore[,1] <- (vertex_TScore[,1] - min(vertex_TScore[,1]))/(max(vertex_TScore[,1]) - min(vertex_TScore[,1]))

  return(vertex_TScore)
}
initial_Weight <- InitialvertexWeightTScoreCalculation(vertex_TScore)
#initial_Weight

p0Calculation <- function(initial_Weight)
{
  # Calculate initial weight of genes in global pathway network
  # input:
  # initial_Weight: the PCT of each gene in the global pathway network
  # output:
  # the initial weight of genes in global pathway network

  Gene_Weight <- initial_Weight[,1]

  return(Gene_Weight)
}
p_0 <- p0Calculation(initial_Weight)
#p_0

entropyProbCalculation <- function(p_0)
{
  # Calculate entropy probability vector for each gene in global pathway network
  # Entropy formula: -(pij*(log2(pij))), pij: probability of each gene (edge weight) across the graph
  # input:
  # p_0: initial weight of genes in global pathway network
  # output:
  # the entropy of genes (probability) in global pathway network

  prob <- p_0/sum(p_0)
  entropy <- -(prob*(log2(prob)))
  entropy[is.nan(entropy)] = 0

  return(entropy)
}
entropy_Prob <- entropyProbCalculation(p_0)
#entropy_Prob

Pidgraph2Adjmatrix<-function(PidGraph)
{
  # Convert igraph to row-normalized adjacency matrix
  # input:
  # PidGraph: the global pathway network (208 PID pathways)
  # output:
  # the row-normalized adjacency matrix of global pathway network

  AdjM <- as_adjacency_matrix(PidGraph)

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
P <- Pidgraph2Adjmatrix(PidGraph)
P

PidEntAdjmatrix <-function(P)
{
  # Convert igraph to entropy of row-normalized adjacency matrix
  # input:
  # P: the row-normalized adjacency matrix
  # output:
  # the entropy row-normalized adjacency matrix of global pathway network
  # the entropy of edge-weighted graph

  # Convert W sparse matrix to matrix
  Net_Entropy <- as.matrix(P)

  # get adjacency information entropy
  for (i in 1:dim(P)[1]){
    Net_Entropy[i,] <- -(P[i,]*(log2(P[i,])))
  }
  # Replace Nan to 0
  Net_Entropy[is.nan(Net_Entropy)] <- 0
  # Convert matrix to sparse matrix
  Net_Entropy <- as(Net_Entropy, "sparseMatrix")

  return(Net_Entropy)
}
Net_Entropy <- PidEntAdjmatrix(P)
#Net_Entropy

rw_direct_pid <- function(Net_Entropy, entropy_Prob, gamma=0.5)
{
  # Calculate the weight of genes based on random walk
  # input:
  # Net_Entropy: the entropy row-normalized row-normalized adjacency matrix
  # entropy_Prob: the entropy of genes (probability) in global pathway network
  # gamma: the restart probability (0.1-0.9)
  # output:
  # the weight of nodes after random walk

  # Add Ground Node, construct new adjacent matrix
  newrow <- matrix(1,1,dim(Net_Entropy)[2])
  rownames(newrow) <- c("GN")
  W1 <- rbind(Net_Entropy, newrow)
  newcol <- matrix(1,dim(W1)[1],1)
  colnames(newcol) <- c("GN")
  WGN <- cbind(W1,newcol) # adjacency matrix after adding ground node

  entropy_Prob <- t(as.matrix(entropy_Prob))
  entropy_Prob <- cbind(entropy_Prob,0) # The initial probability of the ground node is 0
  colnames(entropy_Prob)[dim(entropy_Prob)[2]] = "GN"

  PT <- entropy_Prob

  k <- 0
  delta <- 1

  # iteration
  while(is.na(delta > 1e-10))
  {
    PT1 <- (1-gamma)*WGN
    PT2 <- PT1 %*% t(PT)
    PT3 <- (gamma * entropy_Prob)
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
e_DRW <- rw_direct_pid(Net_Entropy, entropy_Prob, gamma=0.5)
#e_DRW

EntropyWeightCalculation <- function(e_DRW)
{
  # Calculate the entropy weight of each gene in global pathway network
  # input:
  # e_DRW: the weight of each gene after random walk in global pathway network
  # output:
  # the entropy weight of genes in global pathway network

  for (i in 1:length(e_DRW))
  {
    weight <- 1-(e_DRW)
  }
  sumt <- sum(weight)
  for (i in 1:length(weight))
  {
    EntWeight <- weight / sumt
  }
  return(EntWeight)
}
Entropy_Weight <- EntropyWeightCalculation(e_DRW)
#Entropy_Weight

vertexWeightRWCalculation <- function(e_DRW, PidGraph)
{
  # Get the weight of genes based on random walk
  # input:
  # e_DRW: the weight of nodes after random walk
  # PidGraph: the global pathway network
  # output:
  # the weight of genes based on random walk

  vertexWeight <- e_DRW
  names(vertexWeight) <- get.vertex.attribute(PidGraph, "name", index=V(PidGraph))

  return(vertexWeight)
}
vertex_Weight <- vertexWeightRWCalculation(e_DRW, PidGraph)
#vertex_Weight

vertexEntWeightRWCalculation <- function(Entropy_Weight, PidGraph)
{
  # Get the entropy weight of genes based on random walk
  # input:
  # Entropy_Weight: the entropy weight of nodes after random walk
  # PidGraph: the global pathway network
  # output:
  # the entropy weight of genes based on random walk

  vertexEntWeightPid <- Entropy_Weight
  names(vertexEntWeightPid) <- get.vertex.attribute(PidGraph, "name", index=V(PidGraph))

  return(vertexEntWeightPid)
}
vertex_EntWeight <- vertexEntWeightRWCalculation(Entropy_Weight, PidGraph)
vertex_EntWeight

SplitDataFrames <- function(normalized_Matrix)
{
  # Split data frames for classification
  # input:
  # normalized_Matrix: the normalised gene expression profile
  # output:
  # Data frames for training, test set and cross test set

  # Create random training, test set and cross test set
  df1 <- normalized_Matrix

  # Set the fractions of the data frame you want to split into training, test set and cross test set
  fractionTraining   <- 0.60
  fractionValidation <- 0.20
  fractionTest       <- 0.20

  # Compute sample sizes
  sampleSizeTraining   <- floor(fractionTraining   * nrow(df1))
  sampleSizeValidation <- floor(fractionValidation * nrow(df1))
  sampleSizeTest       <- floor(fractionTest       * nrow(df1))

  # Create the randomly-sampled indices for the data frame
  # Use setdiff() to avoid overlapping subsets of indices
  indicesTraining    <- sort(sample(seq_len(nrow(df1)), size=sampleSizeTraining))
  indicesNotTraining <- setdiff(seq_len(nrow(df1)), indicesTraining)
  indicesValidation  <- sort(sample(indicesNotTraining, size=sampleSizeValidation))
  indicesTest        <- setdiff(indicesNotTraining, indicesValidation)

  # Output the three data frames for training, test set and cross test set
  matrix_training <- df1[indicesTraining, ]
  matrix_test <- df1[indicesValidation, ]
  matrix_crosstest <- df1[indicesTest, ]

  return(list(matrix_training, matrix_test, matrix_crosstest))
}
data_frames <- SplitDataFrames(normalized_Matrix)
data_frames

getTrainSetExpr <- function(data_frames)
{
  # Extract training set from data frames
  # input:
  # dataframes: data frames for training, test set and cross test set
  # output:
  # Data frames for training set

  TrainSetExpr <- data_frames[[1]]

  return(TrainSetExpr)
}
mRNA_matrix_training <- getTrainSetExpr(data_frames)
#mRNA_matrix_training

getTestSetExpr <- function(data_frames)
{
  # Extract test set from data frames
  # input:
  # dataframes: data frames for training, test set and cross test set
  # output:
  # Data frames for test set

  TestSetExpr <- data_frames[[2]]

  return(TestSetExpr)
}
mRNA_matrix_test <- getTestSetExpr(data_frames)
#mRNA_matrix_test

getCrossTestSetExpr <- function(data_frames)
{
  # Extract cross test set from data frames
  # input:
  # dataframes: data frames for training, test set and cross test set
  # output:
  # Data frames for cross test set

  CrossTestSetExpr <- data_frames[[3]]

  return(CrossTestSetExpr)
}
mRNA_matrix_crosstest <- getCrossTestSetExpr(data_frames)
#mRNA_matrix_crosstest

InferPathwayRWTrain <- function(gPid2, gPid1, mRNA_matrix_training, vertex_Weight, vertex_TScore, vertex_EntWeight)
{
  # Infer pathway expression profile for training set
  # input:
  # gPid2: the nonMetabolic pathway
  # gPid1: the metabolic pathway
  # mRNA_matrix_training: the expression profile of training set
  # vertex_Weight: the weight of genes from random walk
  # vertex_TScore: the PCT of each gene in global pathway network
  # vertex_EntWeight: the entropy weight of each gene in global pathway network
  # output:
  # the pathway expression profile for training set
  # the t-test statistics of each pathways for training set

  # infer non metabolic pathway expression profile
  pathwayRW_training_gPid2 <- c()
  TValueRW_training_gPid2 <- matrix(NA, nrow=length(gPid2), ncol=1)
  rownames(TValueRW_training_gPid2) <- names(gPid2)
  for (i in 1 : length(gPid2))
  {
    Vpathwayi <- get.vertex.attribute(gPid2[[i]], "name", index=V(gPid2[[i]]))
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
          if (rownames(mRNA_matrix_training)[Idx] %in% names(vertex_Weight))
          {
            if(vertex_TScore[rownames(mRNA_matrix_training)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_training_tmp <- pathway_training_tmp + vertex_TScore[rownames(mRNA_matrix_training)[Idx],1] * sum(vertex_Weight[rownames(mRNA_matrix_training)[Idx]]) * mRNA_matrix_training[Idx,]
              TValueRW_training_tmp <- TValueRW_training_tmp + vertex_Weight[rownames(mRNA_matrix_training)[Idx]] * abs(vertex_TScore[rownames(mRNA_matrix_training)[Idx],1])
              n <- n + 1
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx)
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_training_tmp <- pathway_training_tmp / sqrt(sum(vertex_EntWeight[rownames(mRNA_matrix_training)[Idx_pathwayi]]^2))
        TValueRW_training_gPid2[i, 1] <- TValueRW_training_tmp / sum(vertex_Weight[rownames(mRNA_matrix_training)[Idx_pathwayi]])
        rownames(pathway_training_tmp) <- names(gPid2)[i]
        pathwayRW_training_gPid2 <- rbind(pathwayRW_training_gPid2, pathway_training_tmp)
      }#else
    }#else
  }


  # infer metabolic pathway expression profile
  pathwayRW_training_gPid1 <- c()
  TValueRW_training_gPid1 <- matrix(NA, nrow=length(gPid1), ncol=1)
  rownames(TValueRW_training_gPid1) <- names(gPid1)
  for (i in 1 : length(gPid1))
  {
    Vpathwayi <- get.vertex.attribute(gPid1[[i]], "name", index=V(gPid1[[i]]))
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
          if (rownames(mRNA_matrix_training)[Idx] %in% names(vertex_Weight))
          {
            if(vertex_TScore[rownames(mRNA_matrix_training)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_training_tmp <- pathway_training_tmp + vertex_TScore[rownames(mRNA_matrix_training)[Idx],1] * sum(vertex_Weight[rownames(mRNA_matrix_training)[Idx]]) * mRNA_matrix_training[Idx,]
              TValueRW_training_tmp <- TValueRW_training_tmp + vertex_Weight[rownames(mRNA_matrix_training)[Idx]] * abs(vertex_TScore[rownames(mRNA_matrix_training)[Idx],1])
              n <- n + 1
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx)
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_training_tmp <- pathway_training_tmp / sqrt(sum(vertex_EntWeight[rownames(mRNA_matrix_training)[Idx_pathwayi]]^2))
        TValueRW_training_gPid1[i, 1] <- TValueRW_training_tmp / sum(vertex_Weight[rownames(mRNA_matrix_training)[Idx_pathwayi]])
        rownames(pathway_training_tmp) <- names(gPid1)[i]
        pathwayRW_training_gPid1 <- rbind(pathwayRW_training_gPid1, pathway_training_tmp)
      }#else
    }#else
  }

  # combine the nonMetabolic pathway and metabolic pathway expression profiles
  pathwayRW_training <- rbind(pathwayRW_training_gPid2, pathwayRW_training_gPid1)
  TValueRW_training <- rbind(TValueRW_training_gPid2, TValueRW_training_gPid1)
  TValueRW_training <- TValueRW_training[!is.na(TValueRW_training), ]

  Idx <- sort(TValueRW_training, decreasing = TRUE, index.return=TRUE)$ix
  TValueRW_training <- TValueRW_training[Idx]

  return(list(pathwayRW_training, TValueRW_training))
}
PathwayRW_Train <- InferPathwayRWTrain(gPid2, gPid1, mRNA_matrix_training, vertex_Weight, vertex_TScore, vertex_EntWeight)
#PathwayRW_Train

InferPathwayRWTest <- function(gPid2, gPid1, mRNA_matrix_test, vertex_Weight, vertex_TScore, vertex_EntWeight)
{
  # Infer pathway expression profile for test set
  # input:
  # gPid2: the nonMetabolic pathway
  # gPid1: the metabolic pathway
  # mRNA_matrix_test: the expression profile of test set
  # vertex_Weight: the weight of genes from random walk
  # vertex_TScore: the PCT of each gene in global pathway network
  # vertex_EntWeight: the entropy weight of each gene in global pathway network
  # output:
  # the pathway expression profile for test set
  # the t-test statistics of each pathways for test set

  # infer non metabolic pathway expression profile
  pathwayRW_test_gPid2 <- c()
  TValueRW_test_gPid2 <- matrix(NA, nrow=length(gPid2), ncol=1)
  rownames(TValueRW_test_gPid2) <- names(gPid2);
  for (i in 1 : length(gPid2))
  {
    Vpathwayi <- get.vertex.attribute(gPid2[[i]], "name", index=V(gPid2[[i]]))
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
          if (rownames(mRNA_matrix_test)[Idx] %in% names(vertex_Weight))
          {
            if(vertex_TScore[rownames(mRNA_matrix_test)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_test_tmp <- pathway_test_tmp + vertex_TScore[rownames(mRNA_matrix_test)[Idx],1] * sum(vertex_Weight[rownames(mRNA_matrix_test)[Idx]]) * mRNA_matrix_test[Idx,]
              TValueRW_test_tmp <- TValueRW_test_tmp + vertex_Weight[rownames(mRNA_matrix_test)[Idx]] * abs(vertex_TScore[rownames(mRNA_matrix_test)[Idx],1])
              n <- n + 1;
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx);
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_test_tmp <- pathway_test_tmp / sqrt(sum(vertex_EntWeight[rownames(mRNA_matrix_test)[Idx_pathwayi]]^2))
        TValueRW_test_gPid2[i, 1] <- TValueRW_test_tmp / sum(vertex_Weight[rownames(mRNA_matrix_test)[Idx_pathwayi]])
        rownames(pathway_test_tmp) <- names(gPid2)[i]
        pathwayRW_test_gPid2 <- rbind(pathwayRW_test_gPid2, pathway_test_tmp)
      }#else
    }#else
  }


  # infer metabolic pathway expression profile
  pathwayRW_test_gPid1 <- c()
  TValueRW_test_gPid1 <- matrix(NA, nrow=length(gPid1), ncol=1);
  rownames(TValueRW_test_gPid1) <- names(gPid1)
  for (i in 1 : length(gPid1))
  {
    Vpathwayi <- get.vertex.attribute(gPid1[[i]], "name", index=V(gPid1[[i]]))
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
          if (rownames(mRNA_matrix_test)[Idx] %in% names(vertex_Weight))
          {
            if(vertex_TScore[rownames(mRNA_matrix_test)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_test_tmp <- pathway_test_tmp + vertex_TScore[rownames(mRNA_matrix_test)[Idx],1] * sum(vertex_Weight[rownames(mRNA_matrix_test)[Idx]]) * mRNA_matrix_test[Idx,]
              TValueRW_test_tmp <- TValueRW_test_tmp + vertex_Weight[rownames(mRNA_matrix_test)[Idx]] * abs(vertex_TScore[rownames(mRNA_matrix_test)[Idx],1])
              n <- n + 1;
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx);
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_test_tmp <- pathway_test_tmp / sqrt(sum(vertex_EntWeight[rownames(mRNA_matrix_test)[Idx_pathwayi]]^2))
        TValueRW_test_gPid1[i, 1] <- TValueRW_test_tmp / sum(vertex_Weight[rownames(mRNA_matrix_test)[Idx_pathwayi]])
        rownames(pathway_test_tmp) <- names(gPid1)[i]
        pathwayRW_test_gPid1 <- rbind(pathwayRW_test_gPid1, pathway_test_tmp)
      }#else
    }#else
  }

  # combine the nonMetabolic pathway and metabolic pathway expression profiles
  pathwayRW_test <- rbind(pathwayRW_test_gPid2, pathwayRW_test_gPid1)
  TValueRW_test <- rbind(TValueRW_test_gPid2, TValueRW_test_gPid1)
  TValueRW_test <- TValueRW_test[!is.na(TValueRW_test), ]

  Idx <- sort(TValueRW_test, decreasing = TRUE, index.return=TRUE)$ix
  TValueRW_test <- TValueRW_test[Idx]

  return(list(pathwayRW_test, TValueRW_test))
}
PathwayRW_Test <- InferPathwayRWTest(gPid2, gPid1, mRNA_matrix_test, vertex_Weight, vertex_TScore, vertex_EntWeight)
#PathwayRW_Test

InferPathwayRWCross <- function(gPid2, gPid1, mRNA_matrix_crosstest, vertex_Weight, vertex_TScore, vertex_EntWeight)
{
  # Infer pathway expression profile for cross test set
  # input:
  # gPid2: the nonMetabolic pathway
  # gPid1: the metabolic pathway
  # mRNA_matrix_crosstest: the expression profile of cross test set
  # vertex_Weight: the weight of genes from random walk
  # vertex_TScore: the PCT of each gene in global pathway network
  # vertex_EntWeight: the entropy weight of each gene in global pathway network
  # output:
  # the pathway expression profile for cross test set
  # the t-test statistics of each pathways for cross test set

  # infer non metabolic pathway expression profile
  pathwayRW_crosstest_gPid2 <- c()
  TValueRW_crosstest_gPid2 <- matrix(NA, nrow=length(gPid2), ncol=1)
  rownames(TValueRW_crosstest_gPid2) <- names(gPid2);
  for (i in 1 : length(gPid2))
  {
    Vpathwayi <- get.vertex.attribute(gPid2[[i]], "name", index=V(gPid2[[i]]))
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
          if (rownames(mRNA_matrix_crosstest)[Idx] %in% names(vertex_Weight))
          {
            if(vertex_TScore[rownames(mRNA_matrix_crosstest)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_crosstest_tmp <- pathway_crosstest_tmp + vertex_TScore[rownames(mRNA_matrix_crosstest)[Idx],1] * sum(vertex_Weight[rownames(mRNA_matrix_crosstest)[Idx]]) * mRNA_matrix_crosstest[Idx,]
              TValueRW_crosstest_tmp <- TValueRW_crosstest_tmp + vertex_Weight[rownames(mRNA_matrix_crosstest)[Idx]] * abs(vertex_TScore[rownames(mRNA_matrix_crosstest)[Idx],1])
              n <- n + 1;
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx);
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_crosstest_tmp <- pathway_crosstest_tmp / sqrt(sum(vertex_EntWeight[rownames(mRNA_matrix_crosstest)[Idx_pathwayi]]^2))
        TValueRW_crosstest_gPid2[i, 1] <- TValueRW_crosstest_tmp / sum(vertex_Weight[rownames(mRNA_matrix_crosstest)[Idx_pathwayi]])
        rownames(pathway_crosstest_tmp) <- names(gPid2)[i]
        pathwayRW_crosstest_gPid2 <- rbind(pathwayRW_crosstest_gPid2, pathway_crosstest_tmp)
      }#else
    }#else
  }


  # infer metabolic pathway expression profile
  pathwayRW_crosstest_gPid1 <- c()
  TValueRW_crosstest_gPid1 <- matrix(NA, nrow=length(gPid1), ncol=1);
  rownames(TValueRW_crosstest_gPid1) <- names(gPid1)
  for (i in 1 : length(gPid1))
  {
    Vpathwayi <- get.vertex.attribute(gPid1[[i]], "name", index=V(gPid1[[i]]))
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
          if (rownames(mRNA_matrix_crosstest)[Idx] %in% names(vertex_Weight))
          {
            if(vertex_TScore[rownames(mRNA_matrix_crosstest)[Idx],2] < 0.05) # p-value < 0.05
            {
              pathway_crosstest_tmp <- pathway_crosstest_tmp + vertex_TScore[rownames(mRNA_matrix_crosstest)[Idx],1] * sum(vertex_Weight[rownames(mRNA_matrix_crosstest)[Idx]]) * mRNA_matrix_crosstest[Idx,]
              TValueRW_crosstest_tmp <- TValueRW_crosstest_tmp + vertex_Weight[rownames(mRNA_matrix_crosstest)[Idx]] * abs(vertex_TScore[rownames(mRNA_matrix_crosstest)[Idx],1])
              n <- n + 1;
              Idx_pathwayi <- rbind(Idx_pathwayi, Idx);
            }
          }
        }
      }
      if(n > 0)
      {
        pathway_crosstest_tmp <- pathway_crosstest_tmp / sqrt(sum(vertex_EntWeight[rownames(mRNA_matrix_crosstest)[Idx_pathwayi]]^2))
        TValueRW_crosstest_gPid1[i, 1] <- TValueRW_crosstest_tmp / sum(vertex_Weight[rownames(mRNA_matrix_crosstest)[Idx_pathwayi]])
        rownames(pathway_crosstest_tmp) <- names(gPid1)[i]
        pathwayRW_crosstest_gPid1 <- rbind(pathwayRW_crosstest_gPid1, pathway_crosstest_tmp)
      }#else
    }#else
  }

  # combine the nonMetabolic pathway and metabolic pathway expression profiles
  pathwayRW_crosstest <- rbind(pathwayRW_crosstest_gPid2, pathwayRW_crosstest_gPid1)
  TValueRW_crosstest <- rbind(TValueRW_crosstest_gPid2, TValueRW_crosstest_gPid1)
  TValueRW_crosstest <- TValueRW_crosstest[!is.na(TValueRW_crosstest), ]

  Idx <- sort(TValueRW_crosstest, decreasing = TRUE, index.return=TRUE)$ix
  TValueRW_crosstest <- TValueRW_crosstest[Idx]

  return(list(pathwayRW_crosstest, TValueRW_crosstest))
}
PathwayRW_Cross <- InferPathwayRWCross(gPid2, gPid1, mRNA_matrix_crosstest, vertex_Weight, vertex_TScore, vertex_EntWeight)
#PathwayRW_Cross

predictExprCrossPid <- function(mRNA_matrix_training, mRNA_matrix_test, mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
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

  tScoreTrain <- TScoreCalculation(mRNA_matrix_training, LungNormSample, LungDiseaseSample)
  tScoreTest <- TScoreCalculation(mRNA_matrix_test, LungNormSample_test, LungDiseaseSample_test)
  tScoreCrossTest <- TScoreCalculation(mRNA_matrix_crosstest, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  vertexTScoreTrain <- vertexTScoreCalculation(tScoreTrain, PidGraph, rownames(mRNA_matrix_training))
  vertexTScoreTest <- vertexTScoreCalculation(tScoreTest, PidGraph, rownames(mRNA_matrix_test))
  vertexTScoreCrossTest <- vertexTScoreCalculation(tScoreCrossTest, PidGraph, rownames(mRNA_matrix_crosstest))

  # extract pathway expression profiles
  PathwayRWTrain <- InferPathwayRWTrain(gPid2, gPid1, mRNA_matrix_training, vertex_Weight, vertexTScoreTrain, vertex_EntWeight)
  PathwayRWTest <- InferPathwayRWTest(gPid2, gPid1, mRNA_matrix_test, vertex_Weight, vertexTScoreTest, vertex_EntWeight)
  PathwayRWCross <- InferPathwayRWCross(gPid2, gPid1, mRNA_matrix_crosstest, vertex_Weight, vertexTScoreCrossTest, vertex_EntWeight)

  # calculate the t-test statistic of each pathway in the pathway expression profiles for training set
  tScoreP_pathway_train <- TScoreCalculation(PathwayRWTrain[[1]], LungNormSample, LungDiseaseSample)
  tScore_pathway_train <- tScoreP_pathway_train[ ,1]
  Idx <- sort(tScore_pathway_train, decreasing=TRUE, index.return=TRUE)$ix
  tScore_pathway_train <- tScore_pathway_train[Idx]

  # calculate the t-test statistic of each pathway in the pathway expression profiles for test set
  tScoreP_pathway_test <- TScoreCalculation(PathwayRWTest[[1]], LungNormSample_test, LungDiseaseSample_test)
  tScore_pathway_test <- tScoreP_pathway_test[ ,1]
  Idx <- sort(tScore_pathway_test, decreasing=TRUE, index.return=TRUE)$ix
  tScore_pathway_test <- tScore_pathway_test[Idx]

  # calculate the t-test statistic of each pathway in the pathway expression profiles for cross test set
  tScoreP_pathway_cross <- TScoreCalculation(PathwayRWCross[[1]], LungNormSample_crosstest, LungDiseaseSample_crosstest)
  tScore_pathway_cross <- tScoreP_pathway_cross[ ,1]
  Idx <- sort(tScore_pathway_cross, decreasing=TRUE, index.return=TRUE)$ix
  tScore_pathway_cross <- tScore_pathway_cross[Idx]

  tScore_pathway <- c(tScore_pathway_train, tScore_pathway_test, tScore_pathway_cross)

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
predictExprPid <- predictExprCrossPid(mRNA_matrix_training, mRNA_matrix_test, mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
#predictExprPid
