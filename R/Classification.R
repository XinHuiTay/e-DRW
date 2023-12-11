getpredictKeggPidKNN <- function(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                                 gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                                 LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
{
  # Classification using KNN classifier
  # normalizedMatrix: the normalised gene expression profile (KEGG)
  # normalized_Matrix: the normalised gene expression profile (PID)
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # vertexWeight: the weight of genes from random walk (KEGG)
  # vertexTScore: the t-test statistic of each gene in the global pathway network (KEGG)
  # vertexEntWeight: the entropy weight of each gene in the global pathway network (KEGG)
  # LungNormSample: the index of normal samples in training set
  # LungDiseaseSample: the index of disease samples in training set
  # LungNormSample_test: the index of normal samples in test set
  # LungDiseaseSample_test: the index of disease samples in test set
  # LungNormSample_crosstest: the index of normal samples in cross test set
  # LungDiseaseSample_crosstest: the index of disease samples in cross test set
  # gPid2: the nonMetabolic pathway
  # gPid1: the metabolic pathway
  # vertex_Weight: the weight of genes from random walk (PID)
  # vertex_TScore: the PCT of each gene in global pathway network (PID)
  # vertex_EntWeight: the entropy weight of each gene in global pathway network (PID)
  # output:
  # the predictions results using KNN classifier

  dataframes <- getDataFrames(normalizedMatrix)
  mRNA_matrix_training <- getTrainSet(dataframes)
  mRNA_matrix_test <- getTestSet(dataframes)
  mRNA_matrix_crosstest <- getCrossTestSet(dataframes)
  PathwayRWTrain <- getPathwayRWTrain(gNonMetabolic, gMetabolic, mRNA_matrix_training, vertexWeight, vertexTScore, vertexEntWeight)
  PathwayRWTest <- getPathwayRWTest(gNonMetabolic, gMetabolic, mRNA_matrix_test, vertexWeight, vertexTScore, vertexEntWeight)
  PathwayRWCross <- getPathwayRWCross(gNonMetabolic, gMetabolic, mRNA_matrix_crosstest, vertexWeight, vertexTScore, vertexEntWeight)
  predictExprKegg <- predictExprCross(mRNA_matrix_training, mRNA_matrix_test, mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)

  data_frames <- SplitDataFrames(normalized_Matrix)
  mRNA_matrix_training <- getTrainSetExpr(data_frames)
  mRNA_matrix_test <- getTestSetExpr(data_frames)
  mRNA_matrix_crosstest <- getCrossTestSetExpr(data_frames)
  PathwayRW_Train <- InferPathwayRWTrain(gPid2, gPid1, mRNA_matrix_training, vertex_Weight, vertex_TScore, vertex_EntWeight)
  PathwayRW_Test <- InferPathwayRWTest(gPid2, gPid1, mRNA_matrix_test, vertex_Weight, vertex_TScore, vertex_EntWeight)
  PathwayRW_Cross <- InferPathwayRWCross(gPid2, gPid1, mRNA_matrix_crosstest, vertex_Weight, vertex_TScore, vertex_EntWeight)
  predictExprPid <- predictExprCrossPid(mRNA_matrix_training, mRNA_matrix_test, mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)

  predictExpr <- predictExprKeggPid(predictExprKegg, predictExprPid)
  pathwayExpr <- pathwaySelection(predictExpr)
  expr_training <- pathwaySelectionTrain(pathwayExpr)
  expr_test <- pathwaySelectionTest(pathwayExpr)
  expr_crosstest <- pathwaySelectionCrossTest(pathwayExpr)
  train_Control <- crossValidation(expr_training)
  classExprKNN <- evaluateExprKNN(expr_training, expr_test, train_Control)
  evaluateKNN <- evaluateGreedyKNN(expr_training, expr_test, expr_crosstest, train_Control)

  return(evaluateKNN)
}
predictKeggPidKNN <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                                          gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                                          LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
#predictKeggPidKNN

getPredictCrossKNN <- function(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest,
                               gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight)
{
  # Classification across 10 experiments using KNN classifier
  # normalizedMatrix: the normalised gene expression profile (KEGG)
  # normalized_Matrix: the normalised gene expression profile (PID)
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # vertexWeight: the weight of genes from random walk (KEGG)
  # vertexTScore: the t-test statistic of each gene in the global pathway network (KEGG)
  # vertexEntWeight: the entropy weight of each gene in the global pathway network (KEGG)
  # LungNormSample: the index of normal samples in training set
  # LungDiseaseSample: the index of disease samples in training set
  # LungNormSample_test: the index of normal samples in test set
  # LungDiseaseSample_test: the index of disease samples in test set
  # LungNormSample_crosstest: the index of normal samples in cross test set
  # LungDiseaseSample_crosstest: the index of disease samples in cross test set
  # gPid2: the nonMetabolic pathway
  # gPid1: the metabolic pathway
  # vertex_Weight: the weight of genes from random walk (PID)
  # vertex_TScore: the PCT of each gene in global pathway network (PID)
  # vertex_EntWeight: the entropy weight of each gene in global pathway network (PID)
  # output:
  # sumRiskPathways: list of risk pathways detected across 10 experiments
  # sumAUC: list of AUC across 10 experiments
  # AUCKNN: mean AUC across 10 experiments

  KNN1 <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                               gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                               LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp1 <- KNN1[[3]]
  AUCExp1 <- KNN1[[4]]

  KNN2 <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                               gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                               LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp2 <- KNN2[[3]]
  AUCExp2 <- KNN2[[4]]

  KNN3 <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                               gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                               LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp3 <- KNN3[[3]]
  AUCExp3 <- KNN3[[4]]

  KNN4 <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                               gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                               LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp4 <- KNN4[[3]]
  AUCExp4 <- KNN4[[4]]

  KNN5 <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                               gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                               LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp5 <- KNN5[[3]]
  AUCExp5 <- KNN5[[4]]

  KNN6 <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                               gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                               LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp6 <- KNN6[[3]]
  AUCExp6 <- KNN6[[4]]

  KNN7 <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                               gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                               LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp7 <- KNN7[[3]]
  AUCExp7 <- KNN7[[4]]

  KNN8 <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                               gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                               LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp8 <- KNN8[[3]]
  AUCExp8 <- KNN8[[4]]

  KNN9 <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                               gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                               LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp9 <- KNN9[[3]]
  AUCExp9 <- KNN9[[4]]

  KNN10 <- getpredictKeggPidKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                                gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                                LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp10 <- KNN10[[3]]
  AUCExp10 <- KNN10[[4]]

  sumRiskPathways <- list(riskPathExp1, riskPathExp2, riskPathExp3, riskPathExp4, riskPathExp5, riskPathExp6, riskPathExp7, riskPathExp8, riskPathExp9, riskPathExp10)
  sumAUC <- c(AUCExp1, AUCExp2, AUCExp3, AUCExp4, AUCExp5, AUCExp6, AUCExp7, AUCExp8, AUCExp9, AUCExp10)
  AUCKNN <- mean(sumAUC)
  sumAUC <- as.list(sumAUC)

  return(list(sumRiskPathways, sumAUC, AUCKNN))
}
PredictCrossKNN <- getPredictCrossKNN(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest,
                                      gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight)
#PredictCrossKNN

getpredictKeggPidNB <- function(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                                gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                                LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
{
  # Classification using NB classifier
  # normalizedMatrix: the normalised gene expression profile (KEGG)
  # normalized_Matrix: the normalised gene expression profile (PID)
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # vertexWeight: the weight of genes from random walk (KEGG)
  # vertexTScore: the t-test statistic of each gene in the global pathway network (KEGG)
  # vertexEntWeight: the entropy weight of each gene in the global pathway network (KEGG)
  # LungNormSample: the index of normal samples in training set
  # LungDiseaseSample: the index of disease samples in training set
  # LungNormSample_test: the index of normal samples in test set
  # LungDiseaseSample_test: the index of disease samples in test set
  # LungNormSample_crosstest: the index of normal samples in cross test set
  # LungDiseaseSample_crosstest: the index of disease samples in cross test set
  # gPid2: the nonMetabolic pathway
  # gPid1: the metabolic pathway
  # vertex_Weight: the weight of genes from random walk (PID)
  # vertex_TScore: the PCT of each gene in global pathway network (PID)
  # vertex_EntWeight: the entropy weight of each gene in global pathway network (PID)
  # output:
  # the predictions results using NB classifier

  dataframes <- getDataFrames(normalizedMatrix)
  mRNA_matrix_training <- getTrainSet(dataframes)
  mRNA_matrix_test <- getTestSet(dataframes)
  mRNA_matrix_crosstest <- getCrossTestSet(dataframes)
  PathwayRWTrain <- getPathwayRWTrain(gNonMetabolic, gMetabolic, mRNA_matrix_training, vertexWeight, vertexTScore, vertexEntWeight)
  PathwayRWTest <- getPathwayRWTest(gNonMetabolic, gMetabolic, mRNA_matrix_test, vertexWeight, vertexTScore, vertexEntWeight)
  PathwayRWCross <- getPathwayRWCross(gNonMetabolic, gMetabolic, mRNA_matrix_crosstest, vertexWeight, vertexTScore, vertexEntWeight)
  predictExprKegg <- predictExprCross(mRNA_matrix_training, mRNA_matrix_test, mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)

  data_frames <- SplitDataFrames(normalized_Matrix)
  mRNA_matrix_training <- getTrainSetExpr(data_frames)
  mRNA_matrix_test <- getTestSetExpr(data_frames)
  mRNA_matrix_crosstest <- getCrossTestSetExpr(data_frames)
  PathwayRW_Train <- InferPathwayRWTrain(gPid2, gPid1, mRNA_matrix_training, vertex_Weight, vertex_TScore, vertex_EntWeight)
  PathwayRW_Test <- InferPathwayRWTest(gPid2, gPid1, mRNA_matrix_test, vertex_Weight, vertex_TScore, vertex_EntWeight)
  PathwayRW_Cross <- InferPathwayRWCross(gPid2, gPid1, mRNA_matrix_crosstest, vertex_Weight, vertex_TScore, vertex_EntWeight)
  predictExprPid <- predictExprCrossPid(mRNA_matrix_training, mRNA_matrix_test, mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)

  predictExpr <- predictExprKeggPid(predictExprKegg, predictExprPid)
  pathwayExpr <- pathwaySelection(predictExpr)
  expr_training <- pathwaySelectionTrain(pathwayExpr)
  expr_test <- pathwaySelectionTest(pathwayExpr)
  expr_crosstest <- pathwaySelectionCrossTest(pathwayExpr)
  train_Control <- crossValidation(expr_training)
  classExprNB <- evaluateExprNB(expr_training, expr_test, train_Control)
  evaluateNB <- evaluateGreedyNB(expr_training, expr_test, expr_crosstest, train_Control)

  return(evaluateNB)
}
predictKeggPidNB <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                                        gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                                        LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
#predictKeggPidNB

getPredictCrossNB <- function(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest,
                              gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight)
{
  # Classification across 10 experiments using NB classifier
  # normalizedMatrix: the normalised gene expression profile (KEGG)
  # normalized_Matrix: the normalised gene expression profile (PID)
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # vertexWeight: the weight of genes from random walk (KEGG)
  # vertexTScore: the t-test statistic of each gene in the global pathway network (KEGG)
  # vertexEntWeight: the entropy weight of each gene in the global pathway network (KEGG)
  # LungNormSample: the index of normal samples in training set
  # LungDiseaseSample: the index of disease samples in training set
  # LungNormSample_test: the index of normal samples in test set
  # LungDiseaseSample_test: the index of disease samples in test set
  # LungNormSample_crosstest: the index of normal samples in cross test set
  # LungDiseaseSample_crosstest: the index of disease samples in cross test set
  # gPid2: the nonMetabolic pathway
  # gPid1: the metabolic pathway
  # vertex_Weight: the weight of genes from random walk (PID)
  # vertex_TScore: the PCT of each gene in global pathway network (PID)
  # vertex_EntWeight: the entropy weight of each gene in global pathway network (PID)
  # output:
  # sumRiskPathways: list of risk pathways detected across 10 experiments
  # sumAUC: list of AUC across 10 experiments
  # AUCNB: mean AUC across 10 experiments

  NB1 <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp1 <- NB1[[3]]
  AUCExp1 <- NB1[[4]]

  NB2 <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp2 <- NB2[[3]]
  AUCExp2 <- NB2[[4]]

  NB3 <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp3 <- NB3[[3]]
  AUCExp3 <- NB3[[4]]

  NB4 <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp4 <- NB4[[3]]
  AUCExp4 <- NB4[[4]]

  NB5 <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp5 <- NB5[[3]]
  AUCExp5 <- NB5[[4]]

  NB6 <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp6 <- NB6[[3]]
  AUCExp6 <- NB6[[4]]

  NB7 <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp7 <- NB7[[3]]
  AUCExp7 <- NB7[[4]]

  NB8 <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp8 <- NB8[[3]]
  AUCExp8 <- NB8[[4]]

  NB9 <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp9 <- NB9[[3]]
  AUCExp9 <- NB9[[4]]

  NB10 <- getpredictKeggPidNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                              gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                              LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp10 <- NB10[[3]]
  AUCExp10 <- NB10[[4]]

  sumRiskPathways <- list(riskPathExp1, riskPathExp2, riskPathExp3, riskPathExp4, riskPathExp5, riskPathExp6, riskPathExp7, riskPathExp8, riskPathExp9, riskPathExp10)
  sumAUC <- c(AUCExp1, AUCExp2, AUCExp3, AUCExp4, AUCExp5, AUCExp6, AUCExp7, AUCExp8, AUCExp9, AUCExp10)
  AUCNB <- mean(sumAUC)
  sumAUC <- as.list(sumAUC)

  return(list(sumRiskPathways, sumAUC, AUCNB))
}
PredictCrossNB <- getPredictCrossNB(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest,
                                    gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight)
#PredictCrossNB

getpredictKeggPidLR <- function(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                                gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                                LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
{
  # Classification using LR classifier
  # normalizedMatrix: the normalised gene expression profile (KEGG)
  # normalized_Matrix: the normalised gene expression profile (PID)
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # vertexWeight: the weight of genes from random walk (KEGG)
  # vertexTScore: the t-test statistic of each gene in the global pathway network (KEGG)
  # vertexEntWeight: the entropy weight of each gene in the global pathway network (KEGG)
  # LungNormSample: the index of normal samples in training set
  # LungDiseaseSample: the index of disease samples in training set
  # LungNormSample_test: the index of normal samples in test set
  # LungDiseaseSample_test: the index of disease samples in test set
  # LungNormSample_crosstest: the index of normal samples in cross test set
  # LungDiseaseSample_crosstest: the index of disease samples in cross test set
  # gPid2: the nonMetabolic pathway
  # gPid1: the metabolic pathway
  # vertex_Weight: the weight of genes from random walk (PID)
  # vertex_TScore: the PCT of each gene in global pathway network (PID)
  # vertex_EntWeight: the entropy weight of each gene in global pathway network (PID)
  # output:
  # the predictions results using LR classifier

  dataframes <- getDataFrames(normalizedMatrix)
  mRNA_matrix_training <- getTrainSet(dataframes)
  mRNA_matrix_test <- getTestSet(dataframes)
  mRNA_matrix_crosstest <- getCrossTestSet(dataframes)
  PathwayRWTrain <- getPathwayRWTrain(gNonMetabolic, gMetabolic, mRNA_matrix_training, vertexWeight, vertexTScore, vertexEntWeight)
  PathwayRWTest <- getPathwayRWTest(gNonMetabolic, gMetabolic, mRNA_matrix_test, vertexWeight, vertexTScore, vertexEntWeight)
  PathwayRWCross <- getPathwayRWCross(gNonMetabolic, gMetabolic, mRNA_matrix_crosstest, vertexWeight, vertexTScore, vertexEntWeight)
  predictExprKegg <- predictExprCross(mRNA_matrix_training, mRNA_matrix_test, mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)

  data_frames <- SplitDataFrames(normalized_Matrix)
  mRNA_matrix_training <- getTrainSetExpr(data_frames)
  mRNA_matrix_test <- getTestSetExpr(data_frames)
  mRNA_matrix_crosstest <- getCrossTestSetExpr(data_frames)
  PathwayRW_Train <- InferPathwayRWTrain(gPid2, gPid1, mRNA_matrix_training, vertex_Weight, vertex_TScore, vertex_EntWeight)
  PathwayRW_Test <- InferPathwayRWTest(gPid2, gPid1, mRNA_matrix_test, vertex_Weight, vertex_TScore, vertex_EntWeight)
  PathwayRW_Cross <- InferPathwayRWCross(gPid2, gPid1, mRNA_matrix_crosstest, vertex_Weight, vertex_TScore, vertex_EntWeight)
  predictExprPid <- predictExprCrossPid(mRNA_matrix_training, mRNA_matrix_test, mRNA_matrix_crosstest, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)

  predictExpr <- predictExprKeggPid(predictExprKegg, predictExprPid)
  pathwayExpr <- pathwaySelection(predictExpr)
  expr_training <- pathwaySelectionTrain(pathwayExpr)
  expr_test <- pathwaySelectionTest(pathwayExpr)
  expr_crosstest <- pathwaySelectionCrossTest(pathwayExpr)
  train_Control <- crossValidation(expr_training)
  classExprLR <- evaluateExprLR(expr_training, expr_test, train_Control)
  evaluateLR <- evaluateGreedyLR(expr_training, expr_test, expr_crosstest, train_Control)

  return(evaluateLR)
}
predictKeggPidLR <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                                        gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                                        LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
#predictKeggPidLR

getPredictCrossLR <- function(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest,
                              gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight)
{
  # Classification across 10 experiments using LR classifier
  # normalizedMatrix: the normalised gene expression profile (KEGG)
  # normalized_Matrix: the normalised gene expression profile (PID)
  # gNonMetabolic: nonMetabolic pathway
  # gMetabolic: metabolic pathway
  # vertexWeight: the weight of genes from random walk (KEGG)
  # vertexTScore: the t-test statistic of each gene in the global pathway network (KEGG)
  # vertexEntWeight: the entropy weight of each gene in the global pathway network (KEGG)
  # LungNormSample: the index of normal samples in training set
  # LungDiseaseSample: the index of disease samples in training set
  # LungNormSample_test: the index of normal samples in test set
  # LungDiseaseSample_test: the index of disease samples in test set
  # LungNormSample_crosstest: the index of normal samples in cross test set
  # LungDiseaseSample_crosstest: the index of disease samples in cross test set
  # gPid2: the nonMetabolic pathway
  # gPid1: the metabolic pathway
  # vertex_Weight: the weight of genes from random walk (PID)
  # vertex_TScore: the PCT of each gene in global pathway network (PID)
  # vertex_EntWeight: the entropy weight of each gene in global pathway network (PID)
  # output:
  # sumRiskPathways: list of risk pathways detected across 10 experiments
  # sumAUC: list of AUC across 10 experiments
  # AUCLR: mean AUC across 10 experiments

  LR1 <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp1 <- LR1[[3]]
  AUCExp1 <- LR1[[4]]

  LR2 <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp2 <- LR2[[3]]
  AUCExp2 <- LR2[[4]]

  LR3 <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp3 <- LR3[[3]]
  AUCExp3 <- LR3[[4]]

  LR4 <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp4 <- LR4[[3]]
  AUCExp4 <- LR4[[4]]

  LR5 <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp5 <- LR5[[3]]
  AUCExp5 <- LR5[[4]]

  LR6 <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp6 <- LR6[[3]]
  AUCExp6 <- LR6[[4]]

  LR7 <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp7 <- LR7[[3]]
  AUCExp7 <- LR7[[4]]

  LR8 <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp8 <- LR8[[3]]
  AUCExp8 <- LR8[[4]]

  LR9 <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                             gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                             LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp9 <- LR9[[3]]
  AUCExp9 <- LR9[[4]]

  LR10 <- getpredictKeggPidLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight,
                              gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight,
                              LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest)
  riskPathExp10 <- LR10[[3]]
  AUCExp10 <- LR10[[4]]

  sumRiskPathways <- list(riskPathExp1, riskPathExp2, riskPathExp3, riskPathExp4, riskPathExp5, riskPathExp6, riskPathExp7, riskPathExp8, riskPathExp9, riskPathExp10)
  sumAUC <- c(AUCExp1, AUCExp2, AUCExp3, AUCExp4, AUCExp5, AUCExp6, AUCExp7, AUCExp8, AUCExp9, AUCExp10)
  AUCLR <- mean(sumAUC)
  sumAUC <- as.list(sumAUC)

  return(list(sumRiskPathways, sumAUC, AUCLR))
}
PredictCrossLR <- getPredictCrossLR(normalizedMatrix, normalized_Matrix, gNonMetabolic, gMetabolic, vertexWeight, vertexTScore, vertexEntWeight, LungNormSample, LungDiseaseSample, LungNormSample_test, LungDiseaseSample_test, LungNormSample_crosstest, LungDiseaseSample_crosstest,
                                    gPid2, gPid1, vertex_Weight, vertex_TScore, vertex_EntWeight)
#PredictCrossLR
