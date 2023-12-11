predictExprKeggPid <- function(predictExprKegg, predictExprPid)
{
  # Combine Kegg and Pid pathway expression profiles
  # predictExprKegg: the pathway expression profiles of KEGG pathways
  # predictExprPid: the pathway expression profiles of PID pathways
  # output:
  # the combined pathway expression profiles (KEGG-PID)

  # Combine pathway expression profiles for training set
  predictExpr <- predictExprKegg[[1]]
  x <- predictExpr[-length(predictExpr)]
  expr_training <- cbind.data.frame(x, predictExprPid[[1]])

  # Combine pathway expression profiles for validation set
  predictExpr <- predictExprKegg[[2]]
  y <- predictExpr[-length(predictExpr)]
  expr_test <- cbind.data.frame(y, predictExprPid[[2]])

  # Combine pathway expression profiles for test set
  predictExpr <- predictExprKegg[[3]]
  z <- predictExpr[-length(predictExpr)]
  expr_crosstest <- cbind.data.frame(z, predictExprPid[[3]])

  # Combine t-test statistic of each pathway for training set
  TvalueTrain <- append(predictExprKegg[[4]], predictExprPid[[4]])
  Idx <- sort(TvalueTrain, decreasing=TRUE, index.return=TRUE)$ix
  TvalueTrain <- TvalueTrain[Idx]

  # Combine t-test statistic of each pathway for validation set
  TvalueTest <- append(predictExprKegg[[5]], predictExprPid[[5]])
  Idx <- sort(TvalueTest, decreasing=TRUE, index.return=TRUE)$ix
  TvalueTest <- TvalueTest[Idx]

  # Combine t-test statistic of each pathway for test set
  TvalueCrossTest <- append(predictExprKegg[[6]], predictExprPid[[6]])
  Idx <- sort(TvalueCrossTest, decreasing=TRUE, index.return=TRUE)$ix
  TvalueCrossTest <- TvalueCrossTest[Idx]

  # Create list of pathway expression profiles
  resPredictExpr <- list(expr_training, expr_test, expr_crosstest, TvalueTrain, TvalueTest, TvalueCrossTest)

  return(resPredictExpr)
}
predictExpr <- predictExprKeggPid(predictExprKegg, predictExprPid)
#predictExpr

pathwaySelection <- function(predictExpr)
{
  # Select top 50 pathways of training set, test set and cross test set for classification
  # predictExpr: the combined pathway expression profiles (KEGG-PID)
  # output:
  # the pathway expression profile (top 50 pathways) of training set, test set and cross test set

  # The pathway expression profile of training set, test set and cross test set
  expr_training <- predictExpr[[1]]
  expr_test <- predictExpr[[2]]
  expr_crosstest <- predictExpr[[3]]

  # The t-test statistic of each pathway for training set
  Tvalue <- predictExpr[[4]]

  # Sort pathway expression profile of training set according to t-test statistic (descending order)
  Idx <- which(names(expr_training) %in% names(Tvalue))
  expr_training <- expr_training[, c(names(expr_training[Idx]), "class")]
  expr_training$class <- factor(expr_training$class)
  expr_training <- expr_training[, c(names(Tvalue), "class")]

  # Sort pathway expression profile of validation set according to t-test statistic (descending order)
  Idx <- which(names(Tvalue) %in% names(expr_test))
  names(Tvalue[Idx])
  expr_test <- expr_test[, c(names(Tvalue[Idx]), "class")]
  expr_test$class <- factor(expr_test$class)

  # Sort pathway expression profile of test set according to t-test statistic (descending order)
  Idx <- which(names(Tvalue) %in% names(expr_crosstest))
  names(Tvalue[Idx])
  expr_crosstest <- expr_crosstest[, c(names(Tvalue[Idx]), "class")]
  expr_crosstest$class <- factor(expr_crosstest$class)

  train <- names(expr_training)
  test <- names(expr_test)
  cross <- names(expr_crosstest)

  # Select common values in training set, test set and cross test set
  intersect_data <- intersect(intersect(train, test), cross)

  Idx <- which(names(expr_training) %in% (intersect_data))
  expr_training <- expr_training[, c(names(expr_training[Idx]))]
  expr_training$class <- factor(expr_training$class)

  Idx <- which(names(expr_test) %in% (intersect_data))
  expr_test <- expr_test[, c(names(expr_test[Idx]))]
  expr_test$class <- factor(expr_test$class)

  Idx <- which(names(expr_crosstest) %in% (intersect_data))
  expr_crosstest <- expr_crosstest[, c(names(expr_crosstest[Idx]))]
  expr_crosstest$class <- factor(expr_crosstest$class)

  # Select top 50 pathways from training, validation and test set as candidate features for classification
  expr_training <- expr_training[, c(names(expr_training[1:50]), "class")]
  expr_training$class <- factor(expr_training$class)

  expr_test <- expr_test[, c(names(expr_test[1:50]), "class")]
  expr_test$class <- factor(expr_test$class)

  expr_crosstest <- expr_crosstest[, c(names(expr_crosstest[1:50]), "class")]
  expr_crosstest$class <- factor(expr_crosstest$class)

  return(list(expr_training, expr_test, expr_crosstest))
}
pathwayExpr <- pathwaySelection(predictExpr)
#pathwayExpr

pathwaySelectionTrain <- function(pathwayExpr)
{
  # The pathway expression profile (top 50 pathways) of training set
  # pathwayExpr: the pathway expression profile (top 50 pathways) of training set, test set and cross test set
  # output:
  # the pathway expression profile (top 50 pathways) of training set

  pathwayExpr_train <- pathwayExpr[[1]]

  return(pathwayExpr_train)
}
expr_training <- pathwaySelectionTrain(pathwayExpr)
#expr_training

pathwaySelectionTest <- function(pathwayExpr)
{
  # The pathway expression profile (top 50 pathways) of test set
  # pathwayExpr: the pathway expression profile (top 50 pathways) of training set, test set and cross test set
  # output:
  # the pathway expression profile (top 50 pathways) of test set

  pathwayExpr_test <- pathwayExpr[[2]]

  return(pathwayExpr_test)
}
expr_test <- pathwaySelectionTest(pathwayExpr)
#expr_test

pathwaySelectionCrossTest <- function(pathwayExpr)
{
  # The pathway expression profile (top 50 pathways) of cross test set
  # pathwayExpr: the pathway expression profile (top 50 pathways) of training set, test set and cross test set
  # output:
  # the pathway expression profile (top 50 pathways) of cross test set

  pathwayExpr_crosstest <- pathwayExpr[[3]]

  return(pathwayExpr_crosstest)
}
expr_crosstest <- pathwaySelectionCrossTest(pathwayExpr)
#expr_crosstest
