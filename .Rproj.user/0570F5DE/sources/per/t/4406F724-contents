library(caret)
library(e1071)

crossValidation <- function(expr_training)
{
  # Stratified 10-fold CV for within-dataset experiment
  # input:
  # expr_training: the pathway expression profiles of training set
  # output:
  # the train control of Stratified 10-fold cross validation

  # Define training control (10-fold cross validation)
  folds <- 10
  cvIndex <- createFolds(factor(expr_training$class), folds, returnTrain=T)
  tc <- trainControl(index=cvIndex, method='cv', number=folds, classProbs=TRUE, search="grid")

  return(tc)
}
train_Control <- crossValidation(expr_training)
#train_Control

# Naive Bayes Classifier
evaluateExprNB <- function(expr_training, expr_test, train_Control)
{
  # Train the classifier on the training set and predict on the test set (validation)
  # To evaluate the performance of the model
  # expr_training: the pathway expression profiles of training set
  # expr_test: the pathway expression profiles of test set
  # train_Control: the train control of Stratified 10-fold cross validation
  # output:
  # the prediction results

  # Training model using naive bayes
  rfFit <- train(class ~ ., data=expr_training, method="naive_bayes", trControl=train_Control)

  # Predict test data based on model
  eTestSet <- predict(rfFit, newdata=expr_test)

  return(eTestSet)
}
classExprNB <- evaluateExprNB(expr_training, expr_test, train_Control)
#classExprNB

evaluateGreedyNB <- function(expr_training, expr_test, expr_crosstest, train_Control)
{
  # Feature selection and classification
  # expr_training: the pathway profiles of training set
  # expr_test: the pathway profiles of test set (validation)
  # expr_crosstest: the pathway profiles of cross test set
  # train_Control: the train control of Stratified 10-fold cross validation
  # output:
  # the prediction result

  eTestSet <- evaluateExprNB(expr_training[, c(names(expr_training),"class")], expr_test[, c(names(expr_test),"class")], train_Control)

  # Evaluating model accuracy using confusion matrix
  cm <- table(expr_test$class, eTestSet)
  cm <- confusionMatrix(cm)
  AUC <- cm$overall['Accuracy']

  # Feature selection
  # Pathways in training set were added sequentially to train the model
  # The added pathway markers was kept in the feature set if the AUC increased and was removed otherwise
  flag <- 1
  for (i in 2 : (ncol(expr_training)-1))
  {
    eTestSetTmp <- evaluateExprNB(expr_training[, c(names(expr_training)[c(flag, i)],"class")], expr_test[, c(names(expr_test)[c(flag, i)],"class")], train_Control)
    AUCTmp <- table(expr_test$class, eTestSetTmp)
    AUCTmp <- confusionMatrix(AUCTmp)
    AUCTmp <- AUCTmp$overall['Accuracy']
    if (AUC < AUCTmp)
    {
      flag <- c(flag, i)
      eTestSet <- eTestSetTmp
      AUC <- AUCTmp
    }
  }
  riskPathways <- names(expr_training)[flag]

  # Classification using the selected pathway features
  # The performance of the optimized classifier was evaluated on the test set by using the pathway markers in the best feature set
  model <- train(class ~ ., data=expr_training[, c(names(expr_training)[flag],"class")],
                 method="naive_bayes", trControl=train_Control)

  predRes <- predict(model, newdata=expr_crosstest[, c(names(expr_crosstest)[flag], "class")], type="raw")

  # Evaluating model accuracy using confusion matrix
  cmCross <- table(expr_crosstest$class, predRes)
  cmCross <- confusionMatrix(cmCross)
  AUCCross <- cmCross$overall['Accuracy']


  return(list(predRes, AUC, riskPathways, AUCCross))
}
evaluateNB <- evaluateGreedyNB(expr_training, expr_test, expr_crosstest, train_Control)
evaluateNB
