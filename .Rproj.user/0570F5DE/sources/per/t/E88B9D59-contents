# K-nearest neighbors Classifier
evaluateExprKNN <- function(expr_training, expr_test, train_Control)
{
  # Train the classifier on the training set and predict on the test set (validation)
  # To evaluate the performance of the model
  # expr_training: training set
  # expr_test: test set
  # train_Control: the train control of Stratified 10-fold cross validation
  # output:
  # the prediction results

  # Training model using K-nearest neighbors
  rfFit <- train(class ~ ., data=expr_training, method="knn", trControl=train_Control)

  # Predict test data based on model
  eTestSet <- predict(rfFit, newdata=expr_test, type="raw")

  return (eTestSet)
}
classExprKNN <- evaluateExprKNN(expr_training, expr_test, train_Control)
#classExprKNN

evaluateGreedyKNN <- function(expr_training, expr_test, expr_crosstest, train_Control)
{
  # Feature selection and classification
  # expr_training: the pathway profiles of training set
  # expr_test: the pathway profiles of test set (validation)
  # expr_crosstest: the pathway profiles of cross test set
  # train_Control: the train control of Stratified 10-fold cross validation
  # output:
  # the prediction result

  eTestSet <- evaluateExprKNN(expr_training[, c(names(expr_training),"class")], expr_test[, c(names(expr_test),"class")], train_Control)

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
    eTestSetTmp <- evaluateExprKNN(expr_training[, c(names(expr_training)[c(flag, i)],"class")], expr_test[, c(names(expr_test)[c(flag, i)],"class")], train_Control)
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
  #model <- evaluateExpr(expr_training[, c(names(expr_training)[flag],"class")], expr_crosstest[, c(names(expr_crosstest)[flag],"class")], train_Control)
  model <- train(class ~ ., data=expr_training[, c(names(expr_training)[flag],"class")],
                 method="knn", trControl=train_Control)

  predRes <- predict(model, newdata=expr_crosstest[, c(names(expr_crosstest)[flag], "class")], type="raw")

  # Evaluating model accuracy using confusion matrix
  cmCross <- table(expr_crosstest$class, predRes)
  cmCross <- confusionMatrix(cmCross)
  AUCCross <- cmCross$overall['Accuracy']


  return(list(predRes, AUC, riskPathways, AUCCross))
}
evaluateKNN <- evaluateGreedyKNN(expr_training, expr_test, expr_crosstest, train_Control)
evaluateKNN
