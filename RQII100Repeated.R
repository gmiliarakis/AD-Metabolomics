require(caret)
require(rags2ridges)
library(pROC)
# set y as a factor
data(ADdata)
X <- as.matrix(scale(t(ADmetabolites)))
y <- as.factor(APOE$APOE)

#Encode y so valid R variables can be created from its classes
levels(y) <- make.names(levels(y))
X <- thomson
#Merge X and y into df
df <- cbind.data.frame(X, y)

# Set the random seed for reproducibility
set.seed(666)
times= 10

# Define the training control hyperparameters
ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  sampling = 'smote'
)
auc.lr <- numeric(times)
metrics.lr  <- matrix(nrow = 11, ncol = times)
# Perform 10-fold CV repeated 100 times
for (i in 1:times) {
  set.seed(666 + i)  # Set a different seed for each iteration
  
  # Create a matrix with the training indexes
  split_index <- createDataPartition(df$y, p = 0.7, list = FALSE)
  
  # Create the training set
  train <- df[split_index,]
  
  # Create the test set
  test <- df[-split_index,]
  
  # Store the test y in labels
  labels <- test$y
  labels.n <- as.numeric(labels) - 1
  
  # Train the logistic regression model
  lr <- train(y ~ .,
              train,
              method = "glm",
              metric = "ROC",
              trControl = ctrl)
  
  yhat.lr <- predict(lr, test)
  
  # Create a confusion matrix and get performance metrics from caret
  cm.lr <- confusionMatrix(labels, yhat.lr, positive = 'Class.2')
  
  # Store cumulative ROC curve
  yhat.lr <- as.numeric(yhat.lr) - 1
  auc.lr[i] <- pROC::auc(labels.n, yhat.lr, quiet=TRUE, ci=TRUE, stratified=TRUE, direction='<', plot=TRUE, legacy.axes=T)
  
  # Store the model performance metrics in the cumulative data frame
  metrics.lr[,i] <- cm.lr$byClass
}
metrics.lr <- data.frame("Logistic Regression"=rowMeans(metrics.lr, na.rm=TRUE))
row.names(metrics.lr) <- names(cm.lr$byClass)
metrics.lr[12,] <- mean(unlist(auc.lr))
row.names(metrics.lr)[12] <- 'AUC'

## Decision Tree
auc.tree <- numeric(times)
metrics.tree  <- matrix(nrow = 11, ncol = times)

# Perform 10-fold CV repeated 100 times
for (i in 1:times) {
  set.seed(666 + i)  # Set a different seed for each iteration
  # Create a matrix with the training indexes
  split_index <- createDataPartition(df$y, p = 0.7, list = FALSE)
  
  # Create the training set
  train <- df[split_index,]
  
  # Create the test set
  test <- df[-split_index,]
  
  # Store the test y in labels
  labels <- test$y
  labels.n <- as.numeric(labels) - 1
  
  # Train the Decision Tree
  tree <- train(
    y ~ .,
    train,
    method = "rpart2",
    tuneGrid = expand.grid(maxdepth=5),
    metric = "ROC",
    trControl = ctrl
  )
  
  yhat.tree <- predict(tree, test)
  
  # Create a confusion matrix and get performance metrics from caret
  cm.tree <- confusionMatrix(labels, yhat.tree, positive = 'Class.2')
  
  # Store cumulative ROC curve
  yhat.tree <- as.numeric(yhat.tree) - 1
  auc.tree[i] <- pROC::auc(labels.n, yhat.tree, quiet=TRUE)
  
  # Store the model performance metrics in the cumulative data frame
  metrics.tree[,i] <- cm.tree$byClass
}
metrics.tree <- data.frame('Decision Tree'=rowMeans(metrics.tree, na.rm=TRUE))
metrics.tree[12,] <- mean(unlist(auc.tree))
row.names(metrics.tree)[12] <- 'AUC'

## Random Forest
auc.rf <- list(times)
metrics.rf  <- matrix(nrow = 11, ncol = times)

# Perform 10-fold CV repeated 100 times
for (i in 1:times) {
  set.seed(666 + i)  # Set a different seed for each iteration
  
  # Create a matrix with the training indexes
  split_index <- createDataPartition(df$y, p = 0.7, list = FALSE)
  
  # Create the training set
  train <- df[split_index,]
  
  # Create the test set
  test <- df[-split_index,]
  
  # Store the test y in labels
  labels <- test$y
  labels.n <- as.numeric(labels) - 1
  
  # Train the logistic regression model
  rf <- train(
    y ~ .,
    train,
    method = "rf",
    tuneLength=5,
    metric = "ROC",
    trControl = ctrl, logit=TRUE
  )
  
  yhat.rf <- predict(rf, test)
  
  # Create a confusion matrix and get performance metrics from caret
  cm.rf <- confusionMatrix(labels, yhat.rf, positive = 'Class.2')
  
  # Store cumulative ROC curve
  yhat.rf <- as.numeric(yhat.rf) - 1
  auc.rf[i] <- pROC::auc(labels.n, yhat.rf,
                         # arguments for ci
                         ci=TRUE, direction='<', plot=TRUE, legacy.axes=T)
  
  # Store the model performance metrics in the cumulative data frame
  metrics.rf[,i] <- cm.rf$byClass
}
metrics.rf <- data.frame('Random Forest'=rowMeans(metrics.rf, na.rm=TRUE))
metrics.rf[12,] <- median(unlist(auc.rf))
row.names(metrics.rf)[12] <- 'AUC'
 

### Model comparison

#Create a data frame to store and compare the metrics
metrics <- cbind(metrics.lr, metrics.tree, metrics.rf)

#Display the table of metrics
metrics


fit <- function(X = NULL,
                Y = NULL,
                max_depth= 4,
                mtry= 6,
                seed= 123) {
  
  # set y as a factor
  y <- as.factor(x)
  
  #Encode y so valid R variables can be created from its classes
  levels(y) <- make.names(levels(y))
  
  ### Training-test split
  #Merge X and y into df
  df <- cbind.data.frame(X, y)
  set.seed(seed)
  
  #Create a matrix with the training indexes
  split_index <- createDataPartition(df$y, p = 0.7, list = FALSE)
  
  #Create the training set
  train <- df[split_index,]
  
  #Create the test set
  test <- df[-split_index,]
  
  #Store the test y in labels
  labels <- test$y
  labels.n <- as.numeric(labels) - 1
  
  #Define the training control hyperparameters
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 100,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  ### Logistic Regression
  set.seed(seed)
  
  lr <- train(y ~ .,
              train,
              method = "glm",
              metric = "ROC",
              trControl = ctrl)
  yhat.lr <- predict(lr, test)
  
  #Create a confusion matrix and get performance metrics from caret
  cm.lr <- confusionMatrix(labels, yhat.lr, positive = 'Class.2')
  
  #Plot ROC curve for LR
  yhat.lr <- as.numeric(yhat.lr) -1
  pROC::roc(labels.n, yhat.lr,
            plot=TRUE)
  
  #Create a data frame to store the performance metrics
  metrics.lr <<- data.frame(cm.lr$byClass)
  auc.lr <<- perf.lr@y.values[[1]][2]
  colnames(metrics.lr) <- 'Logistic Regression'
  
  ### Decision Tree
  set.seed(seed)
  # Define the tree training hyperparameters
  tree <- train(
    y ~ .,
    train,
    method = "rpart2",
    tuneLength = max_depth,
    metric = "ROC",
    trControl = ctrl
  )
  #Predict
  yhat.tree <- predict(tree, test)
  
  #Create a confusion matrix and get metrics
  cm.tree <- confusionMatrix(labels, yhat.tree, positive = 'Class.2')
  
  #Create a ROC curve
  yhat.tree <- as.numeric(yhat.tree) - 1
  pred.tree <- prediction(yhat.tree, labels.n)
  perf.tree <<- performance(pred.tree, "tpr", "fpr")
  plot(perf.tree,
       col=2,
       lwd = 2,
       add = TRUE
  )
  
  #Store metrics in a new dataframe
  metrics.tree <<- data.frame(cm.tree$byClass)
  auc.tree <<- perf.tree@y.values[[1]][2]
  colnames(metrics.tree) <- 'Decision Tree'
  
  
  ### Random Forest
  
  set.seed(seed)
  rf <- train(
    y ~ .,
    train,
    method = "rf",
    tuneLength = mtry,
    metric = "ROC",
    trControl = ctrl
  )
  #Predict
  yhat.rf <- predict(rf, test)
  
  #Get the confusion matrix and other metrics
  cm.rf <- confusionMatrix(labels,yhat.rf, positive='Class.2')
  
  #Create a ROC curve
  yhat.rf <- as.numeric(yhat.rf) -1
  pred.rf <- prediction(yhat.rf, labels.n)
  perf.rf <- performance(pred.rf, "tpr", "fpr")
  plot(
    perf.rf,
    col= 3,
    lwd = 2,
    add= TRUE
  )
  legend("bottomright", c("LR", "DT","RF"), lwd=2, 
         col = c(1:3), bty="n")
  #Store metrics in a new dataframe
  metrics.rf <- data.frame(cm.rf$byClass)  
  #Add AUC
  auc.rf <- perf.rf@y.values[[1]][2]
  colnames(metrics.rf) <- 'Random Forest'
  ### Model comparison
  
  #Create a data frame to store and compare the metrics
  metrics <<- cbind(metrics.lr, metrics.tree, metrics.rf)
  metrics[12,] <- auc <<- c(auc.lr, auc.tree, auc.rf)
  row.names(metrics)[12] <- 'AUC' 
  
  #Display the table of metrics
  return(metrics)
}
