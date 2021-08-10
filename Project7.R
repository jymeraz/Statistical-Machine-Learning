################
## Data cleaning
################

library(kernlab); data(spam)
miss.rate <- sapply(spam, function(x)sum(is.na(x)/nrow(spam)))
format <- paste(round((miss.rate)*100,digits=2),"%",sep="")
miss.rate <- structure(format, names=colnames(spam))
miss.rate

# Partition data and set a seed
n <- nrow(spam); ratio <- 2/3
set.seed(123)
id.training <- sample(1:n, size=trunc(n*ratio), replace=FALSE)
train <- spam[id.training, ] # training
test <- spam[-id.training, ] # test

################
## Supervised learning
################

################
## Linear discriminant analysis (LDA)
################

library(MASS)
library(verification)
library(cvAUC)

# Fitting the model
fit.LDA <- lda(type ~ ., data=train)

# Prediction
yhat.LDA <- predict(fit.LDA, newdata=test)$x
yobs <- ifelse(test$type=="spam", 1 ,0)

# ROC curve and AUC
n <- NROW(test)
AUC.LDA <- ci.cvAUC(predictions=yhat.LDA, labels=yobs, folds=1:n, confidence=0.95)
yhat <- scale(yhat.LDA, center=min(yhat.LDA), scale = max(yhat.LDA)-min(yhat.LDA))
mod.glm <- verify(obs=yobs, pred=yhat)
roc.plot(mod.glm, plot.thres = NULL, main="ROC Curve from LDA Fitting")
text(x=0.5, y=0.2, paste("AUC = ", round(AUC.LDA$cvAUC, digits=3), 
                         " with 95% CI (", round(AUC.LDA$ci, digits=3)[1], ",", round(AUC.LDA$ci, digits=3)[2], ").",
                         sep=""), col="blue", cex=1.2)

################
## logistic regression model
################

library(glmnet)
X <- model.matrix(type~.,data=train)
y <- train$type
CV <- cv.glmnet(x=X, y=y, family="binomial", alpha = 1, 
                lambda.min = 1e-4, nlambda = 57, standardize = T, thresh = 1e-07, 
                maxit=1000)

# Selecting the best tuning parameter
b.lambda <- CV$lambda.1se
fit.lasso <- glmnet(x=X, y=y, family="binomial", alpha = 1, 
                    lambda=b.lambda, standardize = T, thresh = 1e-07, 
                    maxit=1000)

# Prediction
X.test <- model.matrix(type~., data=test)
yhat <- predict(fit.lasso, newx=X.test, s=b.lambda, type="response")
y <- ifelse(test$type=="spam",1,0)

# ROC curve and AUC
AUC <- ci.cvAUC(predictions=yhat, labels=y, folds=1:n, confidence=0.95)
mod.glm <- verify(obs=y, pred=yhat)
roc.plot(mod.glm, plot.thres = NULL, main="ROC Curve from Logistic Regression")
text(x=0.7, y=0.2, paste("AUC = ", round(AUC$cvAUC, digits=3),
                         " with 95% CI (", round(AUC$ci, digits=3)[1], ",", round(AUC$ci, digits=3)[2], ").",
                         sep=""), col="blue", cex=1.2)

################
## One single decision tree
################

library(rpart)
library(verification)

# Fitting the model
tre0 <- rpart(type~., data=train,  method="class")
cv.error <- (tre0$cptable)[,4]
a0 <- 1 
SE1 <- min(cv.error) + a0*((tre0$cptable)[,5])[which.min(cv.error)] 
position <- min((1:length(cv.error))[cv.error <= SE1])

# Finding the best tree
best.cp <-  sqrt(tre0$cptable[position,1] *  tre0$cptable[(position-1),1])
best.tree <- prune(tre0, cp=best.cp)

# Prediction
yhat <- predict(best.tree, newdata=test, type="prob")[,2]
y <- ifelse(test$type=="spam",1,0)

# ROC curve and AUC
a.ROC <- roc.area(obs=y, pred=yhat)$A 
AUC <- round(a.ROC, digits=4)
par(mfrow=c(1,1), mar=c(4, 4, 4, 4))
mod.glm <- verify(obs=y, pred=yhat, bins = FALSE)
roc.plot(mod.glm, plot.thres = NULL, main ="ROC Curve from One single decision tree")
text(x=0.6, y=0.3, paste("Area under ROC = ", AUC, sep=""))

################
## Bagging 
################

library(ipred)
library(verification)

# Fit the model
fit.bagging <- bagging(type~., data=train, nbagg=50, coob=TRUE)

# Prediction
yhat <- predict(fit.bagging, newdata=test, type="prob")[, 2]
y <- ifelse(test$type=="spam",1,0)

# ROC curve and AUC
AUC <- roc.area(obs=y, pred=yhat)$A
mod.glm <- verify(obs=y, pred=yhat)
roc.plot(mod.glm, plot.thres = NULL, col="red", main ="ROC Curve from Bagging")
text(x=0.7, y=0.2, paste("Area under ROC =", round(AUC, digits=4), 
                         sep=" "), col="blue", cex=1.2)

################
## Random Forests (RF)
################

library(randomForest)
library(verification)

# Fit the model
fit.rf <- randomForest(type~., data=train, importance=TRUE, proximity=TRUE, ntree=500)

# Prediction
yhat <- predict(fit.rf, newdata=test, type="prob")[, 2]
y <- ifelse(test$type=="spam",1,0)

# ROC curve and AUC
AUC <- roc.area(obs=y, pred=yhat)$A
mod.glm <- verify(obs=y, pred=yhat)
roc.plot(mod.glm, plot.thres = NULL, col="red", main ="ROC Curve from Random Forests")
text(x=0.7, y=0.2, paste("Area under ROC =", round(AUC, digits=4), 
                         sep=" "), col="blue", cex=1.2)

################
## Boosting
################

library(ada)
library(verification)

# Fit model
stump <- rpart.control(cp=-1, maxdepth=1, minsplit=0)
fit.stump <- ada(type~., data=train, iter=2000, loss="e", type="discrete", control=stump)

# Prediction
yhat <- predict(fit.stump, newdata=test, type="probs")[, 2]
y <- ifelse(test$type=="spam",1,0)

# ROC curve and AUC
AUC <- roc.area(obs=y, pred=yhat)$A
mod.glm <- verify(obs=y, pred=yhat)
roc.plot(mod.glm, plot.thres = NULL, col="red", main ="ROC Curve from Boosting")
text(x=0.7, y=0.2, paste("Area under ROC =", round(AUC, digits=4),
                         sep=" "), col="blue", cex=1.2)

################
## Additional RF features
################

library(randomForest)
library(verification)

# Fit model
fit.RF <- randomForest(type~., data=spam, importance=TRUE, proximity=TRUE, ntree=2000)

# Plot of variable importance ranking
varImpPlot(fit.RF, main="Variable Importance")

# Variables in logistic regression
fit.lasso$beta

# Partial dependence plot 
par(mfrow=c(1, 3), mar=rep(4,4))
partialPlot(fit.RF, pred.data=spam, x.var=charExclamation, rug=TRUE)
partialPlot(fit.RF, pred.data=spam, x.var=remove, rug=TRUE)

fit.mds <- cmdscale(1 - fit.RF$proximity, eig=TRUE) 
plot(fit.mds$points, xlab="dim I", ylab="dim II", pch=19, cex=0.8, 
     col=c("red", "green")[as.numeric(spam$type)])










