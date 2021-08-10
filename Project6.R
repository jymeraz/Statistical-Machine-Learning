################
## Data cleaning
################

ILPD <- read.csv(file="ILPD.csv",
                 header=FALSE, col.names=c("age", "gender", "TB", "DB", "alkphos",
                                           "sgpt", "sgot", "TP", "alb", "AGratio", "liver"))
dim(ILPD); head(ILPD)

table(ILPD$liver)
length(which((ILPD$liver == 1)))/nrow(ILPD)

# check for missing data
miss.rate <- sapply(ILPD, function(x)sum(is.na(x)/nrow(ILPD)))
format <- paste(round((miss.rate)*100,digits=2),"%",sep="")
miss.rate <- structure(format, names=colnames(ILPD))
miss.rate

suppressPackageStartupMessages(library(VIM))
dat <- kNN(ILPD, variable=c("AGratio"), k=5)
dat <- subset(dat, select=-c(AGratio_imp))

dat$liver <- ifelse(dat$liver==1, 1, 0)

################
## EDA
################

sapply(dat, typeof)

library(car)
vars.nominal <- c("gender")
cols.x <- 1:(NCOL(dat)-1)

xnames <- names(dat)[cols.x]
y <- dat$liver
OUT <- NULL
for (j in 1:length(cols.x)){
  x <- dat[, cols.x[j]]
  xname <- xnames[j]
  if (is.element(xname, vars.nominal)){
    tbl <- table(x, y)
    # chi square test
    pvalue <- chisq.test(tbl)$p.value
  } else {
    # two-sample t test
    pvalue.equal.var <- (leveneTest(x~factor(y))$"Pr(>F)")[1]
    equal.var <- ifelse(pvalue.equal.var <= 0.05, FALSE, TRUE)
    pvalue <- t.test(x~y, alternative="two.sided",
                     var.equal=equal.var)$p.value
  }
  OUT <- rbind(OUT, cbind(xname=xname, pvalue=pvalue))
}
OUT <- as.data.frame(OUT)
colnames(OUT) <- c("name", "pvalue")
OUT

# exclude predictors that are associated with a p-value larger than 0.20
dat <- subset(dat, select=-c(TP))

################
## Variable selection
################

formula0 <- liver ~ age + factor(gender) + TB + DB + alkphos + sgpt + sgot + alb + AGratio
fit.full <- glm(formula0, data=dat)
summary(fit.full)
BIC(fit.full)

library(MASS)
fit.step <- step(fit.full, direction = "both", k=log(nrow(dat)))
# Refit the model
fit.step <- glm(liver ~ age + DB + alkphos + sgpt, family=binomial, data=dat)
summary(fit.step)
BIC(fit.step)

X <- model.matrix(as.formula(formula0),data=dat)
Y <- dat$liver

library(ncvreg)
cvfit.SCAD <- cv.ncvreg(X=X,y=y, nfolds=5, family="binomial", penalty="SCAD", 
                        lambda.min=.01, nlambda=9, eps=.01, max.iter=1000) 

# Plot the variables selected
plot(cvfit.SCAD)
result.SCAD <- cvfit.SCAD$fit
beta.hat <- as.vector(result.SCAD$beta[-1, cvfit.SCAD$min])
cutoff <- 0
terms <- colnames(X)[abs(beta.hat) > cutoff]

# Take into consideration the factoring for categorical variables
for (i in 1:length(terms)) {
  if(terms[i]=="factor(gender)Male") {
    terms[i] = "gender"
  }
}

# Refit the model
formula.SCAD <- as.formula(paste(c("liver ~ 1", terms), collapse=" + "))
fit.pen <- glm(formula.SCAD, data = dat, family="binomial")
summary(fit.pen)
BIC(fit.pen)

################
## Comparing models
################

library(caret)
n <- NROW(dat)

# Full model
full.p.jk <- rep(0, n)
for (i in 1:n){
  fit.i <- glm(formula(fit.full), data=dat[-i,], family = "binomial")
  full.p.jk[i] <- predict(fit.i, newdata=dat[i,], type="response")
}
confusionMatrix(factor(sign(full.p.jk >= 0.5)), factor(y))

# Stepwise model
step.p.jk <- rep(0, n)
for (i in 1:n){
  fit.i <- glm(formula(fit.step), data=dat[-i,], family = "binomial")
  step.p.jk[i] <- predict(fit.i, newdata=dat[i,], type="response")
}
confusionMatrix(factor(sign(step.p.jk >= 0.5)), factor(y))


# Penalization model
pen.p.jk <- rep(0, n)
for (i in 1:n){
  fit.i <- glm(formula(fit.pen), data=dat[-i,], family = "binomial")
  pen.p.jk[i] <- predict(fit.i, newdata=dat[i,], type="response")
}
confusionMatrix(factor(sign(pen.p.jk >= 0.5)), factor(y))

p.jk <- c(full.p.jk, step.p.jk, pen.p.jk)

library(verification)
library(cvAUC)

# Full model
y <- dat$liver
yhat <- full.p.jk
a.ROC <- roc.area(obs=y, pred=yhat)$A
AUC <- ci.cvAUC(predictions=yhat, labels=y, folds=1:n, confidence=0.95)
mod.glm <- verify(obs=y, pred=yhat)
roc.plot(mod.glm, plot.thres = NULL, main="Full Model")
text(x=0.7, y=0.2, paste("AUC = ", round(AUC$cvAUC, digits=3), 
                         " with 95% CI (", round(AUC$ci, digits=3)[1], ",", round(AUC$ci, digits=3)[2], ").",
                         sep=""), col="blue", cex=1.2)

# Stepwise model
y <- dat$liver
yhat <- step.p.jk
a.ROC <- roc.area(obs=y, pred=yhat)$A
AUC <- ci.cvAUC(predictions=yhat, labels=y, folds=1:n, confidence=0.95)
mod.glm <- verify(obs=y, pred=yhat)
roc.plot(mod.glm, plot.thres = NULL, main="Stepwise Model")
text(x=0.7, y=0.2, paste("AUC = ", round(AUC$cvAUC, digits=3), 
                         " with 95% CI (", round(AUC$ci, digits=3)[1], ",", round(AUC$ci, digits=3)[2], ").",
                         sep=""), col="blue", cex=1.2)

# Penalization model
y <- dat$liver
yhat <- pen.p.jk
a.ROC <- roc.area(obs=y, pred=yhat)$A
AUC <- ci.cvAUC(predictions=yhat, labels=y, folds=1:n, confidence=0.95)
mod.glm <- verify(obs=y, pred=yhat)
roc.plot(mod.glm, plot.thres = NULL, main="Penalty Model")
text(x=0.7, y=0.2, paste("AUC = ", round(AUC$cvAUC, digits=3), 
                         " with 95% CI (", round(AUC$ci, digits=3)[1], ",", round(AUC$ci, digits=3)[2], ").",
                         sep=""), col="blue", cex=1.2)

################
## Confidence intervals
################

library(MASS)
ci <- confint(fit.step, level = 0.95) 
ci
exp(ci)

