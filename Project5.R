################
## Prepare the data
################

dat <- read.csv("crime.csv")
dat0 <- dat[,-c(1:5)]

# Compute and print the missing rate of each variable
miss.rate <- sapply(dat0, function(x)sum(is.na(x)/nrow(dat0)))
format <- paste(round((miss.rate)*100,digits=2),"%",sep="")
miss.rate <- structure(format, names=colnames(dat0))
miss.rate

# Remove the variables that are missing over 60% of data
dat0 <- subset(dat0, select=-c(LemasSwornFT,LemasSwFTPerPop, LemasSwFTFieldOps, 
                               LemasSwFTFieldPerPop, LemasTotalReq, LemasTotReqPerPop, 
                               PolicReqPerOffic, PolicPerPop, RacialMatchCommPol, 
                               PctPolicWhite, PctPolicBlack, PctPolicHisp, PctPolicAsian, 
                               PctPolicMinor, OfficAssgnDrugUnits, NumKindsDrugsSeiz, 
                               PolicAveOTWorked, PolicCars, PolicOperBudg, 
                               LemasPctPolicOnPatr, LemasGangUnitDeploy, PolicBudgPerPop))

suppressPackageStartupMessages(library(VIM))

# Missing values imputed with the k-nearest-neighbor imputation technique
#dat0$OtherPerCap_imp <- ifelse(dat0$OtherPerCap==FALSE, 0, 1)
dat0 <- kNN(dat0, variable="OtherPerCap", k=5)
table(dat0$OtherPerCap_imp)
dat0 <- subset(dat0, select=-c(OtherPerCap_imp))

################
## EDA
################

# Histogram
par(mfrow=c(1,2),mar=c(8,4,8,4))
hist(dat0$ViolentCrimesPerPop, xlab="ViolentCrimesPerPop", col="green4",
     main="Histogram")

# q-q plot
qqnorm(dat0$ViolentCrimesPerPop, pch = 1, frame = FALSE)
qqline(dat0$ViolentCrimesPerPop, col = "steelblue", lwd = 2)

# Perform a log-transormation
dat0$ViolentCrimesPerPop <- ifelse(dat0$ViolentCrimesPerPop == 0, log(dat0$ViolentCrimesPerPop + 1), log(dat0$ViolentCrimesPerPop))

# Histogram
par(mfrow=c(1,2),mar=c(8,4,8,4))
hist(dat0$ViolentCrimesPerPop, xlab="ViolentCrimesPerPop", col="green4",
     main="Histogram")

# q-q plot
qqnorm(dat0$ViolentCrimesPerPop, pch = 1, frame = FALSE)
qqline(dat0$ViolentCrimesPerPop, col = "steelblue", lwd = 2)

# Partition data and set a seed
n <- nrow(dat0); ratio <- 2/3
set.seed(123)
id.training <- sample(1:n, size=trunc(n*ratio), replace=FALSE)
D1 <- dat0[id.training, ] # training
D2 <- dat0[-id.training, ] # test

################
## Predictive modeling
################

################
## Principal Components Regression (PCR)
################

# ======================================
# Fitting the model
# ======================================
library(pls)
fit.PCR <- pcr(ViolentCrimesPerPop ~ ., ncomp=100, data=D1, method = pls.options()$pcralg, 
               validation = "CV", segments = 10, segment.type ="random", scale=TRUE)
CV <- fit.PCR$validation

par(mfrow=c(2,1), mar=rep(4,4))
plot(1:CV$ncomp, CV$PRESS,  xlab="number of PCs", ylab="PRESS", type="b", col="blue",
     lwd=2)

ncomp.best <- which.min(CV$PRESS)
fit.PCR.best <- pcr(ViolentCrimesPerPop ~ ., ncomp=ncomp.best, data=D1, 
                    method = pls.options()$pcralg, segments = 1, scale=TRUE)

# ======================================
# Prediction
# ======================================
yhat.PCR <- predict(fit.PCR.best, newdata=D2, comps=1:ncomp.best)
yobs <- D2$ViolentCrimesPerPop

# Predicted vs. Observed
par(mfrow=c(1,1), mar=rep(4,4))
plot(yobs, yhat.PCR, type="p", pch=18, col="blue", xlab="observed", ylab="predicted", 
     main="PCR")
abline(a=0, b=1, col="orange", lwd=2)

# MSE
MSEP.PCR <- mean((yobs-yhat.PCR)^2)

################
## Partial least squares regression (PLSR)
################

# ======================================
# Fitting the model
# ======================================
library(pls)
fit.PLS <- plsr(ViolentCrimesPerPop ~ ., ncomp=100, data=D1, method = "simpls", 
                validation = "CV", segments = 10, segment.type ="random", scale=TRUE)
CV <- fit.PLS$validation

par(mfrow=c(2,1), mar=rep(4,4))
plot(1:CV$ncomp, CV$PRESS,  xlab="number of PCs", ylab="PRESS", type="b", col="blue", 
     lwd=2)
plot(fit.PLS)
validationplot(fit.PLS, val.type = "MSEP", main="mean squared error of prediction")

ncomp.best <- which.min(CV$PRESS)
fit.PLSR.best <- plsr(ViolentCrimesPerPop ~ ., ncomp=ncomp.best, data=D1, 
                      method = "simpls", validation = "none", scale=F)

# ======================================
# Prediction
# ======================================
yhat.PLSR <- predict(fit.PLSR.best, newdata=D2, comps=1:ncomp.best) 
yobs <- D2$ViolentCrimesPerPop

# Predicted vs. Observed
par(mfrow=c(1,1), mar=rep(4,4))
plot(yobs, yhat.PLSR, type="p", pch=18, col="blue", xlab="observed", ylab="predicted", 
     main="Partial LS")
abline(a=0, b=1, col="orange", lwd=2)

# MSE
MSEP.PLSR <- mean((yobs-yhat.PLSR)^2)

################
## Weighted orthogonal components regression (WOCR)
################

# ======================================
# Fitting the model
# ======================================
library(WOCR)
fit.WOCR <- WOCR(ViolentCrimesPerPop~., data = D1, scale = TRUE, 
                 model = "PCR.gamma.a.c")

# ======================================
# Prediction
# ======================================
yhat.WOCR <- predict(fit.WOCR, newdata=D2)
MSEP.WOCR <- mean((D2$ViolentCrimesPerPop-yhat.WOCR)^2)

################
## Stagewise regression
################

# ======================================
# Fitting the model
# ======================================
library(lars)
# Identify ViolentCrimesPerPop as the dependent variable
X <- as.matrix(D1[, -101]); y <- D1[, 101]
fit.stagewise <- lars(X,y,type="for", trace = T, normalize = TRUE, intercept = TRUE)

par(mfrow=c(3,1), mar=rep(4,4))
plot(fit.stagewise, xvar="norm", breaks=F, plottype="coefficients",omit.zeros = T)
plot(fit.stagewise, xvar="norm", breaks=T, plottype="coefficients", omit.zeros = F)
plot(fit.stagewise, xvar= "step", breaks = T, plottype = "Cp")

# ======================================
# Prediction
# ======================================
b.max.steps <- 20
yobs <- D2$ViolentCrimesPerPop
yhat.stagewise <- predict(fit.stagewise,  s=b.max.steps, newx=as.matrix(D2[, -101]))$fit
MSEP.stagewise <- mean((yobs-yhat.stagewise)^2)

################
## Least angle regression (LAR)
################

# ======================================
# Fitting the model
# ======================================
# Identify ViolentCrimesPerPop as the dependent variable
X <- as.matrix(D1[, -101]); y <- D1[, 101]
fit.lar <- lars(X,y,type="lar", trace = FALSE, normalize = TRUE, intercept = TRUE) 

# Finding L1-norm
BETA <- as.matrix(fit.lar$beta)
L1.norm <- function(x) sum(abs(x))
norm.L1 <- apply(BETA, 1, L1.norm)
norm.L1 <- norm.L1/max(norm.L1)
df <- as.vector(unlist(fit.lar$df))
Cp <- as.vector(fit.lar$Cp)
RSS <- as.vector(fit.lar$RSS)
lambda <- c(as.vector(fit.lar$lambda), 0)
b.lambda <- lambda[Cp==min(Cp)]
b.L1norm <- norm.L1[Cp==min(Cp)]

par(mfrow=c(3,1), mar=rep(4,4))
plot(fit.lar, xvar="norm", breaks=F, plottype="coefficients",omit.zeros = T)
abline(v=b.L1norm, col="gold", lwd=3)
plot(fit.lar, xvar="norm", breaks=T, plottype="coefficients", omit.zeros = F)
plot(fit.lar, xvar= "df", plottype="Cp", omit.zeros = F)

# ======================================
# Prediction
# ======================================
b.max.steps <- 10
yhat.lar <- predict(fit.lar,  s=b.max.steps, newx=as.matrix(D2[, -101]))$fit
MSEP.lar <- mean((yobs-yhat.lar)^2)

################
## Compare MSEP values across models
################

rbind(MSEP.PCR, MSEP.PLSR, MSEP.WOCR, MSEP.stagewise, MSEP.lar)

