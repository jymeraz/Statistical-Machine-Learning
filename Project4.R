################
## Loading the data
################

baseball <- read.table(file=
                         "http://www.amstat.org/publications/jse/datasets/baseball.dat.txt",
                       header = F, col.names=c("salary", "batting.avg", "OBP", "runs", "hits",
                                               "doubles", "triples", "homeruns", "RBI", "walks", "strike.outs",
                                               "stolen.bases", "errors", "free.agency.elig", "free.agent.91",
                                               "arb.elig", "arb.91", "name"))
head(baseball)

################
## Prepare the data
################

# Histogram of salary
par(mfrow=c(2,2),mar=c(4, 4, 4, 4))
hist(baseball$salary, xlab="salary", main="Histogram of Salary")
qqnorm(baseball$salary, main="Q-Q Plot of Salary")
qqline(baseball$salary)

# Histogram of log salary
logsalary <- log(baseball$salary)
hist(logsalary, xlab="log-alary", main="Histogram of Log(Salary)")
qqnorm(logsalary, main="Q-Q Plot of Salary")
qqline(logsalary)

# Proceed with the log-transformed salary
dat <- baseball[,-1]
dat$logsalary <- logsalary

# Display the types of each predictor
head(dat)

# Compute and print the missing rate of each variable
miss.rate <- sapply(dat, function(x)sum(is.na(x)/nrow(dat)))
format <- paste(round((miss.rate)*100,digits=2),"%",sep="")
miss.rate <- structure(format, names=colnames(dat))
miss.rate

################
## Linear Regression
################

# Partition data
n <- NROW(dat); ratio <- 2/3
set.seed(123)
id.training <- sample(1:n, size=trunc(n*ratio), replace=FALSE)
D <- dat[id.training, ] # training
D0 <- dat[-id.training, ] # test

fit.full <- lm(logsalary ~ batting.avg + OBP + runs + hits + doubles + 
                 triples + homeruns + RBI + walks + strike.outs + 
                 stolen.bases + errors + free.agency.elig + 
                 free.agent.91 + arb.elig + arb.91, data=D)

# ======================
# Stepwise Regression
# ======================
library(MASS)
fit.step <- stepAIC(fit.full, direction="backward", k=log(nrow(D))) # Stepwise fit

# ======================
# LASSO
# ======================
formula0 <- logsalary ~ batting.avg + OBP + runs + hits + doubles + 
  triples + homeruns + RBI + walks + strike.outs + stolen.bases + 
  errors + free.agency.elig + free.agent.91 + arb.elig + arb.91 -1		
y <- D[, all.vars(formula0)[1]]
X <- as.data.frame(model.matrix(as.formula(formula0),D))

library(ncvreg)
cvfit.L1 <- cv.ncvreg(X=X,y=y, nfolds=10, family="gaussian",
                      penalty="lasso", lambda.min=.005, nlambda=100, eps=.001, max.iter=1000)
beta.hat <- coef(cvfit.L1)  # LASSO coefficients

# Refit the model with the chosen variables.
cutoff <- 0.0001
terms <- names(beta.hat)[abs(beta.hat) > cutoff]
formula.LASSO <- as.formula(paste(c("logsalary ~ ", terms[-1]), collapse=" + "))
fit.L1 <- lm(formula.LASSO, data=D) # LASSO fit

# =========================
# SCAD 
# =========================
y <- D[, all.vars(formula0)[1]]
X <- model.matrix(as.formula(formula0),D)
X <- as.data.frame(scale(X, center = TRUE, scale = TRUE));   	# Standardize X
y <- scale(y, center = TRUE, scale = FALSE);	# Center y
dat1 <- as.data.frame(cbind(X, y))

library(ncvreg)
cvfit <- cv.ncvreg(X=X,y=y, nfolds=10,
                   family="gaussian", penalty="SCAD") 	
result.SCAD <- cvfit$fit
beta.hat <- as.vector(result.SCAD$beta[-1, cvfit$min])
cutoff <- 0
terms <- colnames(X)[abs(beta.hat) > cutoff]
formula.SCAD <- as.formula(paste(c("logsalary ~ 1", terms), collapse=" + "))
fit.SCAD <- lm(formula.SCAD, data=D) # SCAD fit

# Compute AIC values
require(kableExtra)
Techniques <- c("SCAD", "LASSO(L1)", "STEPWISE")
aic <- c(AIC(fit.SCAD), AIC(fit.L1), AIC(fit.step))
result <- data.frame(Techniques,aic)
result

summary(fit.L1)
summary(fit.step)
summary(fit.SCAD)

################
## Test the models
################

library(DAAG)
formula0 <- logsalary ~ batting.avg + OBP + runs + hits + doubles + triples + 
  homeruns + RBI + walks + strike.outs + stolen.bases + errors + 
  free.agency.elig + free.agent.91 + arb.elig + arb.91 -1			
y <- D0[, all.vars(formula0)[1]]
X <- model.matrix(as.formula(formula0),D0)
DAT <- data.frame(cbind(logsalary=y, X))

CV <- CVlm(data=D0, form.lm=fit.L1, m=10, seed=123, plotit=T, printit=T)
SSPE1 <- sum((CV$logsalary-CV$cvpred)^2)
CV <- CVlm(data=D0, form.lm=fit.step, m=10, seed=123, plotit=T, printit=T)
SSPE2 <- sum((CV$logsalary-CV$cvpred)^2)
CV <- CVlm(data=D0, form.lm=fit.SCAD, m=10, seed=123, plotit=T, printit=T)
SSPE3 <- sum((CV$logsalary-CV$cvpred)^2)

c(SSPE1, SSPE2, SSPE3)

################
## Model diagnostics on final model
################

fit.final <- lm(formula=formula.LASSO, data=dat)
summary(fit.final)

####
## normality
####
library(car)
# Obtain studentized jackknife residuals
r.jack <- rstudent(fit.final)

# histogram
par(mfrow=c(1,2),mar=c(8,4,8,4)) 
hist(r.jack, xlab="Jackknife Residual", col="green4",
     main="(a) Histogram") 

# q-q plot
qqPlot(fit.final, pch=19, cex=.8, col="blue", main="(b) Q-Q Plot") 

# Shapiro-Wilk test for normality
shapiro.test(r.jack) 

####
## Homoscedasticity
####
# Plot Absolute Jackknife Residuals vs. Fitted values 
par(mfrow=c(1,1),mar=c(4, 4, 4, 4)) 
spreadLevelPlot(fit.final, pch=18, cex=0.5, col="blue",
                main="HV Model on Baseball Salary: Heteroscedasticity")

# The Breusch-Pagan Test for Non-Constant Error Variance 
ncvTest(fit.final) 

####
## Independence
####
# Durbin-Watson test
durbinWatsonTest(fit.final)

####
## Linearity
####
# Partial regression plot
leveragePlots(fit.final, main="Partial Regression (Leverage) Plots")

####
## Outlier detection
####
# Number of observations and slopes
n <- NROW(dat); 
p <- length(coef(fit.final))-1

# Using influencePlot
influencePlot(fit.final, id.method="identify", 
              col="blue", 
              main="Influence Plot", 
              sub="Circle size is proportial to Cook's d")

# Computing outliers
infl <- influence.measures(fit.final); 
infl.mat <- as.data.frame(infl$infmat)
h <- infl.mat$hat # leverage hii
cook.d <- infl.mat$cook.d # Cook’s distance di

# Bubble plot of leverage hii, studentized jackknife residuals ri, and Cook’s distance di
par(bg="white", mar=c(5, 5, 5, 5), mfrow=c(1, 1), xaxt="s")
plot(x=c(min(h), max(h)), y=c(min(r.jack), max(r.jack)), xlab=expression(h[ii]), 
     ylab=expression(r[(-i)]), cex.lab=1.5, type="n")
symbols(h, r.jack, circles=cook.d, inches=0.35, fg="white", bg="green", add=T)
abline(h=0, col="black", lwd=1)
abline(h=qt(.975, n-p-2), col="blue", lwd=2)
abline(h=qt(.025, n-p-2), col="blue", lwd=2)
text(x=.12, y=qt(.975, n-p-2)+.3, labels=expression(t[.975]^(n-p-2)),  col="blue")
text(x=.12, y=qt(.025, n-p-2)+.3, labels=expression(t[.025]^(n-p-2)),  col="blue")
abline(v=2*(p+1)/n, lwd=4, col="grey")
text(2*(p+1)/n+.008, 3.0, labels="2(p+1)/n", col="grey")

# Identifying and plotting outliers
t0 <- qt(.975, n-1-(p+1))
subset1 <- (cook.d >=0.065)|(h> 2*(p+1)/n)|(abs(r.jack) > t0) 
symbols(h[subset1], r.jack[subset1], circles=cook.d[subset1], inches=0.35,
        fg="white", bg="orangered", add=T)
text(h[subset1], r.jack[subset1], (1:n)[subset1], cex=0.5)

# Output of the outliers
cbind(id=(1:n)[subset1],  "r.jack"= abs(r.jack[subset1]) > t0, "h"= h[subset1]> 2*(p+1)/n, 
      "cook.d" = cook.d[subset1] >= 0.065)

####
## Multicollinearity
####
# Condition number
fit <- lm(logsalary ~ batting.avg + runs + RBI + walks + strike.outs + stolen.bases + 
            errors + free.agency.elig + free.agent.91 + arb.elig + arb.91,  
          data=dat, x=TRUE)
kappa(fit$x)

# Variance inflation factor
vif(fit)

################
## Final model prediction
################

library(ggplot2)

test <- read.csv("bb92-test.csv")
pred <- predict(fit.final, test, interval="prediction");
dat.plot <- data.frame(player=1:20, exp(pred)); names(dat.plot)
ggplot(dat.plot, aes(x=player, y=fit)) +
  geom_errorbar(aes(ymin=lwr, ymax=upr)) + geom_point()




