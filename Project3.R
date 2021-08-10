################
## Loading the data
################

# Read data into R.
dat <- read.csv("hmeq.csv")

# Compute and print the missing rate of each variable.
miss.rate <- sapply(dat, function(x)sum(is.na(x)/nrow(dat)))
format <- paste(round((miss.rate)*100,digits=2),"%",sep="")

cols <- colnames(dat)
miss.rate <- structure(format, names=cols)
miss.rate

################
## Data Cleaning
################

# Replace missing values for JOB.
replaceMissing <- as.vector(dat$JOB)
replaceMissing[is.na(replaceMissing) | replaceMissing=="" | replaceMissing=="NA"] <- "Unknown"
dat$JOB <- replaceMissing
table(dat$JOB, useNA="ifany")

# Replace missing values for REASON.
replaceMissing <- as.vector(dat$REASON)
replaceMissing[is.na(replaceMissing) | replaceMissing=="" | replaceMissing=="NA"] <- "Unknown"
dat$REASON <- replaceMissing
table(dat$REASON, useNA="ifany")

# log transformation
dat$LOAN <- ifelse(dat$LOAN == 0, log(dat$LOAN + 1), log(dat$LOAN))
dat$VALUE <- ifelse(dat$VALUE == 0, log(dat$VALUE + 1), log(dat$VALUE))
dat$MORTDUE <- ifelse(dat$MORTDUE == 0, log(dat$MORTDUE + 1), log(dat$MORTDUE))
dat$YOJ <- ifelse(dat$YOJ == 0, log(dat$YOJ + 1), log(dat$YOJ))
dat$CLAGE <- ifelse(dat$CLAGE == 0, log(dat$CLAGE + 1), log(dat$CLAGE))

# Impute remaining values
suppressPackageStartupMessages(library(VIM))

# Missing values imputed with the k-nearest-neighbor imputation technique.
dat <- kNN(dat, variable=c("MORTDUE", "VALUE", "YOJ", "DEROG", "DELINQ", "CLAGE", "NINQ", 
                           "CLNO", "DEBTINC"), k=5)

# Obtain a distance matrix
# Exclude the variable BAD.
dat1 <- dat[,-1]

# Obtain the distance matrix.
library(cluster)
numericDat <- model.matrix(~.-1, data=dat1)
distanceMatrix <- daisy(numericDat, metric = "gower", stand = F) 

################
## Clustering Algorithms
################

# Hierarchical Clustering
fit.ward <- hclust(distanceMatrix, method="complete") 
plot(fit.ward , hang = -0.5) 

# Screeplot to determine the number of clusters
K.max <- 30
height <- tail(fit.ward$height, n=K.max)
n.cluster <- tail((nrow(numericDat)-1):1, n=K.max)
plot(n.cluster, height,  type="b", pch=19, cex=.5, xlab="number of clusters", 
     ylab="height", col="blue", lwd=2)
axis(1, seq(0,30,1))

# Cut tree into 2 clusters and plot the dendogram with the two clusters.
# Output the cluster memberships.
groups <- cutree(fit.ward, k=2) 
plot(fit.ward, hang = -0.5)   
rect.hclust(fit.ward, k=2, border="red") 

# Plot the data using PCA.
library(ggplot2)
library(factoextra)
fviz_cluster(list(data = numericDat, cluster = groups), geom="point") + geom_text(aes(label=dat$BAD))

# Screeplot to determine the number of clusters in k-means clustering
wss <- (nrow(numericDat)-1)*sum(apply(numericDat,2,var))
K.max <- 15
for (K in 2:K.max) wss[K] <- sum(kmeans(numericDat, centers=K)$withinss)
plot(1:K.max, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 

# 3 clusters
K <- 3

# K-Means Cluster Analysis
fit.k <- kmeans(numericDat, K) 

# Output the cluster memberships.
aggregate(numericDat,by=list(fit.k$cluster),FUN=mean)

# Plot the data using PCA.
library(ggplot2)
library(factoextra)
fviz_cluster(list(data = numericDat, cluster = fit.k$cluster), geom="point") + geom_text(aes(label=dat$BAD))

# Compute contingency table with clusters from hierarchical and k-means clustering.
contingencyTable <- table(groups, fit.k$cluster)
contingencyTable
chisq.test(contingencyTable)

# Compute and print jaccard and rand similarity index
library(jaccard)
jaccard(groups, fit.k$cluster)

library(fossil)
rand.index(groups, fit.k$cluster)

################
## Post-hoc analysis
################

# Print a summary of the groups to see where the groups are skewed towards.
summary(groups)

# Use the chi-square test and the fisher test to check for associations between the variables.
cols <- colnames(dat)
for (i in 1:length(cols)) {
  variable = cols[i]
  if (typeof(dat[,variable])=="integer" | typeof(dat[,variable])=="double") {
    tab <- table(dat[,variable], groups, useNA="no")
    par(mfrow=c(1,1))
    boxplot(dat[,variable] ~ groups, col = "orange",
            xlab=groups,
            main=paste("Parallel Boxplots of ", variable, sep=""))
  }
  print(variable)
  print(chisq.test(tab))
  print(fisher.test(tab, simulate.p.value =TRUE), end = "\n")
}








