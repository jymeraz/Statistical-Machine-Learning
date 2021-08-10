################
## Loading the data
################

# Load the data
train <- read.table(file=
                      "http://archive.ics.uci.edu/ml/machine-learning-databases/optdigits/optdigits.tra",
                    sep=",", header = FALSE, na.strings = c("NA", "", " "),
                    col.names = c(paste("x", 1:64, sep=""), "digit"))
test <- read.table(file=
                     "http://archive.ics.uci.edu/ml/machine-learning-databases/optdigits/optdigits.tes",
                   sep=",", header = FALSE, na.strings = c("NA", "", " "),
                   col.names = c(paste("x", 1:64, sep=""), "digit"))
dim(train); dim(test)

dat <- rbind(train, test) 
dim(dat)

################
## Exploratory Data Analysis (EDA) 
################

# Obtain a heat map of the data
dat0 <- data.matrix(dat[order(dat$digit), -65])
n <- NROW(dat0)
color <- rainbow(n, alpha = 0.8)
heatmap(dat0, col=color, scale="column", Rowv=NA, Colv=NA,
        labRow=FALSE, margins=c(4,4), xlab="Image Variables", ylab="Samples",
        main="Heatmap of Handwritten Digit Data")

# Print unique values for each column
size <- length(colnames(dat0))
for (i in 1:size) {
  if (length(unique(dat0[,i])) < 5) {
    cat("column:", i, ", unique values:", unique(dat0[,i]), "\n", sep = " ")
  }
}

# Remove the unary values from the dataset
dat0 <- dat0[,-c(1, 40)]

################
## Principal Components Analysis (PCA)  
################

# Running PCA and scaling the data.
pca.res <- prcomp(dat0, retx=TRUE)

# Screeplots.
screeplot(pca.res)
screeplot(pca.res, type="lines")

# Cumulative proportion of variation.
sd.pc <- pca.res$sdev 
var.pc <- sd.pc^2
prop.pc <- var.pc/sum(var.pc)

plot(cumsum(prop.pc), xlab = "Principal Component", col="blue",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b", pch=19)

# First two PC directions.
PC.directions <- pca.res$rotation 
a1.a2 <- pca.res$rotation[,1:2] 
a1.a2

# Plot PC1 vs PC2
plot(pca.res$x[,1:2], pch="", main="PC1 vs PC2")
text(pca.res$x[,1:2], labels=c(1:n), col=color)

################
## Non-metric MDS 
################

library(MASS)
# Euclidean distances between the rows
d <- dist(dat0) 
# Set k=2 for 2 dimensions
mds.fit <- isoMDS(d, k=2) 

# Plot non-metric MDS 
x <- mds.fit$points[,1] 
y <- mds.fit$points[,2]
plot(x, y, main="Nonmetric MDS", type="n")
text(x, y, labels=row.names(dat0), cex=.7, col=color) 

################
## tSNE
################

library(Rtsne)
# Applying tSNE on the data.
tsne <- Rtsne(dat0, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)

# Plot tSNE output
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, col=color)







