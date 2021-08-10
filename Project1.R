################
## Data Input
################

# Read daily weather data in 2021
dat <- NULL
current.month <- 6

# Update the file name based on the month and save its contents
for (i in 1:current.month){
  i0 <- ifelse(i<10, paste("0", i, sep=""), i)
  mth <- paste("2021", i0, sep="")
  bom <- paste("IDCJDW2801.", mth, ".csv", sep="")
  dat.i <- read.csv(bom, skip=6, check.names=FALSE, na.strings = c("NA", "", " "))
  dat.i[, 1] <- toupper(month.abb[i])
  dat <- rbind(dat, dat.i)
}

dim(dat); head(dat)

################
## Data Cleaning and Preparation
################

# Print out the (sorted) unique values or levels of every variable in the data set
vnames <- colnames(dat)
for (i in 1:length(vnames)){
  print(vnames[i])
  print(table(dat[, i], useNA="ifany"))
  #   print(sort(unique(dat[, i])))
  cat("\n")
}

# Remove “Time of maximum wind gust” from the data set
dat <- dat[, -c(6, 7, 10)]

# Rename the data set
names(dat) <- c("Month", "Date", "MinTemp", "MaxTemp", "Rainfall", "WindGustDir", 
                "WindGustSpeed", "Temp9am", "Humidity9am", "Cloud9am", 
                "WindDir9am", "WindSpeed9am", "Pressure9am", "Temp3pm", 
                "Humidity3pm", "Cloud3pm", "WindDir3pm", "WindSpeed3pm", 
                "Pressure3pm")

# Change the "Calm" values of variables to 0
variablesWithCalm <- c("WindSpeed9am", "WindSpeed3pm")

# Go through all the variables that have Calm as a value and replace it with 0
for (i in 1:length(variablesWithCalm)){
  variable <- variablesWithCalm[i]
  column <- dat[variable]
  column[column == "Calm"] <- 0
  dat[variable] <- as.numeric(unlist(column))
  
}

# Simplify the direction values in order to make the analysis easier
useWindDir <- c("WindGustDir", "WindDir9am", "WindDir3pm")

# Go through all the variables that use direction values and simplify them.
for (i in 1:length(useWindDir)) {
  variable <- useWindDir[i]
  column <- dat[variable]
  column[column == "ENE" | column == "NNE"] <-  "NE"
  column[column == "ESE" | column == "SSE"] <-  "SE"
  column[column == "SSW" | column == "WSW"] <-  "SW"
  column[column == "NNW" | column == "WNW"] <-  "NW"
  dat[variable] <- column
}

# If Rainfall is greater than 1, set RainToday to 1, otherwise to 0
RainToday <- ifelse(dat$Rainfall > 1, 1, 0)
dat$RainToday <- RainToday
dat$RainTomorrow <- c(dat$RainToday[2:nrow(dat)], 0)

# Save an Rdata copy of the data
save(dat, file = "preparedData.Rdata")

################
## Exploratory Data Analysis
################

# Function describe.intreval() does eda for quantitative variables.
describe.interval <- function(dat, cols=1:ncol(dat), by.col=0, wilcox.exact=F, plot.it=TRUE) {
  library("moments")
  # Remove the last row. 
  # Verify that the outcome is binary.
  datEdited <- dat[-nrow(dat),]
  n <- nrow(datEdited)
  vnames <- colnames(datEdited)
  result <- NULL
  if (length(unique(datEdited[, by.col]))!=2) stop("The BY variable is not binary?!!")
  
  for (j in cols) {
    x <- datEdited[,j]
    varname <- vnames[j]
    out <- NULL
    # Obtain missing information for test.
    if (by.col !=0) {
      by.var <- datEdited[, by.col]
      for (k in unique(by.var)){
        x1 <- x[by.var==k]
        n.1 <- sum(is.na(x1), na.rm=T)
        n.2 <- sum(x1=="NA", na.rm=T)
        n.3 <- sum(x1=="", na.rm=T)
        n.miss <- n.1 + n.2 + n.3
        n.complete <- length(x1) -n.miss
        out1 <- c(n.complete, mean(x1, na.rm=T), sd(x1, na.rm=T),
                  kurtosis(x1, na.rm=T), skewness(x1, na.rm=T))
        names(out1) <- paste(c("n.", "mean.", "sd.", "kurto.", "skew."),
                             k, sep="")
        out <- c(out, out1)
      }
      
      # Non-parametic Test Wilcoxon rank sum test
      wilcox <- wilcox.test(x ~ by.var, exact=wilcox.exact, correct=T,
                            na.action=na.omit)
      out <- c(out, c(wilcox$statistic, pvalue.wilcoxon=wilcox$p.value))
      
      # Box Plots
      if (plot.it==T)  {
        par(mfrow=c(1,1))
        boxplot(x ~ as.factor(by.var), col = "orange",
                xlab=vnames[by.col],
                main=paste("Parallel Boxplots of ", varname, sep=""))
      }
    }
    result <- rbind(result, out)
  } 
  result <- as.data.frame(result)
  row.names(result) <- c("MinTemp", "MaxTemp", "Rainfall", 
                         "WindGustSpeed", "Temp9am", "Humidity9am", "Cloud9am", 
                         "WindSpeed9am", "Pressure9am", "Temp3pm", 
                         "Humidity3pm", "Cloud3pm", "WindSpeed3pm", 
                         "Pressure3pm", "RainToday")
  return(result)
}

# Include quantitative columns
# Column 21 contains RainTomorrow
cols.interval <- c(3:5, 7:10, 12:16, 18:20)
summary1 <- describe.interval(dat, cols=cols.interval, by.col=21, wilcox.exact=F)
summary1

# Function describe.cat() does eda for categorical variables.
describe.cat <- function(dat, cols=1:ncol(dat), by.col=0, plot.it=TRUE) {
  # Remove the last row. 
  # Verify that the outcome is binary.
  dat <- dat[-nrow(dat),]
  vnames <- colnames(dat)
  if (by.col !=0) {
    by.var <- dat[, by.col]
    n.groups <- length(unique(na.omit(by.var)))
    if (n.groups !=2) stop("The BY variable is not binary?!!")
  }
  
  result <- NULL
  testResults <- NULL
  for (j in cols) {
    #  Two-way Contingency Tables
    status <- gl(2, 1, labels = c(vnames[by.col], vnames[j]), length = length(dat[, j]))
    tbl0 <- table(dat[, j], status, useNA = "ifany")
    print(vnames[j])
    print(margin.table(tbl0, 1))
    print(prop.table(tbl0, 1))
    print(chisq.test(tbl0))
    
    x <- dat[,j]
    varname <- vnames[j]
    
    # Chi-Squared Test and Fisher's Exact Test
    xlevels <- sort(unique(x))
    chi2.stat <- pvalue.chi2 <- fisher <- NA
    if (length(xlevels) > 1) {
      chi2 <- chisq.test(x, by.var)
      fisher <- fisher.test(x, by.var, hybrid =T)$p.value
      chi2.stat <- chi2$statistic
      cat("\n")
      pvalue.chi2 <- chi2$p.value
    }
    testResults <- rbind(testResults, c(varname=varname, chi2.stat=chi2.stat,
                                        pvalue.chi2=pvalue.chi2, fisher=fisher))
    # Bar Plots
    if (plot.it==T)  {
      par(mfrow=c(1,1))
      rainGroups <- table(dat$RainToday, x)
      barplot(rainGroups, xlab=vnames[by.col],
              main=paste("Barplot of ", varname, sep=""),
              legend=TRUE, beside=TRUE, args.legend=list(title="Rain Tomorrow"))
    }
  }
  colnames(testResults) <- c("Variable", "Chi Square Stat", "Chi Square P-value", "Fisher's Exact")
  return(testResults)
}

# Include categorical columns
# Column 21 contains RainTomorrow
cols.cat <- c(1, 6, 11, 17)
summary2 <- describe.cat(dat, cols=cols.cat, by.col=21)
summary2





