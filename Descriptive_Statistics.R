# NCBS Biostatistics

# Exercise 1: descriptive statistics

# Set you working directory (path where you save scripts and data)
setwd ("/Users/kartikmajila/Downloads/Courses/Biostatistics") # replace with your desired path

# read the input data file
df<-read.csv("normal_skewed_distr.csv")
?head # to get help page
head(df)

# plot histogram
?hist # to get help page
hist(df$geneA)

# add x-axis lab and title to the plot
hist(df$geneA, xlab="expression of geneA", main="Histogram of geneA")

# measure of central tendency
# compute mean
?mean()
mean(df$geneA)
# add mean value to the plot
?abline
abline(v=mean(df$geneA), col="red", lw=3) 

# compute median
?median
median(df$geneA)
# add median value to the plot
abline(v=median(df$geneA), col="green", lw=3)

# compute mode
# no default function available for mode
# so, we first compute the frequeny table using table function in R
y <- table(df$geneA)
head(y)
# then find the value which have maximum frequency
print(names(y)[which(y==max(y))])
# add mode to the plot
abline(v=names(y)[which(y==max(y))], col="blue", lw=3)

# compute geometric mean
# log transform the values
?log
log_t <- log(df$geneA[df$geneA>0])
# compute mean of log transformed values
m_log_t <- mean(log_t)
# anti-log
exp(m_log_t)
# add geometric mean to the plot
abline(v=exp(m_log_t), col="orange", lw=3)

# measure of spread
# range
max(df$geneA)-min(df$geneA)

# quantiles
?quantile
quantile(df$geneA)
quantile(df$geneA, probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(df$geneA, probs = 0.5)

# summarize the data
summary(df$geneA)

boxplot(df$geneA)
quantile(df$geneA)

# map the gene exp value corresponding to each quartile value
abline(h=7, lty=2) # 0% - minimum (Q0)
abline(h=40, lty=2) # 25% - lower quartile (Q1)
abline(h=50, lty=2) # 50% - median (Q2)
abline(h=61, lty=2) # 75% - upper quartile (Q3)
abline(h=89, lty=2) # 100% - maximum (Q4)

# Interquartile range = upper quartile - lower quartile
# 61-40 # values taken from above or use IQR function
IQR(df$geneA)

# whiskers
# upper whisker = Q3 + 1.5 IQR
upper_whisker = quantile(df$geneA, probs=0.75) + 1.5 * IQR(df$geneA)
upper_whisker
abline(h=upper_whisker, lty=2, col="red")

# find the value estimated by R
?boxplot.stats
boxplot.stats(df$geneA)

# find the largest value in the table that is lesser than or equal to Q3+1.5*IQR
head(sort(df$geneA[df$geneA<=upper_whisker], decreasing=TRUE))
abline(h=89, lty=1, col="red")

# lower whisker  = Q1 - 1.5 IQR
lower_whisker = quantile(df$geneA, probs=0.25) - 1.5 * IQR(df$geneA)
lower_whisker
abline(h=lower_whisker, lty=2, col="red")

# find the value estimated by R
boxplot.stats(df$geneA)

# find the lowest value in the table that is greater than or equal to Q1-1.5*IQR
head(sort(df$geneA[df$geneA>=lower_whisker]))
abline(h=10, lty=1, col="red")

# standard deviation
sd(df$geneA)

# variance
var(df$geneA)

# coefficient of the variance
sd(df$geneA) / mean(df$geneA)* 100
#####################################
# TODO:
# Run the above exercise for geneB and geneC (starting from plot histogram, then measure of central tendency and spread)
#For geneB
head(df)

# plot histogram
hist(df$geneB)

# add x-axis lab and title to the plot
hist(df$geneB, xlab="expression of geneB", main="Histogram of geneB")

# measure of central tendency
# compute mean
print("The mean expression value for gene B = ") 
print(mean(df$geneB))
# add mean value to the plot
abline(v=mean(df$geneB), col="red", lw=3) 

# compute median
print("The median expression value for gene B = ") 
print(median(df$geneB))
# add median value to the plot
abline(v=median(df$geneB), col="green", lw=3)

# compute mode
# no default function available for mode
# so, we first compute the frequeny table using table function in R
y <- table(df$geneB)
head(y)
# then find the value which have maximum frequency
print("Mode of geneB expression = ")
print(names(y)[which(y==max(y))])
# add mode to the plot
abline(v=names(y)[which(y==max(y))], col="blue", lw=3)

# compute geometric mean
# log transform the values
log_t <- log(df$geneB[df$geneB>0])
# compute mean of log transformed values
m_log_t <- mean(log_t)
# anti-log
exp(m_log_t)
# add geometric mean to the plot
abline(v=exp(m_log_t), col="orange", lw=3)

# measure of spread
# range
max(df$geneB)-min(df$geneB)

# quantiles
quantile(df$geneB)
quantile(df$geneB, probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(df$geneB, probs = 0.5)

# summarize the data
summary(df$geneB)

boxplot(df$geneB)
quantile(df$geneB)

# map the gene exp value corresponding to each quartile value
abline(h=0, lty=2) # 0% - minimum (Q0)
abline(h=4, lty=2) # 25% - lower quartile (Q1)
abline(h=13, lty=2) # 50% - median (Q2)
abline(h=24, lty=2) # 75% - upper quartile (Q3)
abline(h=72, lty=2) # 100% - maximum (Q4)

# Interquartile range = upper quartile - lower quartile
# 61-40 # values taken from above or use IQR function
IQR(df$geneB)

# whiskers
# upper whisker = Q3 + 1.5 IQR
upper_whisker = quantile(df$geneB, probs=0.75) + 1.5 * IQR(df$geneB)
upper_whisker
abline(h=upper_whisker, lty=2, col="red")

# find the value estimated by R
boxplot.stats(df$geneB)

# find the largest value in the table that is lesser than or equal to Q3+1.5*IQR
head(sort(df$geneB[df$geneB<=upper_whisker], decreasing=TRUE))
abline(h=72, lty=1, col="red")

# lower whisker  = Q1 - 1.5 IQR
lower_whisker = quantile(df$geneB, probs=0.25) - 1.5 * IQR(df$geneB)
lower_whisker
abline(h=lower_whisker, lty=2, col="red")

# find the value estimated by R
boxplot.stats(df$geneB)

# find the lowest value in the table that is greater than or equal to Q1-1.5*IQR
head(sort(df$geneB[df$geneB>=lower_whisker]))
abline(h=0, lty=1, col="red")

# standard deviation
print("sd(geneB")
print(sd(df$geneB))

# variance
print("var(geneB)")
print(var(df$geneB))

# coefficient of the variance
print("Coefficient of Variance for geneB = ")
print(sd(df$geneB) / mean(df$geneB)* 100)

#####################################
#For geneC
head(df)

# plot histogram
hist(df$geneC)

# add x-axis lab and title to the plot
hist(df$geneC, xlab="expression of geneB", main="Histogram of geneB")

# measure of central tendency
# compute mean
print("The mean expression value for gene C = ") 
print(mean(df$geneC))
# add mean value to the plot
abline(v=mean(df$geneC), col="red", lw=3) 

# compute median
print("The median expression value for gene C = ")
print(median(df$geneC)))
# add median value to the plot
abline(v=median(df$geneC), col="green", lw=3)

# compute mode
# no default function available for mode
# so, we first compute the frequeny table using table function in R
y <- table(df$geneC)
head(y)
# then find the value which have maximum frequency
print("Mode for geneC expression = ")
print(names(y)[which(y==max(y))])
# add mode to the plot
abline(v=names(y)[which(y==max(y))], col="blue", lw=3)

# compute geometric mean
# log transform the values
log_t <- log(df$geneC[df$geneC>0])
# compute mean of log transformed values
m_log_t <- mean(log_t)
# anti-log
exp(m_log_t)
# add geometric mean to the plot
abline(v=exp(m_log_t), col="orange", lw=3)

# measure of spread
# range
max(df$geneC)-min(df$geneC)

# quantiles
quantile(df$geneC)
quantile(df$geneC, probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(df$geneC, probs = 0.5)

# summarize the data
summary(df$geneC)

boxplot(df$geneC)
quantile(df$geneC)

# map the gene exp value corresponding to each quartile value
abline(h=15, lty=2) # 0% - minimum (Q0)
abline(h=60, lty=2) # 25% - lower quartile (Q1)
abline(h=73, lty=2) # 50% - median (Q2)
abline(h=84, lty=2) # 75% - upper quartile (Q3)
abline(h=99, lty=2) # 100% - maximum (Q4)

# Interquartile range = upper quartile - lower quartile
# 61-40 # values taken from above or use IQR function
IQR(df$geneC)

# whiskers
# upper whisker = Q3 + 1.5 IQR
upper_whisker = quantile(df$geneC, probs=0.75) + 1.5 * IQR(df$geneC)
upper_whisker
abline(h=upper_whisker, lty=2, col="red")

# find the value estimated by R
boxplot.stats(df$geneC)

# find the largest value in the table that is lesser than or equal to Q3+1.5*IQR
head(sort(df$geneC[df$geneC<=upper_whisker], decreasing=TRUE))
abline(h=99, lty=1, col="red")

# lower whisker  = Q1 - 1.5 IQR
lower_whisker = quantile(df$geneC, probs=0.25) - 1.5 * IQR(df$geneC)
lower_whisker
abline(h=lower_whisker, lty=2, col="red")

# find the value estimated by R
boxplot.stats(df$geneC)

# find the lowest value in the table that is greater than or equal to Q1-1.5*IQR
head(sort(df$geneC[df$geneC>=lower_whisker]))
abline(h=25, lty=1, col="red")

# standard deviation
print("sd(geneC)")
print(sd(df$geneC))

# variance
print("var(geneC)")
print(var(df$geneC))

# coefficient of the variance
print("Coefficient of Variance for geneC = ")
print(sd(df$geneC) / mean(df$geneC)* 100)



