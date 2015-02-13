# Date: 2014-01-30
# Purpose: To do the solution for Biostat III exercises in R
# Author: Johan Zetterqvist
###############################################################################


###############################################################################
#Exercise 6
###############################################################################


#Install needed packages only need to be done once
install.packages("foreign")
install.packages("epiR")

require(survival)
require(epiR)

require(foreign) #Needed to read data set from Stata

#source("http://biostat3.net/download/irEst.R")
source("irEst.R")


# Read diet data
#diet <- data.frame(read.dta("diet.dta"))
diet <- data.frame(read.dta("http://biostat3.net/download/diet.dta"))

# Look at data
head(diet)
summary(diet)

#6a
# Create a variable with unit 1000 person-years
diet$y1k <- diet$y/1000

diet.ir6a <- ir.est("chd","y1k","hieng",diet)
diet.ir6a 

diet.ir6a["high","IR"]/diet.ir6a["low","IR"]

#6b
poisson6b <- glm( chd ~ hieng + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6b )
cbind(coef(poisson6b),confint(poisson6b))
exp(cbind(coef(poisson6b),confint(poisson6b)))


#6c
hist6c <- hist(diet$energy, breaks=25, probability=TRUE )
curve(dnorm(x, mean=mean(diet$energy), sd=sd(diet$energy)), add=TRUE)
quantile(diet$energy, probs=c(0.01,0.05,0.1,0.25,0.5,0.75,0.90,0.95,0.99))
summary(diet$energy)
# For kurtosis and skewness, see package e1071

#6d
# Create a factor
diet$eng3 <- cut(diet$energy, breaks=c(1500,2500,3000,4500),labels=c("low","medium","high"))
summary(diet$eng3)

#6e
diet.ir6e <- ir.est("chd","y1k","eng3",diet)
diet.ir6e

diet.ir6e["medium","IR"]/diet.ir6e["low","IR"]
diet.ir6e["high","IR"]/diet.ir6e["low","IR"]

#6f
dummies6f <- model.matrix(~eng3,data=diet)
diet$X2 <- dummies6f[,"eng3medium"]
diet$X3 <- dummies6f[,"eng3high"]
diet$X1 <- (1-diet$X2)*(1-diet$X3)
colSums(diet[c("X1","X2","X3")])

#6g
head(diet[which(diet$eng3=="low"),c("energy","eng3","X1","X2","X3")])
head(diet[which(diet$eng3=="medium"),c("energy","eng3","X1","X2","X3")])
head(diet[which(diet$eng3=="high"),c("energy","eng3","X1","X2","X3")])

#6h
poisson6h <- glm( chd ~ X2 + X3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6h )

cbind(coef(poisson6h),confint(poisson6h))
exp(cbind(coef(poisson6h),confint(poisson6h)))

#6i
poisson6i <- glm( chd ~ X1 + X3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6i )

cbind(coef(poisson6i),confint(poisson6i))
exp(cbind(coef(poisson6i),confint(poisson6i)))

#6j
poisson6j <- glm( chd ~ eng3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6j )

cbind(coef(poisson6j),confint(poisson6j))
exp(cbind(coef(poisson6j),confint(poisson6j)))

#6k
sums6k <- colSums(diet[c("chd","y")])
sums6k[1]/sums6k[2]
