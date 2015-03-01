## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist, 2014-01-30
## Edited: Andreas Karlsson, 2015-03-01
###############################################################################

## ## Install needed packages only need to be done once
## install.packages("foreign")  #Needed to read data set from Stata
## install.packages("epiR")
## install.packages("survival")
## install.packages("ggplot2")
## install.packages("RCurl")

###############################################################################
## Exercise 6
###############################################################################
## @knitr loadDependecies
require(survival)
require(epiR)
require(foreign) #Needed to read data set from Stata
require(ggplot2)
require(RCurl)
eval(expr=parse(text=getURL("https://raw.githubusercontent.com/andreasakarlsson/BiostatIII/master/R/helpFuncBioIII.R")))

## @knitr loadPreprocess
diet <- data.frame(read.dta("http://biostat3.net/download/diet.dta"))

## Look at data
head(diet)
summary(diet)


## @knitr 6a_ir
diet$y1k <- diet$y/1000
diet.ir6a <- ir.est(diet, grp="hieng", event="chd", p.time="y1k")
print(diet.ir6a)

## IRR
diet.ir6a[diet.ir6a$hieng=="high","IR"] / diet.ir6a[diet.ir6a$hieng=="low","IR"]


## @knitr 6b_ir
poisson6b <- glm( chd ~ hieng + offset( log( y1k ) ), family=poisson, data=diet)
summary(poisson6b)
exp(cbind(coef(poisson6b),confint(poisson6b)))


## @knitr 6c_energyDist
ggplot(diet, aes(energy)) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    stat_function(fun = dnorm, args = list(mean = mean(diet$energy), sd = sd(diet$energy)),
                  geom='area', alpha=.2, fill="#FF6666",colour="black")
quantile(diet$energy, probs=c(0.01,0.05,0.1,0.25,0.5,0.75,0.90,0.95,0.99))
# For kurtosis and skewness, see package e1071

## @knitr 6d_engCat
diet$eng3 <- cut(diet$energy, breaks=c(1500,2500,3000,4500),labels=c("low","medium","high"))
diet %>% group_by(eng3) %>%
    summarise(Freq = n(), Percent = n()/dim(diet)[1]) %>%
    mutate(Cum = cumsum(Percent))

## @knitr 6e_irEng
diet.ir6e <- ir.est(diet, grp="eng3", event="chd", p.time="y1k")
print(diet.ir6e)
diet.ir6e[diet.ir6e$eng3=="medium","IR"]/diet.ir6e[diet.ir6e$eng3=="low","IR"]
diet.ir6e[diet.ir6e$eng3=="high","IR"]/diet.ir6e[diet.ir6e$eng3=="low","IR"]

## @knitr 6f_irEng
dummies6f <- model.matrix(~eng3,data=diet)
diet$X2 <- dummies6f[,"eng3medium"]
diet$X3 <- dummies6f[,"eng3high"]
diet$X1 <- (1-diet$X2)*(1-diet$X3)
colSums(diet[c("X1","X2","X3")])

## @knitr 6g_irEng
head(diet[which(diet$eng3=="low"),c("energy","eng3","X1","X2","X3")])
head(diet[which(diet$eng3=="medium"),c("energy","eng3","X1","X2","X3")])
head(diet[which(diet$eng3=="high"),c("energy","eng3","X1","X2","X3")])

## @knitr 6h
poisson6h <- glm( chd ~ X2 + X3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6h)
exp(cbind(coef(poisson6h),confint(poisson6h)))

## @knitr 6i
poisson6i <- glm( chd ~ X1 + X3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6i )
exp(cbind(coef(poisson6i),confint(poisson6i)))

## @knitr 6j
poisson6j <- glm( chd ~ eng3 + offset( log( y1k ) ), family=poisson, data=diet )
summary( poisson6j )
exp(cbind(coef(poisson6j),confint(poisson6j)))

## @knitr 6k
sums6k <- colSums(diet[c("chd","y")])
sums6k[1]/sums6k[2]

