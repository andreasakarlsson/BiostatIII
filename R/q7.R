# Date: 2014-01-30
# Purpose: To do the solution for Biostat III exercises in R
# Author: Johan Zetterqvist
###############################################################################

###############################################################################
#Exercise 7
###############################################################################

#Install needed packages only need to be done once
install.packages("survival")
install.packages("foreign")
install.packages("epiR")
install.packages("muhaz")
install.packages("Epi")
install.packages("car")

require(survival)
require(epiR)
require(muhaz) #for hazard estimates
require(foreign) #Needed to read data set from Stata
require(Epi)
require(car)

source("http://biostat3.net/download/R/irEst.R")
#source("irEst.R")

# Read melanoma data
melanoma <- subset.data.frame(read.dta("http://biostat3.net/download/melanoma.dta"), stage=="Localised")

# Look at data
head(melanoma)
summary(melanoma)

# Lexis requires time variables to be numeric and have the same units
# we will use years
melanoma$dx<-1970+as.numeric(melanoma$dx)/365.24
melanoma$exit1<-1970+as.numeric(melanoma$exit)/365.24
melanoma$bdate<-1970+as.numeric(melanoma$bdate)/365.24

#Create a new event indicator
melanoma$dead_cancer1 <- as.numeric(melanoma$status=="Dead: cancer")

# Create new variables where
# patients who survived more than 10 years are censored at 10 years
melanoma$dead_cancer2 <- ifelse(melanoma$surv_mm>120,0,melanoma$dead_cancer1)
melanoma$exit2 <- ifelse(melanoma$surv_mm>120,melanoma$dx + 10,melanoma$exit1)
melanoma$surv_yy2 <- ifelse(melanoma$surv_mm>120,10,melanoma$surv_yy)
melanoma$surv_mm2 <- ifelse(melanoma$surv_mm>120,120,melanoma$surv_mm)

# Check data again
head(melanoma)

#Create a new variable base on the number of months
# translated to years (more accurate than surv_yy)
melanoma$surv_yy1 <- melanoma$surv_mm/12

head( melanoma)

#7a

# Estimate Kaplan-Meier curve
#sfit7a1 <- survfit(Surv(surv_yy, event=dead_cancer1)~year8594, data = melanoma)
#plot(sfit7a1,mark.time=F,xlab="Years since diagnose",ylab="S(t)",col=c("blue","red"),lty=c("solid","dashed"))
#legend("bottomleft",legend=levels(melanoma$year8594),col=c("blue","red"),lty=c("solid","dashed"))

sfit7a1 <- survfit(Surv(surv_mm, event=dead_cancer1)~year8594, data = melanoma)
plot(sfit7a1,mark.time=F,xscale=12,xlab="Years since diagnose",ylab="S(t)",col=c("blue","red"),lty=c("solid","dashed"))
legend("bottomleft",legend=levels(melanoma$year8594),col=c("blue","red"),lty=c("solid","dashed"))

# Calculate smoothed hazard functions per diagnose period
# The plot.muhaz cannot rescale so we create a new variable
# based on the number of months
# translated to years (more accurate than surv_yy)
melanoma$surv_yy1 <- melanoma$surv_mm/12

hazfit7a_early <- muhaz(times=melanoma$surv_yy1,delta=melanoma$dead_cancer1,subset=melanoma$year8594=="Diagnosed 75-84")
hazfit7a_late <- muhaz(times=melanoma$surv_yy1,delta=melanoma$dead_cancer1,subset=melanoma$year8594=="Diagnosed 85-94",max.time=10.5)

# Create a new window so we can compare survival curves with hazards
# windows() # On Windows
# X11() # On linux
quartz() # On Mac

# Plot smoothed hazards
plot(hazfit7a_early,xlab="Years since diagnose",col="blue",lty="solid")
lines(hazfit7a_late$est.grid,hazfit7a_late$haz.est,col="red",lty="dashed")
legend("topright",legend=levels(melanoma$year8594),col=c("blue","red"),lty=c("solid","dashed"))


#7b
# Create a variable with unit 1000 person-years
melanoma$y1k1 <- melanoma$surv_yy1/1000

ir7b <- ir.est("dead_cancer1","y1k1","year8594",melanoma)
ir7b

#7c

# Create a new death indicator (only count deaths within 120 months)
melanoma$dead_cancer2 <- melanoma$dead_cancer1*(melanoma$surv_mm<=120)
# Create a new time variable
melanoma$surv_mm2 <- ifelse(melanoma$surv_mm<=120,melanoma$surv_mm,120)
melanoma$surv_yy2 <- melanoma$surv_mm2/12
# Create a variable with unit 1000 person-years
melanoma$y1k2 <- melanoma$surv_yy2/1000

ir7c <- ir.est("dead_cancer2","y1k2","year8594",melanoma)
ir7c

ir7c["Diagnosed 85-94","IR"]/ir7c["Diagnosed 75-84","IR"]

poisson7c <- glm( dead_cancer2 ~ year8594 + offset( log( y1k2 ) ), family=poisson, data=melanoma )
summary( poisson7c )
# IR
exp(coef(poisson7c))
# with confidence intervals (takes some time to compute)
#exp(cbind(coef(poisson7c),confint(poisson7c)))

#7d

# Define a Lexis object

melanoma.lex <- Lexis( entry = list( per=dx, fu=0), 
			exit = list( per=exit2 ),
 			exit.status = dead_cancer2,
 			data = melanoma )

# Check the new variables
head(melanoma.lex)

# Split follow-up time into annual intervals
melanoma.spl <- splitLexis(melanoma.lex, breaks=seq(0,10,1), time.scale="fu")

# Look at data
# lex.Xst is the event indicator and
# fu is the year since diagnose
# lex.dur is the observation time in each interval
head( melanoma.spl, 20 )
 
#7e
# Recalculate the time variable in units 1000 years 
melanoma.spl$y1k2 <- melanoma.spl$lex.dur/1000
# Calculate incidence rates
ir7e <- ir.est("lex.Xst","y1k2","fu",melanoma.spl)
ir7e

# Plot by year
matplot(0:9,ir7e[,c("IR","CI lower","CI upper")],lty=c("solid","dashed","dashed"),col=c("black","gray","gray"),type="l",main="Cancer deaths per 1000 person-year by years since diagnose",ylab="IR",xlab="Years since diagnose")

#7f

# Calculate smoothed hazard functions 
hazfit7e <- muhaz(times=melanoma$surv_yy2,delta=melanoma$dead_cancer2,max.time=9.5)

# Plot smoothed hazards
quartz()
plot(hazfit7e,xlab="Years since diagnose")

#7g

# recode fu as a factor
melanoma.spl$fu <- as.factor(melanoma.spl$fu)
# Run Poisson regression
summary(poisson7g <- glm( melanoma.spl$lex.Xst ~ fu + offset( log(lex.dur) ), family=poisson, data=melanoma.spl ))
# IR
exp(coef(poisson7g))
# with confidence intervals (takes some time to compute)
#exp(cbind(coef(poisson7g),confint(poisson7g)))

#7h
summary(poisson7h <- glm( lex.Xst ~ fu + year8594 + offset( log(lex.dur) ), family=poisson, data=melanoma.spl ))
# IR
exp(coef(poisson7h))
# with confidence intervals (takes some time to compute)
#exp(cbind(coef(poisson7h),confint(poisson7h)))

# Add interaction term
summary(poisson7h2 <- glm( lex.Xst ~ fu*year8594 + offset( log(lex.dur) ), family=poisson, data=melanoma.spl ))
# IR
exp(coef(poisson7h2))


#7i
summary(poisson7i <- glm( lex.Xst ~ fu + agegrp + year8594 + sex + offset( log(lex.dur) ), family=poisson, data=melanoma.spl ))
# IR
exp(coef(poisson7i))
# with confidence intervals (takes some time to compute)
#exp(cbind(coef(poisson7i),confint(poisson7i)))

# Test if the effect of age is significant
linearHypothesis(poisson7i,c("agegrp45-59 = 0","agegrp60-74 = 0","agegrp75+ = 0"))
# ADVANCED:
# Alternative by comparing deviances
# poisson7i_2 <- update(poisson7i,. ~ . - agegrp)
# anova(poisson7i_2,poisson7i,test="Chisq")

#7j
summary(poisson7j <- glm( lex.Xst ~ fu + agegrp + year8594*sex + offset( log(lex.dur) ), family=poisson, data=melanoma.spl ))

#7k
# 'hand' calculations
hz7k <- exp(coef(poisson7j))
hz7k["sexFemale"]
hz7k["sexFemale"]*hz7k["year8594Diagnosed 85-94:sexFemale"]

# ADVANCED:
# 'hand' calculations with confidence intervals
# use estimates with covariance matrix from glm 
# to obtain estimates with ci's
length(coef(poisson7j))
linvec <- c(rep(0,14),c(1,1))
coef.Female8594 <- crossprod(coef(poisson7j),linvec) 
var.Female8594 <- t(linvec)%*%vcov(poisson7j)%*%linvec
ci.d <- 1.96*sqrt(var.Female8594)
ci.lower <- coef.Female8594 - ci.d
ci.upper <- coef.Female8594 + ci.d
exp(c(coef.Female8594,ci.lower,ci.upper))
exp(c(coef(poisson7j)["sexFemale"],confint(poisson7j)["sexFemale",]))

# Create dummies and Poisson regression
#levels(melanoma.spl$sex)
#levels(melanoma.spl$year8594)
melanoma.spl$femaleEarly <- melanoma.spl$sex=="Female" & melanoma.spl$year8594=="Diagnosed 75-84"
melanoma.spl$femaleLate <- melanoma.spl$sex=="Female" & melanoma.spl$year8594=="Diagnosed 85-94"
summary(poisson7k <- glm( lex.Xst ~ fu + agegrp + year8594 + femaleEarly + femaleLate + offset( log(lex.dur) ), family=poisson, data=melanoma.spl ))
# IR
exp(coef(poisson7k))
# with confidence intervals (takes some time to compute)
#exp(cbind(coef(poisson7k),confint(poisson7k)))

#7l

poisson7l.list <- by(melanoma.spl, melanoma.spl$year8594, 
					function(x) {glm( lex.Xst ~ fu + agegrp + sex + offset( log(lex.dur) ), 
									family=poisson, data=x )
								}
					)

poisson7l.early <- poisson7l.list[["Diagnosed 75-84"]]
# IR
exp(coef(poisson7l.early))
# with confidence intervals (takes some time to compute)
#exp(cbind(coef(poisson7h),confint(poisson7h)))

poisson7l.late <- poisson7l.list[["Diagnosed 85-94"]]
# IR
exp(coef(poisson7l.late))
# with confidence intervals (takes some time to compute)
#exp(cbind(coef(poisson7l.late),confint(poisson7l.late)))

# compare with results in i
exp(coef(poisson7i))

# compare with results in j
exp(coef(poisson7j))

summary(poisson7l.early)

# Poisson-regression with effects specific for diagnose period
summary(glm( lex.Xst ~ fu + fu:year8594 + agegrp + agegrp:year8594 + sex*year8594 + offset( log(lex.dur) ), 
									family=poisson, data=melanoma.spl ))
