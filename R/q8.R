## Purpose: To do the solution for Biostat III exercises in R
## Author: Johan Zetterqvist, 2014-01-30
## Edited: Andreas Karlsson, 2015-03-01
##         Benedicte Delcoigne, 2015-03-03
###############################################################################

## ## Install needed packages only need to be done once
## install.packages("foreign")  #Needed to read data set from Stata
## install.packages("epiR")
## install.packages("survival")
## install.packages("ggplot2")
## install.packages("RCurl")
## install.packages("dplyr")

###############################################################################
## Exercise 8
###############################################################################
## @knitr loadDependecies
require(survival)
require(foreign) #Needed to read data set from Stata
require(ggplot2)
require(dplyr)
require(muhaz) #for hazard estimates


###########################################
### A help function to calculate ###
### and print incidence (hazard) ratios
### from a fitted poisson regression
### from glm
###########################################

IRR <- function(fit){
    summfit <- summary(fit )$coefficients
    IRfit <- exp( cbind( summfit[, 1:2], summfit[, 1] - 1.96*summfit[, 2], summfit[, 1] +
                        1.96*summfit[, 2] ) )
    colnames(IRfit) <- c("IRR", "Std. err", "CI_lower", "CI_upper")
    print(IRfit)
}


## @knitr loadPreprocess
diet <- read.dta("http://biostat3.net/download/diet.dta")
head(diet)
summary(diet)

## @knitr 8a1_first_step_att_age
## Creating a variable for attained age
diet <- diet %>% mutate(att_age = as.numeric(dox - dob) / 365.24)
summary(diet$att_age)


## @knitr 8a1_Haz_att_age
## Run muhaz smoother for Kaplan-Meier hazards ones per strata, note that this is left censured data and the muhaz smoother is not made for this.
dietHiengAge <- diet %>%  group_by(hieng) %>%
    do( h = muhaz(times = .$att_age, delta = .$chd, min.time = min(.$att_age[as.logical(.$chd)]), max.time = max(.$att_age[as.logical(.$chd)])-2)) %>%
    do( data.frame(Hazard = .$h$haz.est, Time = .$h$est.grid, Strata = .$hieng))
ggplot(dietHiengAge, aes(x=Time, y=Hazard, colour= Strata)) +
    geom_line() + ggtitle("Smoothed Hazard estimates") + theme(legend.position="bottom")

## Have a look at the raw hazard if you dare...
## hazard <- data.frame(with( diet, kphaz.fit(time=att_age, status=chd, strata=hieng)))
## ggplot(hazard, aes(time, haz, group=strata, color=factor(strata))) + geom_point()

## @knitr 8a2_Haz_time_entry
diet <- diet %>% mutate(t_entry = as.numeric(dox - doe) / 365.24)
summary(diet$t_entry)

dietHiengEntry <- diet %>%  group_by(hieng) %>%
    do( h = muhaz(times = .$t_entry, delta = .$chd,
            min.time = min(.$t_entry[as.logical(.$chd)]), max.time = max(.$t_entry[as.logical(.$chd)]))) %>%
    do( data.frame(Hazard = .$h$haz.est, Time = .$h$est.grid, Strata = .$hieng))
ggplot(dietHiengEntry, aes(x=Time, y=Hazard, colour= Strata)) +
    geom_line() + ggtitle("Smoothed Hazard estimates") + theme(legend.position="bottom")


## @knitr 8b_ir
diet <- mutate(diet, y1k = y / 1000)
poisson8b <- glm( chd ~ hieng + offset( log( y1k ) ), family=poisson, data=diet)
summary(poisson8b)
exp(cbind(coef(poisson8b),confint(poisson8b)))


## @knitr 8c_ir
## Create BMI variable
diet$bmi <- diet$weight/((diet$height/100)^2)

## Create orderly varable instead of categorical, start at zero
diet <- diet %>% mutate(jobNumber = match(job, c("driver","conductor","bank")) -1)

poisson8c <- glm( chd ~ hieng + jobNumber + bmi + offset( log( y1k ) ), family=poisson, data=diet)
summary(poisson8c)
exp(cbind(coef(poisson8c),confint(poisson8c)))


## @knitr 8d
#timeBand <- survSplit(diet, cut=c(0,30,60,72), end="doxNum", start="doeNum", event="chd")


timeBand <- survSplit(diet, cut=c(0,30,60,72), end="att_age", start="start", event="chd")

## recode start time as a factor
timeBand <- mutate(diet, fu = as.factor(start))

## Output the first three individuals
timeBand %>% select(id, ) %>% filter(id<=3) %>% arrange(id, att_age)


## This was correct ysester day!!!
## @knitr 8e

diet.spl.t_entry <- survSplit(diet, cut=c(0, 5, 10, 15, 22), end="t_entry", start="start", event="chd")
diet.spl.t_entry %>% filter(id<=3) %>% arrange(id, t_entry)

diet.spl.t_entry <- mutate(diet.spl.t_entry,
                           fu = as.factor(start) ,
                           risk_time = (t_entry-start))

poisson8e1 <- glm( chd ~ fu + hieng + offset( log( risk_time) ),
                 family=poisson,
                 data=diet.spl.t_entry )

summary( poisson8e1 )

## IRR

IRR(poisson8e1)

poisson8e2 <- glm( chd ~ fu + hieng + jobNumber + bmi + offset( log( t_entry) ),
                 family=poisson,
                 data=diet.spl.t_entry )

summary( poisson8e2 )

## IRR
IRR(poisson8e2)



tmp <- survSplit(diet, cut=c(0,30,60,72), end="att_age", start="start",
                           event="chd")
