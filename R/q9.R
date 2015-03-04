## Date: 2015-03-03
## Purpose: To do the solution for Biostat III exercise 9 in R
## Author: Annika Tillander,  Johan Zetterqvist
###############################################################################


###############################################################################
## Exercise 9
###############################################################################

## Install needed packages only need to be done once
## install.packages("foreign")
## install.packages("muhaz")
## install.packages("car")

require(survival)
require(dplyr)
require(readstata13)


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

###############################################################################
#Exercise 9
###############################################################################


## Read melanoma data
## and select subcohorts

melanoma.l <-

    tbl_df( read.dta13("http://biostat3.net/download/melanoma.dta") ) %>%

    filter(stage=="Localised") %>%

    mutate(## Create a death indicator
           death_cancer = as.numeric(status=="Dead: cancer"))




melanoma.l2 <-

    mutate(melanoma.l,
           ## Create a death indicator (only count deaths within 120 months)
           death_cancer = death_cancer * as.numeric( surv_mm <= 120),
           ## Create a new time variable
           surv_mm = pmin(surv_mm, 120))


## 9.a

summary( coxfit9a <- coxph(Surv(surv_mm, death_cancer) ~ year8594,
                           data = melanoma.l2) )

## 9.b

summary( coxfit3c <- coxph(Surv(surv_mm, death_cancer) ~ year8594,
                           data = melanoma.l) )

## 9.c

summary( coxfit9c <- coxph(Surv(surv_mm, death_cancer) ~ year8594 + sex + agegrp,
                           data = melanoma.l2) )

## 9.d
## Test if the effect of age is significant
## Use a Wald test with the car package
require(car)
linearHypothesis(coxfit9c,c("agegrp45-59 = 0","agegrp60-74 = 0","agegrp75+ = 0"))

## 9.d

## Test if the effect of age is significant
## Use an LR test with the anova function
## The model without agegrp
summary( coxfit9d <- coxph(Surv(surv_mm, death_cancer) ~ year8594 + sex,
                           data = melanoma.l2) )
## LR test
anova(coxfit9c, coxfit9d)

## 9.e

summary( coxfit9e <- coxph(Surv(surv_mm, death_cancer) ~ year8594 + sex + agegrp,
                           data = melanoma.l2) )

## Split follow up by year
melanoma.spl <- survSplit(melanoma.l2, cut=12*(0:10), end="surv_mm", start="start",
                           event="death_cancer")

## recode start time as a factor
melanoma.spl <- mutate(melanoma.spl,
                       fu = as.factor(start) )

## Run Poisson regression
summary(poisson9e <- glm( death_cancer ~ fu + year8594 + sex + agegrp + offset( log(surv_mm) ),
                         family = poisson,
                         data = melanoma.spl ))

## IRR
IRR(poisson9e)

summary(coxfit9e)

## 9.f

sfit9f <- survfit(Surv(surv_mm, event=death_cancer) ~ 1,
                  data = melanoma.l2 )

## Have a look at the fitted object
str(sfit9f, 1)

head(sfit9f$time)

## Split follow up by year
melanoma.spl2 <- survSplit(melanoma.l2, cut=sfit9f$time, end="surv_mm", start="start",
                           event="death_cancer")

## recode start time as a factor
melanoma.spl2 <- mutate(melanoma.spl2,
                       fu = as.factor(start) )

## Run Poisson regression
poisson9f <- glm( death_cancer ~ fu + year8594 + sex + agegrp + offset( log(surv_mm) ),
                         family = poisson,
                         data = melanoma.spl2 )

## IRR
coef9f <- summary(poisson9f)$coefficients
## We are not interested in nuisance parameters fu1, fu2, etc
npar <- dim(coef9f)[1]
pars <- (npar-4):npar
IRfit9f <- exp( cbind( coef9f[pars, 1:2], coef9f[pars, 1] - 1.96*coef9f[pars, 2], coef9f[pars, 1] +
                        1.96*coef9f[pars, 2] ) )
colnames(IRfit9f) <- c("IRR", "Std. err", "CI_lower", "CI_upper")
IRfit9f

summary(coxfit9e)

## 9.g



