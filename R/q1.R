## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-17
###############################################################################


## Install needed packages only need to be done once
## install.packages("foreign")
## install.packages("survival")
## install.packages("KMsurv")
## install.packages("dplyr")
## install.packages("ggplot2")
## install.packages("GGally")
## install.packages("gridExtra")

###############################################################################
## Exercise 1b
###############################################################################

## @knitr loadDependecies
require(foreign) #Needed to read data set from Stata
require(survival)
require(KMsurv)
require(dplyr)
require(ggplot2)
require(GGally)
require(gridExtra)
setwd(Sys.getenv("PWD")); source("ggkm.R") #For number at risk plot 


## Get the data for exercise 1
colon_sample <- read.dta("http://biostat3.net/download/colon_sample.dta")

## Create 0/1 outcome variable
colon <- mutate(colon_sample,
                death_cancer = ifelse( status == "Dead: cancer", 1, 0))

## @knitr eventsPerMonth
colonByYear <- colon %>%
    mutate(year = floor(surv_yy)) %>%
    group_by(year) %>%
    summarise(nevent = sum(death_cancer),
              nlost = length(death_cancer)-sum(death_cancer))

## @knitr lifeTable
with(colonByYear, lifetab(c(year,tail(year,1)+1), 35, nlost, nevent))[,1:8]

## @knitr KaplanMeier
mfit <- survfit(Surv(surv_mm, death_cancer) ~1, data = colon)
summary(mfit)
ggsurv(mfit) + ylab("S(t)") + xlab("Time since diagnosis in months") +
    ggtitle("Kaplan−Meier estimates of cause−specific survival")

## @knitr KaplanMeierNumberAtRisk
sfit <- survfit(Surv(surv_mm, death_cancer) ~sex, data = colon)
ggkm(sfit, table=T, pval=F, timeby = 10)
