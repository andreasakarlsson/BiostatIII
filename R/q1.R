## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-17, 2016-03-01
###############################################################################



###############################################################################
## Exercise 1b
###############################################################################

## @knitr loadDependecies
require(foreign)  # needed to read data set from Stata
require(dplyr)    # for data manipulation
require(KMsurv)   # for life-tables
require(survival) # for Kaplan-Meier

## Get the data for exercise 1
colon_sample <- read.dta("http://biostat3.net/download/colon_sample.dta")

## Create 0/1 outcome variable
colon <-colon_sample %>%
    mutate(death_cancer = ifelse( status == "Dead: cancer", 1, 0))

## @knitr eventsPerYear
colonByYear <- colon %>%
    mutate(year = floor(surv_yy)) %>%     # floor to whole years
    group_by(year) %>%                    # summarise within each year
    summarise(nevent = sum(death_cancer), # number of events
              nlost = n() - nevent)       # number lost to follow-up


## @knitr lifeTable
with(colonByYear,                           # using the colonByYear data
     lifetab(tis = c(year, tail(year,1)+1), # should be one element longer for the intervals
             ninit = nrow(colon),           # number of individuals at the start
             nlost = nlost,                 # number lost for each interval
             nevent = nevent)) %>%          # number of events for each interval
    round(2)

## @knitr KaplanMeier
mfit <- survfit(Surv(surv_mm, death_cancer) ~ 1, data = colon) # make Kaplan-Meier estimates
summary(mfit)                                                  # print Kaplan-Meier table
plot(mfit,                                                     # plot Kaplan-Meier curve
     ylab="S(t)",
     xlab="Time since diagnosis in months",
     main = "Kaplan−Meier estimates of cause−specific survival")
