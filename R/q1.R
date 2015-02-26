## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-17
###############################################################################


## Install needed packages only need to be done once
install.packages("foreign")
install.packages("survival")
install.packages("KMsurv")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("GGally")

###############################################################################
## Exercise 1b
###############################################################################

## change to git-url
source(paste0("~/src/ki/BiostatIII/newSolutions/R/","ggkm.R"))
## Load the needed packages needs to be done each session
require(foreign) #Needed to read data set from Stata
require(survival)
require(KMsurv)
#require(nlme) #for gsummary
require(dplyr)
require(ggplot2)
require(GGally)

## Get the data for exercise 1
colon_sample <- read.dta("http://biostat3.net/download/colon_sample.dta")

## Create 0/1 outcome variable, create start time variable, 
colon <- mutate(colon_sample,
                death_cancer = ifelse( status == "Dead: cancer", 1, 0),
                start_time = 0)

## Number of events and number lost summarised by month
colon12 <- colon %>%
    mutate(t12m = floor(surv_mm / 12)) %>%
    group_by(t12m) %>%
    summarise(nevent = sum(death_cancer),
              nlost = length(death_cancer)-sum(death_cancer))

##b Life table (NA is added since lifetab wants the end of the last interval)
lifetab(c(colon12$t12m,NA), 35, colon12$nlost, colon12$nevent)

## Kaplan-Meier esimates
mfit <- survfit(Surv(start_time, surv_mm, death_cancer) ~1, data = colon)
summary(mfit)
ggsurv(mfit) + ylab("S(t)") + xlab("Time since diagnosis in months") +
    ggtitle("Kaplan−Meier estimates of cause−specific survival")

## Kaplan-Meier esimates with number at risk
sfit <- survfit(Surv(start_time, surv_mm, death_cancer) ~sex, data = colon)
ggkm(sfit, table=T, pval=F, timeby = 10)





colon_sample <- read.dta("http://biostat3.net/download/colon_sample.dta")


### Construct life table with Actuarial approach ###
attach(colon_sample)
#Recode the status variabel 
death_cancer <- NA 
death_cancer[status =="Dead: cancer"] <- 1
death_cancer[status =="Dead: other" | status == "Alive"] <- 0
#Set 12 month interval
t12m <- floor(surv_mm/12) 
#Combine the two variabels
tall <- data.frame(t12m, death_cancer)
#Calculate the number of death in each interval 
die <-gsummary(tall, sum, groups=t12m)
total <-gsummary(tall, length, groups=t12m)
rm(t12m) #Remove
rm(death_cancer) #Remove
ltab.data <-cbind(die[,1:2], total[,2]) 
detach(colon_sample)
attach(ltab.data)
lt=length(t12m)
t12m[lt+1]=NA
nevent = death_cancer
nlost = total[,2] - death_cancer
mytable <- lifetab(t12m, 35, nlost, nevent)
mytable[,c(1:5, 8)]
detach(ltab.data)

### Construct life table with Kaplan-Meier approach ###
#Construct a variable for starting time
start_time <- rep(0, length(colon_sample$status))
#Recode the status variabel 
death_cancer <- NA 
death_cancer[colon_sample$status =="Dead: cancer"] <- 1
death_cancer[colon_sample$status =="Dead: other" | colon_sample$status == "Alive"] <- 0
#Add the start time and cancer death to the data
colon2 <- data.frame(cbind(colon_sample, death_cancer, start_time))

#Kaplan-Meier estimate
mfit <- survfit(Surv(start_time, surv_mm, death_cancer==1)~1, data = colon2)
summary(mfit)
plot(mfit, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of cause specific survival")

mfit <- survfit(Surv(start_time, surv_mm, death_cancer==1)~sex, data = colon2)
summary(mfit)
plot(mfit, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of cause specific survival")


