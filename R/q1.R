# Date: 2014-01-30
# Purpose: To do the solution for Biostat III exercises in R
# Author: Annika Tillander
###############################################################################

#Install needed packages only need to be done once
install.packages("survival")
install.packages("KMsurv")
install.packages("foreign")
install.packages("lattice")
install.packages("muhaz")
install.packages("nlmez")



###############################################################################
#Exercise 1
###############################################################################

#Load the needed packages needs to be done each session
require(survival)
require(KMsurv)
require(foreign) #Needed to read data set from Stata
require(nlme) #for gsummary

#Get the data for exercise 1
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


