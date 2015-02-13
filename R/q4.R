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
#Exercise 4
###############################################################################

require(foreign) #Needed to read data set from Stata
require(survival) #for Surv and survfit
require(muhaz) #for hazard estimates
require(lattice) #for grouped plots
require(KMsurv)
require(nlme) #for gsummary
#The data
melanoma <- read.dta("http://biostat3.net/download/melanoma.dta")
melanoma_sub <- melanoma[melanoma$stage=="Localised",]
#a
#Construct a variable for starting time
start_time <- rep(0, length(melanoma_sub$status))
#Recode the status variabel 
death_cancer <- NA 
death_cancer[melanoma_sub$status =="Dead: cancer"] <- 1
death_cancer[melanoma_sub$status =="Dead: other" | melanoma_sub$status == "Alive" | melanoma_sub$status == "Lost to follow-up"] <- 0
death_all <- NA 
death_all[melanoma_sub$status =="Dead: cancer"|melanoma_sub$status =="Dead: other"] <- 1
death_all[melanoma_sub$status == "Alive" | melanoma_sub$status == "Lost to follow-up"] <- 0
#Add the start time and cancer death to the data
melanoma3 <- data.frame(cbind(melanoma_sub, death_cancer, death_all, start_time))
#a
#Kaplan-Meier estimate
mfit_years <- survfit(Surv(start_time, surv_yy, death_cancer==1)~1, data = melanoma3)
summary(mfit_years)

mfit_months <- survfit(Surv(start_time, surv_mm, death_cancer==1)~1, data = melanoma3)
summary(mfit_months)

### Construct life table with Actuarial approach ###
#For years
attach(melanoma3)
#Combine the two variabels
tall <- data.frame(surv_yy, death_cancer)
#Calculate the number of death in each interval 
die <-gsummary(tall, sum, groups=surv_yy)
total <-gsummary(tall, length, groups=surv_yy)
rm(surv_yy) #Remove
rm(death_cancer) #Remove
ltab.data <-cbind(die[,1:2], total[,2]) 
detach(melanoma3)
attach(ltab.data)
lt=length(surv_yy)
surv_yy[lt+1]=NA
nevent = death_cancer
nlost = total[,2] - death_cancer
mytable <- lifetab(surv_yy, 5318, nlost, nevent)
mytable[,c(1:5, 8)]
detach(ltab.data)

#For months
attach(melanoma3)
#Combine the two variabels
tall <- data.frame(surv_mm, death_cancer)
#Calculate the number of death in each interval 
die <-gsummary(tall, sum, groups=surv_mm)
total <-gsummary(tall, length, groups=surv_mm)
rm(surv_mm) #Remove
rm(death_cancer) #Remove
ltab.data <-cbind(die[,1:2], total[,2]) 
detach(melanoma3)
attach(ltab.data)
lt=length(surv_mm)
surv_mm[lt+1]=NA
nevent = death_cancer
nlost = total[,2] - death_cancer
mytable <- lifetab(surv_mm, 5318, nlost, nevent)
mytable[,c(1:5, 8)]
detach(ltab.data)