# Date: 2014-01-30
# Purpose: To do the solution for Biostat III exercises 1-4 in R
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


###############################################################################
#Exercise 2
###############################################################################
require(foreign) #Needed to read data set from Stata
require(survival) #for Surv and survfit
require(muhaz) #for hazard estimates
require(lattice) #for grouped plots

#The data
melanoma <- read.dta("http://biostat3.net/download/melanoma.dta")
#a
#Construct a variable for starting time
start_time <- rep(0, length(melanoma$status))
#Recode the status variabel 
death_cancer <- NA 
death_cancer[melanoma$status =="Dead: cancer"] <- 1
death_cancer[melanoma$status =="Dead: other" | melanoma$status == "Alive" | melanoma$status == "Lost to follow-up"] <- 0
death_all <- NA 
death_all[melanoma$status =="Dead: cancer"|melanoma$status =="Dead: other"] <- 1
death_all[melanoma$status == "Alive" | melanoma$status == "Lost to follow-up"] <- 0
no_lost <- NA
no_lost[melanoma$status == "Lost to follow-up"] <- 1
no_lost[melanoma$status =="Dead: cancer"|melanoma$status =="Dead: other"|melanoma$status =="Alive"] <- 0
#Add the start time and cancer death to the data
melanoma2 <- data.frame(cbind(melanoma, death_cancer, death_all, no_lost, start_time))
#Kaplan-Meier estimate
mfit <- survfit(Surv(start_time, surv_mm, death_cancer==1)~stage, data = melanoma2)
plot(mfit, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of cause specific survival")
#Kaplan-Meier hazard non-smoothed
kfit <- kphaz.fit(time=melanoma2$surv_mm,status=melanoma2$death_cancer,strata=melanoma2$stage)
#Smoothed hazard
hfitU <- muhaz(times=melanoma2$surv_mm, delta=melanoma2$death_cancer, subset=melanoma2$stage=="Unknown")
hfitL <- muhaz(times=melanoma2$surv_mm, delta=melanoma2$death_cancer, subset=melanoma2$stage=="Localised")
hfitR <- muhaz(times=melanoma2$surv_mm, delta=melanoma2$death_cancer, subset=melanoma2$stage=="Regional")
hfitD <- muhaz(times=melanoma2$surv_mm, delta=melanoma2$death_cancer, subset=melanoma2$stage=="Distant")
#Combine the data for the different stages
hazard <- c(hfitU$haz.est, hfitL$haz.est, hfitR$haz.est, hfitD$haz.est)
time <- c(hfitU$est.grid, hfitL$est.grid, hfitR$est.grid, hfitD$est.grid)
grupp <- c(rep("Unknown", 101), rep("Localised", 101), rep("Regional", 101), rep("Distant", 101))
#Plot
xyplot(kfit$haz~kfit$time, groups=kfit$strata, type = "a", auto.key = list(points = FALSE, lines = TRUE), main="Kaplan-Meier Hazard estimates", xlab="Time since diagnosis in months", ylab="")
xyplot(hazard~time, groups=grupp, type = "a", auto.key = list(points = FALSE, lines = TRUE), main="Smoothed Hazard estimates", xlab="Time since diagnosis in months", ylab="")

#b
attach(melanoma2)
tall <- data.frame(death_cancer, surv_yy, surv_mm)
#Get the sums for the different stages
summa <- gsummary(tall, sum, groups=stage)
#Number of cancer deaths
D <- summa$death_cancer
#Total number of years
Y <- summa$surv_yy
#Total number of months
M <- summa$surv_mm
#Mortality rate for months
Rate_month <- round(D/M, 4)
#Confidence interval
Upper <- round(Rate_month+1.96*(Rate_month/sqrt(summa$death_cancer)), 4)
Lower <- round(Rate_month-1.96*(Rate_month/sqrt(summa$death_cancer)), 4)
Table_month <- cbind(D, M, Rate_month, Lower, Upper)
rownames(Table_month) <- rownames(summa)
Table_month
#Mortality rate for year
Rate_year <- round(D/Y, 4)
#Confidence interval
Upper <- round(Rate_year+1.96*(Rate_year/sqrt(summa$death_cancer)), 4)
Lower <- round(Rate_year-1.96*(Rate_year/sqrt(summa$death_cancer)), 4)
Table_year <- cbind(D, Y, Rate_year, Lower, Upper)
rownames(Table_year) <- rownames(summa)
Table_year
detach(melanoma2)

#c
attach(melanoma2)
tall <- data.frame(death_cancer, surv_mm)
#Get the sums for the different stages
summa <- gsummary(tall, sum, groups=stage)
#Number of cancer deaths
D <- summa$death_cancer
#Total number of months
Y <- round((summa$surv_mm/12)/1000, 3)
#Mortality rate for 1000 person year
Rate_year <- round(D/Y, 3)
#Confidence interval
Upper <- round(Rate_year+1.96*(Rate_year/sqrt(summa$death_cancer)), 3)
Lower <- round(Rate_year-1.96*(Rate_year/sqrt(summa$death_cancer)), 3)
Table_1000year <- cbind(D, Y, Rate_year, Lower, Upper)
rownames(Table_1000year) <- rownames(summa)
Table_1000year
detach(melanoma2)

#d
#Kaplan-Meier estimate
sfit <- survfit(Surv(start_time, surv_mm, death_cancer==1)~sex, data = melanoma2)
plot(sfit, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of cause specific survival")
#Kaplan-Meier hazard non-smoothed
kfit <- kphaz.fit(time=melanoma2$surv_mm,status=melanoma2$death_cancer,strata=melanoma2$sex)
#Smoothed hazard
hfitF <- muhaz(times=melanoma2$surv_mm, delta=melanoma2$death_cancer, subset=melanoma2$sex=="Female")
hfitM <- muhaz(times=melanoma2$surv_mm, delta=melanoma2$death_cancer, subset=melanoma2$sex=="Male")

#Combine the data for the different sex
hazard <- c(hfitF$haz.est, hfitM$haz.est)
time <- c(hfitF$est.grid, hfitM$est.grid)
grupp <- c(rep("Female", 101), rep("Male", 101))
#Plot
xyplot(kfit$haz~kfit$time, groups=kfit$strata, type = "a", auto.key = list(points = FALSE, lines = TRUE), main="Kaplan-Meier Hazard estimates", xlab="Time since diagnosis in months", ylab="")
xyplot(hazard~time, groups=grupp, type = "a", auto.key = list(points = FALSE, lines = TRUE), main="Smoothed Hazard estimates", xlab="Time since diagnosis in months", ylab="")

#e
attach(melanoma2)
#Table
table(status, agegrp)
detach(melanoma2)

#f 
#Kaplan-Meier estimate
afit <- survfit(Surv(start_time, surv_mm, death_all==1)~stage, data = melanoma2)
par(mfrow=c(1,2))
plot(mfit, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of cause specific survival")
plot(afit, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of non specific survival")

#g
#Kaplan-Meier estimate
mfit75 <- survfit(Surv(start_time, surv_mm, death_cancer==1)~stage, data = melanoma2[melanoma2$agegrp=="75+",])
afit75 <- survfit(Surv(start_time, surv_mm, death_all==1)~stage, data = melanoma2[melanoma2$agegrp=="75+",])
par(mfrow=c(1,2))
plot(mfit75, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of cause specific survival")
plot(afit75, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of non specific survival")

#h
#Kaplan-Meier estimate
mfitage <- survfit(Surv(start_time, surv_mm, death_cancer==1)~agegrp, data = melanoma2)
afitage <- survfit(Surv(start_time, surv_mm, death_all==1)~agegrp, data = melanoma2)
par(mfrow=c(1,2))
plot(mfitage, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of cause specific survival")
plot(afitage, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of non specific survival")


###############################################################################
#Exercise 3
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
no_lost <- NA
no_lost[melanoma_sub$status == "Lost to follow-up"] <- 1
no_lost[melanoma_sub$status =="Dead: cancer"|melanoma_sub$status =="Dead: other"|melanoma_sub$status =="Alive"] <- 0
#Add the start time and cancer death to the data
melanoma3 <- data.frame(cbind(melanoma_sub, death_cancer, death_all, no_lost, start_time))
#a
#Kaplan-Meier estimate
mfityear8594 <- survfit(Surv(start_time, surv_mm, death_cancer==1)~year8594, data = melanoma3)
plot(mfityear8594, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of cause specific survival")

#b
#Kaplan-Meier hazard non-smoothed
kfityear8594 <- kphaz.fit(time=melanoma3$surv_mm,status=melanoma3$death_cancer,strata=melanoma3$year8594)
#min and max time
min(melanoma3$surv_mm[melanoma3$year8594=="Diagnosed 85-94"])
max(melanoma3$surv_mm[melanoma3$year8594=="Diagnosed 85-94"])
min(melanoma3$surv_mm)
max(melanoma3$surv_mm)
#Smoothed hazard
hfit7584 <- muhaz(times=melanoma3$surv_mm, delta=melanoma3$death_cancer, subset=melanoma3$year8594=="Diagnosed 75-84", min.time=0.5, max.time=251.5)
hfit8594 <- muhaz(times=melanoma3$surv_mm, delta=melanoma3$death_cancer, subset=melanoma3$year8594=="Diagnosed 85-94", min.time=0.5, max.time=131.5)
#Combine the data for the different diagnosis year
hazard <- c(hfit7584$haz.est, hfit8594$haz.est)
time <- c(hfit7584$est.grid, hfit8594$est.grid)
grupp <- c(rep("Diagnosed 75-84", 101), rep("Diagnosed 85-94", 101))
#Plot
xyplot(kfityear8594$haz~kfityear8594$time, groups=kfityear8594$strata, type = "a", auto.key = list(points = FALSE, lines = TRUE), main="Kaplan-Meier Hazard estimates", xlab="Time since diagnosis in months", ylab="")
xyplot(hazard~time, groups=grupp, type = "a", auto.key = list(points = FALSE, lines = TRUE), main="Smoothed Hazard estimates", xlab="Time since diagnosis in months", ylab="")

#c
#Test
survdiff(Surv(surv_mm, death_cancer==1) ~ year8594, data=melanoma3)
#Equivalent to the Peto & Peto modfication of the Gehan-Wilcoxon test
survdiff(Surv(surv_mm, death_cancer==1) ~ year8594, data=melanoma3, rho=1)

#d
attach(melanoma3)
tall <- data.frame(death_cancer, surv_mm)
#Get the sums for the different stages
summa <- gsummary(tall, sum, groups=agegrp)
#Number of cancer deaths
D <- summa$death_cancer
#Total number of months
Y <- round(summa$surv_mm/1000, 3)
#Mortality rate for 1000 person year
Rate_year <- round(D/Y, 3)
#Confidence interval
Upper <- round(Rate_year+1.96*(Rate_year/sqrt(summa$death_cancer)), 3)
Lower <- round(Rate_year-1.96*(Rate_year/sqrt(summa$death_cancer)), 3)
Table_1000months <- cbind(D, Y, Rate_year, Lower, Upper)
rownames(Table_1000months) <- rownames(summa)
Table_1000months
detach(melanoma3)

#Kaplan-Meier estimate
mfitagegrp <- survfit(Surv(start_time, surv_mm, death_cancer==1)~agegrp, data = melanoma3)
plot(mfitagegrp, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of cause specific survival")

#e
attach(melanoma3)
tall <- data.frame(death_cancer, surv_mm)
#Get the sums for the different stages
summa <- gsummary(tall, sum, groups=agegrp)
#Number of cancer deaths
D <- summa$death_cancer
#Total number of months
Y <- round((summa$surv_mm/12)/1000, 3)
#Mortality rate for 1000 person year
Rate_year <- round(D/Y, 3)
#Confidence interval
Upper <- round(Rate_year+1.96*(Rate_year/sqrt(summa$death_cancer)), 3)
Lower <- round(Rate_year-1.96*(Rate_year/sqrt(summa$death_cancer)), 3)
Table_1000years <- cbind(D, Y, Rate_year, Lower, Upper)
rownames(Table_1000years) <- rownames(summa)
Table_1000years
detach(melanoma3)

#f
#Kaplan-Meier estimate
mfitsex <- survfit(Surv(start_time, surv_mm, death_cancer==1)~sex, data = melanoma3)
plot(mfitsex, xlab="Time since diagnosis in months", ylab="S(t)", main="Kaplan-Meier estimates of cause specific survival")

#Kaplan-Meier hazard non-smoothed
kfitsex <- kphaz.fit(time=melanoma3$surv_mm,status=melanoma3$death_cancer,strata=melanoma3$sex)

#Smoothed hazard
hfitFemale <- muhaz(times=melanoma3$surv_mm, delta=melanoma3$death_cancer, subset=melanoma3$sex=="Female")
hfitMale <- muhaz(times=melanoma3$surv_mm, delta=melanoma3$death_cancer, subset=melanoma3$sex=="Male")
#Combine the data for the different diagnosis year
hazard <- c(hfitFemale$haz.est, hfitMale$haz.est)
time <- c(hfitFemale$est.grid, hfitMale$est.grid)
grupp <- c(rep("Female", 101), rep("Male", 101))
#Plot
xyplot(kfitsex$haz~kfitsex$time, groups=kfitsex$strata, type = "a", auto.key = list(points = FALSE, lines = TRUE), main="Kaplan-Meier Hazard estimates", xlab="Time since diagnosis in months", ylab="")
xyplot(hazard~time, groups=grupp, type = "a", auto.key = list(points = FALSE, lines = TRUE), main="Smoothed Hazard estimates", xlab="Time since diagnosis in months", ylab="")

#Test
survdiff(Surv(surv_mm, death_cancer==1) ~ sex, data=melanoma3)
#Equivalent to the Peto & Peto modfication of the Gehan-Wilcoxon test
survdiff(Surv(surv_mm, death_cancer==1) ~ sex, data=melanoma3, rho=1)

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