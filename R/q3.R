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

