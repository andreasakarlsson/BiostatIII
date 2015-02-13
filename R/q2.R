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


