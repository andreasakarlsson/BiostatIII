# Date: 2014-02-04
# Purpose: To do the solution for Biostat III exercise 9 in R
# Author: Annika Tillander
###############################################################################

#Load the needed packages needs to be done each session
require(survival)
require(KMsurv)
require(foreign) #Needed to read data set from Stata
require(nlme) #for gsummary
require(muhaz) #for hazards
###############################################################################
#Exercise 9
###############################################################################

#The data
melanoma <- read.dta("http://biostat3.net/download/melanoma.dta")
melanoma_sub <- melanoma[melanoma$stage=="Localised",]

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


# Patients who survived more than 10 years are censored at 10 years
melanoma3$death_cancer[melanoma3$surv_mm>120] <- 0
melanoma3$exit[melanoma3$surv_mm>120] <- melanoma3$dx[melanoma3$surv_mm>120] + 10
melanoma3$surv_yy[melanoma3$surv_mm>120] <- 10
melanoma3$surv_mm[melanoma3$surv_mm>120] <- 120

# Cox regression with time-since-entry as the timescale
# Note that R uses the Efron method for approximating the likelihood in the presence
# whereas Stata (and most other software) use the Breslow method
cox1 <- coxph(Surv(surv_mm, death_cancer==1) ~ 1 + year8594, method=c("breslow"), data=melanoma3)
summary(cox1)

#a
#b

#c
#Cox regression
cox2 <- coxph(Surv(surv_mm, death_cancer==1) ~ 1 + sex + year8594 + agegrp, method=c("breslow"), data=melanoma3)
summary(cox2)
#iii ???? I have not found a Wald test for the global effect of age only how to test for each category of age ...
#Test for each covariate
cox.zph(cox2)
#Tests if there is a difference between two or more survival curves
survdiff(Surv(surv_mm, death_cancer==1) ~ 1 + agegrp, data=melanoma3)
survdiff(Surv(surv_mm, death_cancer==1) ~ 1 + sex + year8594 + agegrp, data=melanoma3)
survdiff(Surv(surv_mm, death_cancer==1) ~ 1 + sex + year8594, data=melanoma3)

#d
#The model without agegrp
cox3 <- coxph(Surv(surv_mm, death_cancer==1) ~ 1 + sex + year8594, method=c("breslow"), data=melanoma3)
#Likelihood ratio test for overall effect of age
anova(cox2, cox3)

#e
#Given a survival data set and a set of specified cut times, split each record into multiple subrecords at each cut time. The new data set will be in ‘counting process’ format, with a start time, stop time, and event status for each record.
melanoma4 <- survSplit(melanoma3, cut=c(0, 12*(1:9)), end="surv_mm", start="start_time", event="death_cancer", episode="i")
melanoma4$start_time <- as.factor(melanoma4$start_time)
#Poisson regression
pos_fit <- glm(death_cancer ~ sex + year8594 + agegrp + start_time, family=poisson, data=melanoma4)
summary(pos_fit)


#Table for comparison between Cox and Poisson
hazratio_pois <- exp(coef(pos_fit))
#Confidence interval
CI_p <- exp(confint(pos_fit))
CI_c <- exp(confint(cox2))
cox_tab <- cbind(summary(cox2)$coefficients[,c(2,4,5)], CI_c)
poisson_tab <- cbind(hazratio_pois, summary(pos_fit)$coefficients[,3:4], CI_p)
round(rbind(cox_tab, poisson_tab[2:6,]), 3)



###############################################################################
#Exercise 10
###############################################################################

#The data
melanoma <- read.dta("http://biostat3.net/download/melanoma.dta")
melanoma_sub <- melanoma[melanoma$stage=="Localised",]

#Construct a variable for starting time
start_time <- rep(0, length(melanoma_sub$status))
#Recode the status variabel 
death_cancer <- NA 
death_cancer[melanoma_sub$status =="Dead: cancer"] <- 1
death_cancer[melanoma_sub$status =="Dead: other" | melanoma_sub$status == "Alive" | melanoma_sub$status == "Lost to follow-up"] <- 0
death_all <- NA 
death_all[melanoma_sub$status =="Dead: cancer"|melanoma_sub$status =="Dead: other"] <- 1
death_all[melanoma_sub$status == "Alive" | melanoma_sub$status == "Lost to follow-up"] <- 0
# Patients who survived more than 10 years are censored at 10 years
death_cancer[melanoma_sub$surv_mm>120] <- 0
melanoma_sub$exit[melanoma_sub$surv_mm>120] <- melanoma_sub$dx[melanoma_sub$surv_mm>120] + 10
melanoma_sub$surv_yy[melanoma_sub$surv_mm>120] <- 10
melanoma_sub$surv_mm[melanoma_sub$surv_mm>120] <- 120
surv_my <- melanoma_sub$surv_mm/12
#Add the start time and cancer death to the data
melanoma3 <- data.frame(cbind(melanoma_sub, death_cancer, death_all, surv_my, start_time))

#a
#min and max time
min(melanoma3$surv_my[melanoma3$year8594=="Diagnosed 85-94"])
max(melanoma3$surv_my[melanoma3$year8594=="Diagnosed 85-94"])
min(melanoma3$surv_my)
max(melanoma3$surv_my)
#Smoothed hazard
hfit75 <- muhaz(times=melanoma3$surv_my, delta=melanoma3$death_cancer, subset=melanoma3$year8594=="Diagnosed 75-84", min.time=min(melanoma3$surv_my), max.time=10)
hfit85 <- muhaz(times=melanoma3$surv_my, delta=melanoma3$death_cancer, subset=melanoma3$year8594=="Diagnosed 85-94", min.time=min(melanoma3$surv_my), max.time=10)

#Combine the data for the different stages
hazard <- c(hfit75$haz.est, hfit85$haz.est)
time <- c(hfit75$est.grid, hfit85$est.grid)
grupp <- c(rep("Diagnosed 75-84", 101), rep("Diagnosed 85-94", 101))
#Plot
require(lattice)
xyplot(hazard~time, groups=grupp, type = "a", auto.key = list(points = FALSE, lines = TRUE), main="Smoothed Hazard estimates", xlab="Time since diagnosis in years", ylab="")

#b
#min and max time
min(log(melanoma3$surv_my[melanoma3$year8594=="Diagnosed 85-94"]))
max(log(melanoma3$surv_my[melanoma3$year8594=="Diagnosed 85-94"]))
min(log(melanoma3$surv_my))
max(log(melanoma3$surv_my))
#Smoothed hazard
hfit75 <- muhaz(times=log(melanoma3$surv_my), delta=melanoma3$death_cancer, subset=melanoma3$year8594=="Diagnosed 75-84", min.time=min(log(melanoma3$surv_my)), max.time=log(10))
hfit85 <- muhaz(times=melanoma3$surv_my, delta=melanoma3$death_cancer, subset=melanoma3$year8594=="Diagnosed 85-94", min.time=min(melanoma3$surv_my), max.time=10)

plot(hfit75)
xyplot(log(hazard)~time, groups=grupp, type = "a", auto.key = list(points = FALSE, lines = TRUE), main="Smoothed Hazard estimates", xlab="Time since diagnosis in years", ylab="")


