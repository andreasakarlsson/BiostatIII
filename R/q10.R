## Purpose: To do the solution for Biostat III exercises in R
## Author: Andreas Karlsson, 2015-03-02
###############################################################################

## Install needed packages only need to be done once
## install.packages("survival")
## install.packages("readstata13")
## install.packages("dplyr")
## install.packages("ggplot2")
## install.packages("muhaz")

###############################################################################
## Exercise 10
###############################################################################
## @knitr loadDependecies
require(readstata13) #Needed to read data set from Stata 13
require(survival) #for Surv and survfit
require(muhaz) #for hazard estimates
require(dplyr)
require(ggplot2)


## @knitr loadPreprocess
melanoma_raw<- read.dta13("http://biostat3.net/download/melanoma.dta")
melanoma <- melanoma_raw %>%
    filter(stage == "Localised") %>%
    mutate(death_cancer = ifelse( status == "Dead: cancer" & surv_mm <= 120, 1, 0), #censuring for > 120 monts
           trunk_yy = ifelse(surv_mm <=  120, surv_mm/12, 10))  #scale to years and trunkate to 10 years 

## @knitr 10.a
# Using muhaz to smooth the kaplan-meier hazards by strata the time and bandwidth options were selected based on smoother performance. There must be better way of doing this. Sorry about the warnings.
hazDiaDate <- melanoma %>%  group_by(year8594) %>%
    do( h = muhaz(times = .$trunk_yy, delta = .$death_cancer, min.time = min(.$trunk_yy[.$death_cancer==1]), max.time = max(.$trunk_yy[.$death_cancer==1]), bw.method="g", bw.grid=5)) %>%
    do( data.frame(Hazard = .$h$haz.est, Time = .$h$est.grid, Strata = .$year8594))

## Max hazard ratio
maxHaz <- hazDiaDate %>% group_by(Strata) %>%
    summarise(stratMax=max(Hazard))
print(maxHaz)
maxHaz$stratMax[2]/maxHaz$stratMax[1]

## Comparing hazards
ggplot(hazDiaDate, aes(x=Time, y=Hazard, colour= Strata)) +
    geom_line() + ggtitle("Smoothed Hazard estimates")

## @knitr 10.b
## Comparing hazards on a log scales
ggplot(hazDiaDate, aes(x=Time, y=Hazard, colour= Strata)) +
    geom_line() + ggtitle("Smoothed Hazard estimates") + scale_y_continuous(trans='log')


## @knitr 10.c
## Calculating cumulative hazard per strata
hazDiaDate <- group_by(hazDiaDate, Strata) %>%
    mutate(cumHaz=cumsum(Hazard))

ggplot(hazDiaDate, aes(log(Time), -log(cumHaz), color=Strata)) + geom_line() + geom_point()




## @knitr 10.d
# Cox regression with time-since-entry as the timescale
# Note that R uses the Efron method for approximating the likelihood in the presence
# whereas Stata (and most other software) use the Breslow method
cox1 <- coxph(Surv(trunk_yy, death_cancer==1) ~ year8594, method=c("breslow"), data=melanoma)
summary(cox1)


## @knitr 10.e

cox2 <- coxph(Surv(trunk_yy, death_cancer==1) ~ sex + year8594 + agegrp, method=c("breslow"), data=melanoma)
summary(cox2)

## Plot of the scaled Schoenfeld residuals for calendar period 1985–94.
## The smooth line shows the estimated log hazard ratio as a function of time.
cox2.phtest <- cox.zph(cox2, transform="identity") #Stata appears to be using 'identity'
plot(cox2.phtest[2],resid=TRUE, se=TRUE, main="Schoenfeld residuals", ylim=c(-4,4))

## @knitr 10.g
## The results from the previous proportional hazards assumption test
print(cox2.phtest)

## @knitr 10.h
melanoma2p8Split <- survSplit(melanoma, cut=c(2), end="trunk_yy", start="start", event="death_cancer")
melanoma2p8Split <- mutate(melanoma2p8Split, fu = as.factor(start))

##Tabulate ageband including risk_time
melanoma2p8Split %>% select(id, start, trunk_yy) %>% filter(id<=3) %>% arrange(id, trunk_yy)

head(melanoma2p8Split)

cox2p8Split1 <- coxph(Surv(start, trunk_yy, death_cancer) ~ sex + year8594 + agegrp*fu, method=c("breslow"), data=melanoma2p8Split)
summary(cox2p8Split1)


## @knitr 10.i
cox2p8Split2 <- coxph(Surv(start, trunk_yy, death_cancer) ~ sex + year8594 + fu + fu:agegrp, method=c("breslow"), data=melanoma2p8Split)
summary(cox2p8Split2)


## @knitr 10.hother
## Ref: http://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
melanoma2p8 <- tmerge(melanoma, melanoma,  id=id, tstart=rep(0, nrow(melanoma)), tstop=trunk_yy,  death = event(death_cancer), after2 = tdc(rep(2, nrow(melanoma))))

head(melanoma2p8)

## This is not look correct!? Shocking!
cox2p8 <- coxph(Surv(tstart, tstop, death) ~ sex + year8594 + agegrp + agegrp:after2  + cluster(id), method=c("breslow"), data=melanoma2p8)
summary(cox2p8)

## @knitr 10.iother
## This is not look correct!? Shocking!
cox2p8_2<- coxph(Surv(tstart, tstop, death) ~  sex + year8594 + after2 + after2:agegrp + cluster(id), method=c("breslow"), data=melanoma2p8)
summary(cox2p8_2)


## cox2p8Split <- coxph(Surv(start, trunk_yy, death_cancer) ~ sex + year8594 + agegrp + agegrp:fu, method=c("breslow"), data=melanoma2p8Split)
## ## I think stata is using another reference category
## summary(cox2p8Split)


