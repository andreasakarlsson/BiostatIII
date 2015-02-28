# Purpose: To do the solution for Biostat III exercises in R
# Author: Annika Tillander, 2014-01-30
# Edited: Andreas Karlsson, 2015-02-24 
###############################################################################

## Install needed packages only need to be done once
## install.packages("survival")
## install.packages("KMsurv")
## install.packages("readstata13")
## install.packages("dplyr")
## install.packages("ggplot2")
## install.packages("Gally")
## install.packages("muhaz")
## install.packages("gridExtra")


###############################################################################
## Exercise 2
###############################################################################
## @knitr loadDependecies
require(readstata13) #Needed to read data set from Stata 13
require(survival) #for Surv and survfit
require(muhaz) #for hazard estimates
require(dplyr)
require(ggplot2)
require(GGally)
require(gridExtra)

## wrapper function for muhaz hazard smoother by strata
smoothHazard <- function(df, strat){
    tmp <- lapply(as.list(sort(levels(df[,strat]))),
                  function(x) c(strat=x, with(df[df[strat]==x,],
                                    muhaz(times=surv_mm, delta=death_cancer, max.time = max(surv_mm)))))
    out <- do.call("rbind", lapply(tmp,function(obj)
                                   data.frame(obj$haz.est, obj$est.grid, obj$strat)))
    colnames(out) <- c("Hazard","Time", strat)
    return(out)}

## @knitr loadPreprocess
melanoma_raw <- read.dta13("http://biostat3.net/download/melanoma.dta")
head(melanoma_raw[c("age", "sex", "stage", "status", "surv_mm", "surv_yy")])
melanoma <- mutate(melanoma_raw,
                   death_cancer = ifelse( status == "Dead: cancer", 1, 0),
                   death_all = ifelse( status == "Dead: cancer" | status == "Dead: other", 1, 0))

## @knitr a_tabulate
melanoma %>%
    group_by(stage) %>%
    summarise(Freq = n(), Percent = n()/dim(melanoma)[1]) %>%
    mutate(Cum = cumsum(Percent))

## @knitr a_plotSurv
mfit <- survfit(Surv(surv_mm, death_cancer) ~ stage, data = melanoma)
p1 <- ggsurv(mfit, cens.shape="|", cens.col = "black") +
    ggtitle("Kaplan-Meier survival estimates") + theme(legend.position="bottom")
p2 <- ggplot(smoothHazard(melanoma,"stage"), aes(x=Time, y=Hazard, colour= stage)) +
    geom_line() + ggtitle("Smoothed Hazard estimates") + theme(legend.position="bottom")
grid.arrange(p1, p2, ncol=2)

## @knitr b_crudeRates
melanoma %>%
    select(death_cancer, surv_mm, stage) %>%
    group_by(stage) %>%
    summarise(D = sum(death_cancer), M = sum(surv_mm), Rate = D/M) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

melanoma %>%
    select(death_cancer, surv_yy, stage) %>%
    group_by(stage) %>%
    summarise(D = sum(death_cancer), Y = sum(surv_yy), Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

## @knitr c_crudeRates1000 
melanoma %>%
    select(death_cancer, surv_yy, stage) %>%
    group_by(stage) %>%
    summarise(D = sum(death_cancer), Y = sum(surv_yy)/1000, Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

## @knitr d_crudeRates1000_sex
melanoma %>%
    select(death_cancer, surv_yy, sex) %>%
    group_by(sex) %>%
    summarise(D = sum(death_cancer), Y = sum(surv_yy)/1000, Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))

## @knitr d_plotSurv_sex
sfit <- survfit(Surv(surv_mm, death_cancer) ~ sex, data = melanoma)
p3 <- ggsurv(sfit, cens.shape="|", cens.col = "black") +
    ggtitle("Kaplan-Meier survival estimates") + theme(legend.position="bottom")
p4 <- ggplot(smoothHazard(melanoma, "sex"), aes(x=Time, y=Hazard, colour= sex)) +
    geom_line() + ggtitle("Smoothed Hazard estimates") + theme(legend.position="bottom")
grid.arrange(p3, p4, ncol=2)

## @knitr e_tabByAge
with(melanoma,table(status, agegrp))

## @knitr f_survStage
afit <- survfit(Surv(surv_mm, death_all) ~ stage, data = melanoma)
ggsurv(afit, cens.shape="|", cens.col = "black") + theme(legend.position="bottom") +
    ggtitle("Kaplan-Meier survival estimates\nAll-cause")

## @knitr g_allCa75p
mfit75 <- survfit(Surv(surv_mm, death_cancer) ~ stage, data = subset(melanoma,agegrp=="75+"))
afit75 <- survfit(Surv(surv_mm, death_all) ~ stage, data = subset(melanoma,agegrp=="75+"))
p6 <- ggsurv(mfit75, cens.shape="|", cens.col = "black") + theme(legend.position="bottom") +
    ggtitle("Kaplan-Meier survival estimates\nCancer | Age 75+")
p7 <- ggsurv(afit75, cens.shape="|", cens.col = "black") + theme(legend.position="bottom") +
    ggtitle("Kaplan-Meier survival estimates\nAll-cause | Age 75+")
grid.arrange(p6, p7, ncol=2)

## @knitr h_allCaAgeGrp 
mfitage <- survfit(Surv(surv_mm, death_cancer) ~ agegrp, data = melanoma)
afitage <- survfit(Surv(surv_mm, death_all) ~ agegrp, data = melanoma)
p8 <- ggsurv(mfitage, cens.shape="|", cens.col = "black") + theme(legend.position="bottom") +
    ggtitle("Kaplan-Meier estimates of\ncancer survival by age group")
p9 <- ggsurv(afitage, cens.shape="|", cens.col = "black") + theme(legend.position="bottom") +
    ggtitle("Kaplan-Meier estimates of\nall-cause survival by age group")
grid.arrange(p8, p9, ncol=2)

