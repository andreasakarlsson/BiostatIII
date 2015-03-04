## Purpose: To do the solution for Biostat III exercises in R
## Author: Annika Tillander, 2014-01-30
## Edited: Andreas Karlsson, 2015-02-28
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
## install.packages("RCurl")


###############################################################################
## Exercise 3
###############################################################################
## @knitr loadDependecies
require(readstata13) #Needed to read data set from Stata 13
require(survival) #for Surv and survfit
require(muhaz) #for hazard estimates
require(dplyr)
require(ggplot2)
require(GGally)
require(gridExtra)
require(RCurl)
helperFunctions <- eval(expr=parse(text=getURL("https://raw.githubusercontent.com/andreasakarlsson/BiostatIII/master/R/helpFuncBioIII.R", .opts = list(ssl.verifypeer = FALSE))))

## @knitr loadPreprocess
melanoma_raw<- read.dta13("http://biostat3.net/download/melanoma.dta")
melanoma <- melanoma_raw %>%
    filter(stage=="Localised") %>%
    mutate(death_cancer = ifelse( status == "Dead: cancer", 1, 0),
           death_all = ifelse( status == "Dead: cancer" | status == "Dead: other", 1, 0))

## @knitr a_survDiaDate
mfityear8594 <- survfit(Surv(surv_mm, death_cancer==1) ~ year8594, data = melanoma)
ggsurv(mfityear8594, cens.shape="|", cens.col = "black") +
    ggtitle("Kaplan-Meier survival estimates")


## @knitr b_hazDiaDate
ggplot(smoothHazard(melanoma,"year8594"), aes(x=Time, y=Hazard, colour= year8594)) +
    geom_line() + ggtitle("Smoothed Hazard estimates")

## @knitr c_testDiaDate
## Log-rank test for equality of survivor functions
survdiff(Surv(surv_mm, death_cancer==1) ~ year8594, data=melanoma)
## Equivalent to the Peto & Peto modfication of the Gehan-Wilcoxon test
survdiff(Surv(surv_mm, death_cancer==1) ~ year8594, data=melanoma, rho=1)

## @knitr d_crudeRates1000_agegrp
melanoma %>%
    select(death_cancer, surv_mm, agegrp) %>%
    group_by(agegrp) %>%
    summarise(D = sum(death_cancer), Y = sum(surv_mm)/1000, Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))


mfit_agegrp <- survfit(Surv(surv_mm, death_cancer==1) ~ agegrp, data = melanoma)
ggsurv(mfit_agegrp, cens.shape="|", cens.col = "black") +
    ggtitle("Kaplan-Meier survival estimates") + xlab("Months since diagnosis")

## @knitr e_crudeRates1000_agegrp
mfit_agegrp_year<- survfit(Surv(surv_mm/12, death_cancer==1) ~ agegrp, data = melanoma)
ggsurv(mfit_agegrp_year, cens.shape="|", cens.col = "black") + 
    ggtitle("Kaplan-Meier survival estimates") + xlab("Years since diagnosis") 

melanoma %>%
    select(death_cancer, surv_mm, agegrp) %>%
    group_by(agegrp) %>%
    summarise(D = sum(death_cancer), Y = sum(surv_mm)/12/1000, Rate = D/Y) %>%
    mutate(CI_low = Rate + qnorm(0.025) * Rate / sqrt(D),
           CI_high = Rate + qnorm(0.975) * Rate / sqrt(D))


## @knitr f_sexDiff
mfit_sex <- survfit(Surv(surv_mm, death_cancer) ~ sex, data = melanoma)
p1 <- ggsurv(mfit_sex, cens.shape="|", cens.col = "black") +
    ggtitle("Kaplan-Meier survival estimates") + theme(legend.position="bottom")
p2 <- ggplot(smoothHazard(melanoma,"sex"), aes(x=Time, y=Hazard, colour= sex)) +
    geom_line() + ggtitle("Smoothed Hazard estimates") + theme(legend.position="bottom")
grid.arrange(p1, p2, ncol=2)


## Log-rank test for equality of survivor functions
survdiff(Surv(surv_mm, death_cancer==1) ~ sex, data=melanoma)
