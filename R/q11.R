## Date: 2015-03-04
## Purpose: To do the solution for Biostat III exercise 11 in R
## Author: Johan Zetterqvist
###############################################################################

## Install needed packages only need to be done once
## install.packages("survival")
## install.packages("dplyr")
## install.packages("readstata13")

###############################################################################
## Exercise 11
###############################################################################
## @knitr loadDependecies

require(survival)
require(dplyr)
require(readstata13)

## @knitr loadPreprocess

## Read melanoma data
## and select subcohorts
melanoma.l <-
  
  tbl_df( read.dta13("http://biostat3.net/download/melanoma.dta") ) %>%
  
  filter(stage=="Localised") %>%
  
  mutate(
    ## Create a death indicator
    death_cancer = as.numeric(status=="Dead: cancer"),
    death_any = as.numeric(status=="Dead: cancer" | status=="Dead: other") )


## Truncate follow-up time

melanoma.l2 <-
  
  mutate(melanoma.l,
         ## Create new death indicators (only count deaths within 120 months)
         death_cancer = death_cancer * as.numeric( surv_mm <= 120),
         death_any = death_any * as.numeric( surv_mm <= 120),
         ## Create a new time variable
         surv_mm = pmin(surv_mm, 120))

## @knitr 11.a

summary( coxfit11a <- coxph(Surv(surv_mm, death_any) ~ sex + year8594 + agegrp,
                           data = melanoma.l2,
                           ties = "breslow") )

## @knitr 11.b

summary( coxfit11b <- coxph(Surv(surv_mm, death_cancer) ~ sex + year8594 + agegrp,
                           data = melanoma.l2,
                           ties = "breslow") )

