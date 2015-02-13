# A first attempt at analysing the localised melanoma data
#
# Please help!
#
# Paul Dickman, November 2010
#
# setwd("c:/survival/")

library( survival )
library(foreign)
library(Epi)

# localised (stage==1) melanoma diagnosed in Finland
melanoma <- subset.data.frame(read.dta("http://biostat3.net/download/melanoma.dta",convert.factors=FALSE), stage==1)

# Lexis requires time variables to be numeric and have the same units
# we will use years
melanoma$dx<-1970+as.numeric(melanoma$dx)/365.24
melanoma$exit<-1970+as.numeric(melanoma$exit)/365.24
melanoma$bdate<-1970+as.numeric(melanoma$bdate)/365.24

head( melanoma)

# Patients who survived more than 10 years are censored at 10 years
melanoma$status <- ifelse(melanoma$surv_mm>120,0,melanoma$status)
melanoma$exit <- ifelse(melanoma$surv_mm>120,melanoma$dx + 10,melanoma$exit)
melanoma$surv_yy <- ifelse(melanoma$surv_mm>120,10,melanoma$surv_yy)
melanoma$surv_mm<- ifelse(melanoma$surv_mm>120,120,melanoma$surv_mm)

# Kaplan-Meir estimates of survival for each period
sf <- survfit( Surv( surv_mm, status==1 ) ~ 1 + year8594, data=melanoma )
plot(sf)

# Cox regression with time-since-entry as the timescale
# Note that R uses the Efron method for approximating the likelihood in the presence
# whereas Stata (and most other software) use the Breslow method
cox1 <- coxph( Surv( surv_mm, status==1 ) ~ 1 + year8594, method=c("breslow"), data=melanoma )
summary( cox1 )

# Cox regression again but using different variables to define time
cox2 <- coxph( Surv( exit-dx, status==1 ) ~ 1 + year8594, method=c("breslow"), data=melanoma )
summary( cox2 )

# Poisson regression assuming constant mortality over time
poisson1 <- glm( status==1 ~ 1 + year8594 + offset( log( surv_mm ) ), family=poisson, data=melanoma )
summary( poisson1 )

# Define a Lexis object
L <- Lexis( entry = list( per=dx, attage=dx-bdate, fu=0 ), exit = list( per=exit ),
 exit.status = factor( status, labels=c("alive","dead:cancer","dead:other","lost") ),
 data = melanoma )
str( L )
sample <- L[sample(5), ]

# Default plot of Lexis object
plot( sample )

# With a grid and deaths as endpoints
plot( sample, grid=0:10*10, col="black" )
points( sample, pch=c(NA,16)[sample$lex.Xst] )

# With a lot of bells and whistles:
plot( sample, grid=0:20*5, col="black", xaxs="i", yaxs="i",
      xlim=c(1975,2000), ylim=c(60,100), lwd=3, las=1 )
points( sample, pch=c(NA,16)[sample$lex.Xst], col="red", cex=1.5 )

# Split follow-up time into annual intervals
head( L, 3)
split <- splitLexis(L, breaks=seq(0,10,1), time.scale="fu")
head( split, 8 )

# Rerun the Poisson regression model assuming constant hazards on the original data
poisson1b <- glm( status==1 ~ 1 + year8594 + offset( log(exit-dx) ), family=poisson, data=melanoma )
summary( poisson1b )

# Now rerun the Poisson regression model assuming constant hazards on the split data (results should be the same)
poisson2 <- glm( as.numeric(split$lex.Xst)==2 ~ 1 + year8594 + offset( log(lex.dur) ), family=poisson, data=split )
summary( poisson2 )

# Now adjust for follow-up time
poisson3 <- glm( as.numeric(split$lex.Xst)==2 ~ 1 + as.factor(fu) + year8594 + offset( log(lex.dur) ), family=poisson, data=split )
summary( poisson3 )







