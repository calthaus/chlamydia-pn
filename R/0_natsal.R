# Modelling partner notification for chlamydia
# Christian L. Althaus, 24 May 2019

# Derive sexual behavior parameters from Natsal-3

library(stats4)
library(mvtnorm)

# Load Natsal-3 data set
load("../data/natsal3urine.RData")

# Age group
age_lo <- 16
age_up <- 34

# Subset of data
urine <- within(subset(natsal3,
					subset = ct_posconfirmed %in% c('negative', 'positive') & !(hetnonew %in% c(-1, 995, 999)) & dage >= age_lo & dage <= age_up & !is.na(urine_wt), 
					select = c('ct_posconfirmed', 'urine_wt', 'hetnonew')),
					ct <- as.integer(ifelse(ct_posconfirmed == 'positive', 1, 0)))

# MLE of two Poisson processes
f <- function(x,a,m1,m2) {
    (1-a)*dpois(x,m1)+a*dpois(x,m2)
}
nll <- function(a,m1,m2) {
	a <- plogis(a)
	m1 <- exp(m1)
	m2 <- exp(m2)
    -sum(urine$urine_wt*log(f(urine$hetnonew,a,m1,m2)))
}

fit <- mle(nll,start=list(a=qlogis(0.05),m1=log(0.5),m2=log(10)))
saveRDS(fit, file = "../out/fit_natsal.rds")
fit_natsal <- readRDS("../out/fit_natsal.rds")
