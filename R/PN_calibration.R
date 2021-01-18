# Heterosexual chlamydia transmission model with partner notification
# Christian L. Althaus, 30 September 2018

# Model calibration using Approximate Bayesian Computation (ABC)

rm(list=ls())

library(deSolve)
library(stats4)
library(mvtnorm)

source('PN_model.R')
fit_natsal <- readRDS("fit_natsal.rds")

# Catch command from bash
args <- commandArgs(trailingOnly = TRUE)
batch <- as.numeric(args[1])
filename <- paste0("sets",batch,".RData")

# Create parameter sets using prior distributions
set.seed(862445+batch)
nsets <- 2e4
sets <- as.data.frame(matrix(data = NA, nrow = nsets, ncol = 40))
names(sets) <- c("N11","N21","N12","N22","c11","c21","c12","c22","b11","b21","b12","b22",
				 "epsilon","delta","f1","f2","sigma1","sigma2","gamma1","gamma2","tau1","tau2","f_P1","f_P2","treat1","treat2","omega1","omega2",
				 "pre11","pre21","pre12","pre22","pre1","pre2","post11","post21","post12","post22","post1","post2")

m <- coef(fit_natsal)
sigma <- vcov(fit_natsal)
pars <- data.frame(rmvnorm(nsets, mean = m, sigma = sigma))
pars[, 1] <- plogis(pars[, 1])
pars[, 2:3] <- exp(pars[, 2:3])

sets$N11 <- 1 - pars[, 1]
sets$N21 <- 1 - pars[, 1]
sets$N12 <- pars[, 1]
sets$N22 <- pars[, 1]
sets$c11 <- pars[, 2]
sets$c21 <- pars[, 2]
sets$c12 <- pars[, 3]
sets$c22 <- pars[, 3]
sets$b12 <- runif(nsets,0,1)
sets$b22 <- runif(nsets,0,1)
sets$b11 <- runif(nsets,sets$b12,1)
sets$b21 <- runif(nsets,sets$b22,1)
sets$epsilon <- runif(nsets,0,1)
sets$delta <- rep(1/19,nsets)
sets$f1 <- runif(nsets,0,1)
sets$f2 <- runif(nsets,0,1)
sets$sigma1 <- 365/runif(nsets,14,42)
sets$sigma2 <- 365/runif(nsets,14,42)
sets$gamma1 <- rgamma(nsets, 100, 100/0.74)
sets$gamma2 <- rgamma(nsets, 2, 4.8)
sets$tau1 <- rgamma(nsets, 1, 1/0.3)
sets$tau2 <- rgamma(nsets, 1, 1/0.1)
sets$f_P1 <- rbeta(nsets, 5, 5*(1 - 0.53)/0.53)
sets$f_P2 <- rbeta(nsets, 5, 5*(1 - 0.53)/0.53)
sets$treat1 <- rep(365/3.2,nsets)
sets$treat2 <- rep(365/3.2,nsets)
sets$omega1 <- 1/rgamma(nsets, 2, 2/5)
sets$omega2 <- 1/rgamma(nsets, 2, 2/5)

# Run simulations
times  <- c(0,1e2)
for(i in 1:nsets) {
    # Set parameters
	temp <- as.list(sets[i,1:28])
	parms <- list(
		N = matrix(c(temp$N11,temp$N21,temp$N12,temp$N22),nrow=2),
		c = matrix(c(temp$c11,temp$c21,temp$c12,temp$c22),nrow=2),
		epsilon = temp$epsilon,
		delta = temp$delta,
		b = matrix(c(temp$b11,temp$b21,temp$b12,temp$b22),nrow=2),
		f = c(temp$f1,temp$f2),
		sigma = c(temp$sigma1,temp$sigma2),
		gamma = c(temp$gamma1,temp$gamma2),
		tau = c(temp$tau1,temp$tau2),
		f_P = c(temp$f_P1,temp$f_P2),
		treat = c(temp$treat1,temp$treat2),
		omega = c(temp$omega1,temp$omega2))
	parms <- c(parms,list(rho=mixing(parms$N[1,],parms$c[1,],parms$epsilon)))
	parms <- c(parms,list(beta=array(c(parms$b[1,1],parms$b[2,1],sqrt(prod(parms$b[1,])),sqrt(prod(parms$b[2,])),sqrt(prod(parms$b[1,])),sqrt(prod(parms$b[2,])),parms$b[1,2],parms$b[2,2]),c(2,2,2))))

    # Set initial variables
    prev <- 0.03 # Asymptomatic female/male only
    init <- rep(0,24)
    init[1:4] <- (1-prev)*parms$N
    init[9:12] <- prev*parms$N
    
    # Run simulations
    parms_orig <- parms
    parms$tau <- c(0,0)
    parms$f_P <- c(0,0)
    simulation <- as.data.frame(ode(init, times, model, parms))

    sets[i,"pre1"] <- sum(simulation[2,c(5,9,13,17,7,11,15,19)+1])
    sets[i,"pre11"] <- sum(simulation[2,c(5,9,13,17)+1])/parms$N[1,1]
    sets[i,"pre12"] <- sum(simulation[2,c(7,11,15,19)+1])/parms$N[1,2]
    
    sets[i,"pre2"] <- sum(simulation[2,c(6,10,14,18,8,12,16,20)+1])
    sets[i,"pre21"] <- sum(simulation[2,c(6,10,14,18)+1])/parms$N[2,1]
    sets[i,"pre22"] <- sum(simulation[2,c(8,12,16,20)+1])/parms$N[2,2]
    
    init <- as.numeric(simulation[2,-1])
    parms <- parms_orig
    simulation <- as.data.frame(ode(init, times, model, parms))
    
    sets[i,"post1"] <- sum(simulation[2,c(5,9,13,17,7,11,15,19)+1])
    sets[i,"post11"] <- sum(simulation[2,c(5,9,13,17)+1])/parms$N[1,1]
    sets[i,"post12"] <- sum(simulation[2,c(7,11,15,19)+1])/parms$N[1,2]
    
    sets[i,"post2"] <- sum(simulation[2,c(6,10,14,18,8,12,16,20)+1])
    sets[i,"post21"] <- sum(simulation[2,c(6,10,14,18)+1])/parms$N[2,1]
    sets[i,"post22"] <- sum(simulation[2,c(8,12,16,20)+1])/parms$N[2,2]
    
	if(i%%(nsets/10) == 0) save(sets,file=filename)
}

q(save="no")
