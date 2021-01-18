# Modelling partner notification for chlamydia
# Christian L. Althaus, 25 May 2019

# Create model output for manuscript

set.seed(93512)

# Load libraries and functions
library(doMC)
library(deSolve)
library(xtable)
source("PN_model.R")

# Function to create ouptut for each scenario
results <- function(scenario) {
	nsets <- 1e3
	load(paste0("../out/posterior", scenario, ".RData"))
	post <- post[1:nsets, ]
	
	### 1. Prevalence
	prev <- apply(post[, 29:40], 2, quantile, probs = c(0.025, 0.5, 0.975))
	prev <- t(round(100*prev, 2))
	write.table(prev, file = paste0("../out/prevalence", scenario, ".csv"), sep = "\t\t")
	
	### 2. WNW matrix
	results_who <- array(NA, c(2, 4, 7, nsets))
	
	### 3. APT interventions
	par_range1 <- seq(0, 0.30, 0.05)
	par_range2 <- seq(3.2, 0.2, -1)
	results_apt1 <- array(NA, c(2, 3, nsets, length(par_range1)))
	results_apt2 <- array(NA, c(2, 3, nsets, length(par_range2)))
	
	for(i in 1:nsets) {
		# Set parameters
		temp <- as.list(post[i, 1:28])
		parms <- list(
			N = matrix(c(temp$N11, temp$N21, temp$N12, temp$N22), nrow = 2),
			c = matrix(c(temp$c11, temp$c21, temp$c12, temp$c22), nrow = 2),
			epsilon = temp$epsilon,
			delta = temp$delta,
			b = matrix(c(temp$b11, temp$b21, temp$b12, temp$b22), nrow = 2),
			f = c(temp$f1, temp$f2),
			sigma = c(temp$sigma1, temp$sigma2),
			gamma = c(temp$gamma1, temp$gamma2),
			tau = c(temp$tau1, temp$tau2),
			f_P = c(temp$f_P1, temp$f_P2),
			treat = c(temp$treat1, temp$treat2),
			omega = c(temp$omega1, temp$omega2))
		parms <- c(parms, list(rho = mixing(parms$N[1, ], parms$c[1, ], parms$epsilon)))
		parms <- c(parms, list(beta = array(c(parms$b[1, 1], parms$b[2, 1], sqrt(prod(parms$b[1, ])), sqrt(prod(parms$b[2,])), sqrt(prod(parms$b[1, ])), sqrt(prod(parms$b[2, ])), parms$b[1, 2], parms$b[2, 2]), c(2, 2, 2))))
		
		# Set initial variables
		prev <- 0.03
		init <- rep(0, 24)
		init[1:4] <- (1 - prev)*parms$N
		init[9:12] <- prev*parms$N
		
		# Run simulation
		times  <- c(0, 1e2)
		simulation <- as.data.frame(ode(init, times, model, parms))
		init <- as.numeric(simulation[2, -1])
		
		# Calculate WNW matrix
		for(j in 1:2) {
			who <- notify(j, init, parms)[, c(1:4, 9, 10)]
			who <- cbind(who, who[, 5]*who[, 6]/sum(who[, 5]*who[, 6]))
			who[, 6] <- who[, 6]/sum(who[, 6])
			results_who[j, , , i] <- who
		}
		
		# Simulate intervention 1 (increase in number of treated partners)
		times  <- c(0,5)
		base_f_P <- parms$f_P
		for(j in 1:length(par_range1)) {
			parms$f_P <- base_f_P + par_range1[j]
			if(parms$f_P[1] > 1) parms$f_P[1] <- 1
			if(parms$f_P[2] > 1) parms$f_P[2] <- 1
			
			simulation <- as.data.frame(ode(init, times, model, parms))
			
			results_apt1[1, 1, i, j] <- 1 - sum(simulation[2,c(5,9,13,17,7,11,15,19)+1])/post[i,"post1"]
			results_apt1[1, 2, i, j] <- 1 - sum(simulation[2,c(5,9,13,17)+1])/parms$N[1,1]/post[i,"post11"]
			results_apt1[1, 3, i, j] <- 1 - sum(simulation[2,c(7,11,15,19)+1])/parms$N[1,2]/post[i,"post12"]
			
			results_apt1[2, 1, i, j] <- 1 - sum(simulation[2,c(6,10,14,18,8,12,16,20)+1])/post[i,"post2"]
			results_apt1[2, 2, i, j] <- 1 - sum(simulation[2,c(6,10,14,18)+1])/parms$N[2,1]/post[i,"post21"]
			results_apt1[2, 3, i, j] <- 1 - sum(simulation[2,c(8,12,16,20)+1])/parms$N[2,2]/post[i,"post22"]
		}
		
		# Simulate intervention 2 (reduction in time to partner treatment)
		times  <- c(0,5)
		parms$f_P <- base_f_P
		for(j in 1:length(par_range2)) {
			parms$treat <- rep(365/par_range2[j], 2)
			
			simulation <- as.data.frame(ode(init, times, model, parms))
			
			results_apt2[1, 1, i, j] <- 1 - sum(simulation[2,c(5,9,13,17,7,11,15,19)+1])/post[i,"post1"]
			results_apt2[1, 2, i, j] <- 1 - sum(simulation[2,c(5,9,13,17)+1])/parms$N[1,1]/post[i,"post11"]
			results_apt2[1, 3, i, j] <- 1 - sum(simulation[2,c(7,11,15,19)+1])/parms$N[1,2]/post[i,"post12"]
			
			results_apt2[2, 1, i, j] <- 1 - sum(simulation[2,c(6,10,14,18,8,12,16,20)+1])/post[i,"post2"]
			results_apt2[2, 2, i, j] <- 1 - sum(simulation[2,c(6,10,14,18)+1])/parms$N[2,1]/post[i,"post21"]
			results_apt2[2, 3, i, j] <- 1 - sum(simulation[2,c(8,12,16,20)+1])/parms$N[2,2]/post[i,"post22"]
		}
	}
	
	# Create Who-Notified-Whom (WNW) matrix
	for(i in 1:2) {
		m <- array(NA, c(4, 7))
		index <- ifelse(i == 1, "female", "male")
		partner <- ifelse(i == 1, "male", "female")
		for(j in 1:4) {
			for(k in 1:7) {
				q <- quantile(results_who[i, j, k, ], probs = c(0.25, 0.5, 0.75))
				m[j, k] <- paste0(round(100*q[2], 0), "% (", round(100*q[1], 0), "%--", round(100*q[3], 0), "%)")	
			}
		}
		m <- t(m)
		m <- rbind(c(rep(c("Low", "High"), 2)), m)
		m <- rbind(c("Symptomatic", "", "Asymptomatic", ""), m)
		m <- rbind(c(paste0("Index (", index, ")"), rep("", 3)), m)
		m <- cbind(c(rep("", 3), "Low", "High", "Low", "High", "", "", ""), m)
		m <- cbind(c("", "", paste0("Partner (", partner, ")"), "Symptomatic", "", "Asymptomatic", "", "Chlamydia positivity in partner", "Proportion of notified partners", "Proportion positive among notified"), m)
		m <- rbind(m[1:7, ], rep("", 6), m[8:10, ])
		m <- xtable(m,
					caption = paste0("Who-Notifies-Whom (WNW) matrix for ", index, " index cases (scenario ", scenario, "). Numbers are given as median and interquartile range."),
					label = paste0("wnw_", index, scenario),
					align = "lll|cccc")
		m <- print(m, hline.after = 3, size = "scriptsize", include.rownames = FALSE, include.colnames = FALSE)
		write(m, file = paste0("../out/wnw_", index, scenario, ".tex"))
	}
	
	# Plot effect of APT intervention
	pdf(file = paste0("../figures/apt", scenario, ".pdf"), width = 8, height = 8)
	par(mfrow = c(3, 2))
	for(group in 1:3) {
		if(group == 1) main_label = "Population overall"
		if(group == 2) main_label = "Sexual activity group (low)"
		if(group == 3) main_label = "Sexual activity group (high)"
		
		boxplot(1e2*results_apt1[1, group, , ],
				at = 1:length(par_range1) - 0.2, boxwex = 0.3,
				ylim = c(0, 100), col = rgb(1, 0, 0, alpha = 0.3),
				xlab = "Increase in number of treated partners (%)", ylab = "Relative reduction in prevalence (%)",
				main = main_label,
				axes = FALSE, frame = FALSE)
		boxplot(1e2*results_apt1[2, group, , ],
				at = 1:length(par_range1) + 0.2, boxwex = 0.3, 
				col = rgb(0, 0, 1, alpha = 0.3),
				axes = FALSE, frame = FALSE, add = TRUE)
		axis(1, 1:length(par_range1), 1e2*par_range1)
		axis(2)
		
		boxplot(1e2*results_apt2[1, group, , ],
				at = 1:length(par_range2) - 0.2, boxwex = 0.3,
				ylim = c(0, 5), col = rgb(1, 0, 0, alpha = 0.3),
				xlab = "Time to partner treatment (days)", ylab = "Relative reduction in prevalence (%)",
				main = main_label,
				axes = FALSE, frame = FALSE)
		boxplot(1e2*results_apt2[2, group, , ],
				at = 1:length(par_range2) + 0.2, boxwex = 0.3, 
				col = rgb(0, 0, 1, alpha = 0.3),
				axes = FALSE, frame = FALSE, add = TRUE)
		axis(1, 1:length(par_range2), par_range2)
		axis(2)
	}
	dev.off()
	
	# Store results of APT interventions
	saveRDS(results_apt1, file = paste0("../out/scenario", scenario, "_apt1.rds"))
	saveRDS(results_apt2, file = paste0("../out/scenario", scenario, "_apt2.rds"))
}

# Loop over all scenarios
registerDoMC(3)
foreach(scenario = 1:3) %dopar% results(scenario)

rm(list = ls())

q(save = "no")
