# Modelling partner notification for chlamydia
# Christian L. Althaus, 21 May 2019

# Plot posterior distributions

set.seed(27545)

# Load libraries
library(plotrix)

# Load and shorten full parameter set
load("../out/sets.RData")
sets <- sets[1:1e4, ]

# Load chlamydia prevalence data
load("../out/prevalence.RData")

# Loop over all parameter subsets
for(i in 1:3) {
	load(paste0("../out/posterior", i, ".RData"))

	# Plot prior and posterior prevalence
	pdf(file = paste0("../figures/prevalence", i, ".pdf"), width = 9, height = 6)
	par(mfrow = c(2, 3))
	plot(density(post$post1),xlim=c(0,0.05),main="Female (overall)",xlab="Chlamydia prevalence",frame=FALSE)
	polygon(density(post$post1),border=NA,col=rgb(1, 0, 0, alpha=0.3))
	lines(density(post$pre1))
	polygon(density(post$pre1),border=NA,col=rgb(1, 0, 0, alpha=0.1))
	plotCI(prev[1, 1], max(density(post$post1)$y)/10, li = prev[1, 2], ui = prev[1, 3], pch = 19, err = "x", add = TRUE)
	plot(density(post$post11),xlim=c(0,0.05),main="Female (low activity)",xlab="Chlamydia prevalence",frame=FALSE)
	polygon(density(post$post11),border=NA,col=rgb(1, 0, 0, alpha=0.3))
	lines(density(post$pre11))
	polygon(density(post$pre11),border=NA,col=rgb(1, 0, 0, alpha=0.1))
	plot(density(post$post12),xlim=c(0,0.5),main="Female (high activity)",xlab="Chlamydia prevalence",frame=FALSE)
	polygon(density(post$post12),border=NA,col=rgb(1, 0, 0, alpha=0.3))
	lines(density(post$pre12))
	polygon(density(post$pre12),border=NA,col=rgb(1, 0, 0, alpha=0.1))
	plot(density(post$post2),xlim=c(0,0.05),main="Male (overall)",xlab="Chlamydia prevalence",frame=FALSE)
	polygon(density(post$post2),border=NA,col=rgb(0, 0, 1, alpha=0.3))
	lines(density(post$pre2))
	polygon(density(post$pre2),border=NA,col=rgb(0, 0, 1, alpha=0.1))
	plotCI(prev[2, 1], max(density(post$post2)$y)/10, li = prev[2, 2], ui = prev[2, 3], pch = 19, err = "x", add = TRUE)
	plot(density(post$post21),xlim=c(0,0.05),main="Male (low activity)",xlab="Chlamydia prevalence",frame=FALSE)
	polygon(density(post$post21),border=NA,col=rgb(0, 0, 1, alpha=0.3))
	lines(density(post$pre21))
	polygon(density(post$pre21),border=NA,col=rgb(0, 0, 1, alpha=0.1))
	plot(density(post$post22),xlim=c(0,0.5),main="Male (high activity)",xlab="Chlamydia prevalence",frame=FALSE)
	polygon(density(post$post22),border=NA,col=rgb(0, 0, 1, alpha=0.3))
	lines(density(post$pre22))
	polygon(density(post$pre22),border=NA,col=rgb(0, 0, 1, alpha=0.1))
	dev.off()
	
	# Plot prior and posterior of parameters
	par_sets <- sets
	par_post <- post
	par_sets["delta"] <- 1/par_sets["delta"]
	par_sets["sigma1"] <- 1/par_sets["sigma1"]*365
	par_sets["gamma1"] <- 1/par_sets["gamma1"]*365
	par_sets["treat1"] <- 1/par_sets["treat1"]*365
	par_sets["omega1"] <- 1/par_sets["omega1"]
	par_sets["sigma2"] <- 1/par_sets["sigma2"]*365
	par_sets["gamma2"] <- 1/par_sets["gamma2"]*365
	par_sets["treat2"] <- 1/par_sets["treat2"]*365
	par_sets["omega2"] <- 1/par_sets["omega2"]
	par_post["delta"] <- 1/par_post["delta"]
	par_post["sigma1"] <- 1/par_post["sigma1"]*365
	par_post["gamma1"] <- 1/par_post["gamma1"]*365
	par_post["treat1"] <- 1/par_post["treat1"]*365
	par_post["omega1"] <- 1/par_post["omega1"]
	par_post["sigma2"] <- 1/par_post["sigma2"]*365
	par_post["gamma2"] <- 1/par_post["gamma2"]*365
	par_post["treat2"] <- 1/par_post["treat2"]*365
	par_post["omega2"] <- 1/par_post["omega2"]
	
	pdf(file = paste0("../figures/posterior", i, ".pdf"), width = 15, height = 8)
	par(mfcol = c(3, 5))
	name <- "Proportion high sexual activity"
	j <- 3
	y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = NA, main = name, frame = FALSE,
		 sub = paste0(q1[2], " (", q1[1], "-", q1[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Heterosexual partner change rate (low)"
	j <- 5
	y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = "Per year", main = name, frame = FALSE,
		 sub = paste0(q1[2], " (", q1[1], "-", q1[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Heterosexual partner change rate (high)"
	j <- 7
	y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = "Per year", main = name, frame = FALSE,
		 sub = paste0(q1[2], " (", q1[1], "-", q1[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Transmission probability (low-low)"
	j <- 9
	y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = NA, main = name, frame = FALSE,
		 sub = paste0("f: ", q1[2], " (", q1[1], "-", q1[3], "), m: ", q2[2], " (", q2[1], "-", q2[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 0, alpha = 0.2))
	lines(density(par_post[, j + 1]))
	polygon(density(par_post[, j + 1]), border = NA, col = rgb(0, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Transmission probability (low-high)"
	j <- 9
	y_lim <- max(c(density(sqrt(par_post[, j]*par_post[, j + 2]))$y, density(sqrt(par_post[, j + 1]*par_post[, j + 3]))$y))
	q1 <- signif(quantile(sqrt(par_post[, j]*par_post[, j + 2]))[2:4], 2)
	q2 <- signif(quantile(sqrt(par_post[, j + 1]*par_post[, j + 3]))[2:4], 2)
	plot(density(sqrt(par_post[, j]*par_post[, j + 2])),
		 xlim = range(sqrt(par_post[, j]*par_post[, j + 2])), ylim = c(0, y_lim), xlab = NA, main = name, frame = FALSE,
		 sub = paste0("f: ", q1[2], " (", q1[1], "-", q1[3], "), m: ", q2[2], " (", q2[1], "-", q2[3], ")"))
	polygon(density(sqrt(par_post[, j]*par_post[, j + 2])), border = NA, col = rgb(1, 0, 0, alpha = 0.2))
	lines(density(sqrt(par_post[, j + 1]*par_post[, j + 3])))
	polygon(density(sqrt(par_post[, j + 1]*par_post[, j + 3])), border = NA, col = rgb(0, 0, 1, alpha = 0.2))
	lines(density(sqrt(par_sets[, j]*par_sets[, j + 2])))
	polygon(density(sqrt(par_sets[, j + 1]*par_sets[, j + 3])), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Transmission probability (high-high)"
	j <- 11
	y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = NA, main = name, frame = FALSE,
		 sub = paste0("f: ", q1[2], " (", q1[1], "-", q1[3], "), m: ", q2[2], " (", q2[1], "-", q2[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 0, alpha = 0.2))
	lines(density(par_post[, j + 1]))
	polygon(density(par_post[, j + 1]), border = NA, col = rgb(0, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Sexual mixing coefficient"
	j <- 13
	y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = NA, main = name, frame = FALSE,
		 sub = paste0(q1[2], " (", q1[1], "-", q1[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Duration of immunity"
	j <- 27
	y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = "Years", main = name, frame = FALSE,
		 sub = paste0("f: ", q1[2], " (", q1[1], "-", q1[3], "), m: ", q2[2], " (", q2[1], "-", q2[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 0, alpha = 0.2))
	lines(density(par_post[, j + 1]))
	polygon(density(par_post[, j + 1]), border = NA, col = rgb(0, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Frequency symptomatic"
	j <- 15
	y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = NA, main = name, frame = FALSE,
		 sub = paste0("f: ", q1[2], " (", q1[1], "-", q1[3], "), m: ", q2[2], " (", q2[1], "-", q2[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 0, alpha = 0.2))
	lines(density(par_post[, j + 1]))
	polygon(density(par_post[, j + 1]), border = NA, col = rgb(0, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Symptomatic duration"
	j <- 17
	y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = "Days", main = name, frame = FALSE,
		 sub = paste0("f: ", q1[2], " (", q1[1], "-", q1[3], "), m: ", q2[2], " (", q2[1], "-", q2[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 0, alpha = 0.2))
	lines(density(par_post[, j + 1]))
	polygon(density(par_post[, j + 1]), border = NA, col = rgb(0, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Asymptomatic duration (females)"
	j <- 19
	y_lim <- max(density(par_post[, j])$y)
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = "Days", main = name, frame = FALSE,
		 sub = paste0("f: ", q1[2], " (", q1[1], "-", q1[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 0, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Asymptomatic duration (males)"
	j <- 20
	y_lim <- max(density(par_post[, j])$y)
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = c(0, 2e3), ylim = c(0, y_lim), xlab = "Days", main = name, frame = FALSE,
		 sub = paste0("m: ", q1[2], " (", q1[1], "-", q1[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(0, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j], to = 2e3))
	polygon(density(par_sets[, j], to = 2e3), border = NA, col = rgb(0, 1, 0, alpha = 0.2))	
	
	name <- "Treatment rate (females)"
	j <- 21
	y_lim <- max(density(par_post[, j])$y)
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = c(0, 1), ylim = c(0, y_lim), xlab = "Per year", main = name, frame = FALSE,
		 sub = paste0("f: ", q1[2], " (", q1[1], "-", q1[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 0, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	name <- "Treatment rate (males)"
	j <- 22
	y_lim <- max(density(par_post[, j])$y)
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = c(0, 1), ylim = c(0, y_lim), xlab = "Per year", main = name, frame = FALSE,
		 sub = paste0("m: ", q1[2], " (", q1[1], "-", q1[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(0, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))	
	
	name <- "Treated partners"
	j <- 23
	y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	q1 <- signif(quantile(par_post[, j])[2:4], 2)
	q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	plot(density(par_post[, j]),
		 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = NA, main = name, frame = FALSE,
		 sub = paste0("f: ", q1[2], " (", q1[1], "-", q1[3], "), m: ", q2[2], " (", q2[1], "-", q2[3], ")"))
	polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 0, alpha = 0.2))
	lines(density(par_post[, j + 1]))
	polygon(density(par_post[, j + 1]), border = NA, col = rgb(0, 0, 1, alpha = 0.2))
	lines(density(par_sets[, j]))
	polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	
	# name <- "Time to partner treatment"
	# j <- 25
	# y_lim <- max(c(density(par_post[, j])$y, density(par_post[, j + 1])$y))
	# q1 <- signif(quantile(par_post[, j])[2:4], 2)
	# q2 <- signif(quantile(par_post[, j + 1])[2:4], 2)
	# plot(density(par_post[, j]),
	# 	 xlim = range(par_sets[, j]), ylim = c(0, y_lim), xlab = "Days", main = name, frame = FALSE,
	# 	 sub = paste0("f: ", q1[2], " (", q1[1], "-", q1[3], "), m: ", q2[2], " (", q2[1], "-", q2[3], ")"))
	# polygon(density(par_post[, j]), border = NA, col = rgb(1, 0, 0, alpha = 0.2))
	# lines(density(par_post[, j + 1]))
	# polygon(density(par_post[, j + 1]), border = NA, col = rgb(0, 0, 1, alpha = 0.2))
	# lines(density(par_sets[, j]))
	# polygon(density(par_sets[, j]), border = NA, col = rgb(0, 1, 0, alpha = 0.2))
	dev.off()
}

rm(list = ls())
