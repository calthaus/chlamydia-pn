# Modelling partner notification for chlamydia
# Christian L. Althaus, 21 May 2019

# Create parameter subsets from UBELIX simulations

set.seed(75247)

# Load full data set
sim_name <- "50x2e4_sets_v4"
sim_length <- 50
load(paste0("../data/", sim_name, "/sets1.RData"))
all_sets <- sets
for(i in 2:sim_length) {
	load(paste0("../data/", sim_name, "/sets", i, ".RData"))
	all_sets <- rbind(all_sets,sets)
}
sets <- all_sets
rm(all_sets)
save(sets,file="../out/sets.RData")

# Chlamydia prevalence for 16-34 year olds in Natsal-3
load("../data/natsal3urine.RData")

age_lo <- 16
age_up <- 34

prev <- matrix(NA, 2, 3)

dct <- within(subset(natsal3, !is.na(urine_wt) & urintested == 'yes' & dage >= age_lo & dage <= age_up & rsex == "Female", 
			select = c('ct_posconfirmed', 'urine_wt')), 
		ct <- as.integer(ifelse(ct_posconfirmed == 'positive', 1, 0)))
prev[1, 1] <- weighted.mean(dct$ct, dct$urine_wt)
p <- numeric(1e4)
for(i in 1:length(p)) {
	p[i] <- mean(sample(dct$ct,length(dct$ct),replace=TRUE,prob=dct$urine_wt/max(dct$urine_wt)))
}
prev[1, 2:3] <- quantile(p, c(0.025, 0.975))

dct <- within(subset(natsal3, !is.na(urine_wt) & urintested == 'yes' & dage >= age_lo & dage <= age_up & rsex == "Male", 
					 select = c('ct_posconfirmed', 'urine_wt')), 
			  ct <- as.integer(ifelse(ct_posconfirmed == 'positive', 1, 0)))
prev[2, 1] <- weighted.mean(dct$ct,dct$urine_wt)
p <- numeric(1e4)
for(i in 1:length(p)) {
	p[i] <- mean(sample(dct$ct,length(dct$ct),replace=TRUE,prob=dct$urine_wt/max(dct$urine_wt)))
}
prev[2, 2:3] <- quantile(p, c(0.025, 0.975))

save(prev, file = "../out/prevalence.RData")

# Calibration

nsets <- numeric(3)

# post1.RData: Post prevalence overall
w <- which(sets$post1 > prev[1, 2]
		& sets$post1 < prev[1, 3]
		& sets$post2 > prev[2, 2]
		& sets$post2 < prev[2, 3])
post <- sets[w,]
nsets[1] <- dim(post)[1]
save(post,file="../out/posterior1.RData")

# post2.RData: Pre and post prevalence overall
w <- which(sets$pre1 < 2*prev[1, 1]
		& sets$pre2 < 2*prev[2, 1]
		& sets$post1 > prev[1, 2]
		& sets$post1 < prev[1, 3]
		& sets$post2 > prev[2, 2]
		& sets$post2 < prev[2, 3])
post <- sets[w,]
nsets[2] <- dim(post)[1]
save(post,file="../out/posterior2.RData")

# post3.RData: Pre prevalence overall, post prevalence per activity group and overall
w <- which(sets$pre1 < 2*prev[1, 1]
		& sets$pre2 < 2*prev[2, 1]
		& sets$post1 > prev[1, 2]
		& sets$post1 < prev[1, 3]
		& sets$post2 > prev[2, 2]
		& sets$post2 < prev[2, 3]
		& sets$post12 < 0.15
		& sets$post22 < 0.15)
post <- sets[w,]
nsets[3] <- dim(post)[1]
save(post,file="../out/posterior3.RData")

write.table(nsets, file = "../out/nsets.csv")

rm(list = ls())
