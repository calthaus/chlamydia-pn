# Modelling partner notification for chlamydia
# Christian L. Althaus, 21 May 2019

# Definition of functions

# Sexual mixing matrix
mixing <- function(N,c,epsilon) {
    max <- length(N)
    cN <- sum(c*N)
    f <- matrix(0,max,max)
    for(i in 1:max) f[i,i] <- epsilon
    rho <- matrix(NA,max,max)
    for(i in 1:max) for(j in 1:max) rho[i,j] <- f[i,j] + (1-epsilon)*c[j]*N[j]/cN
    rho
}

# Partner notification matrix
notify <- function(ii, x, parms) { # ii is index case
	with(as.list(c(parms,x)),{
		# Assign variables
		S	<- matrix(x[1:4],nrow=2)
		I_S <- matrix(x[5:8],nrow=2)
		I_A <- matrix(x[9:12],nrow=2)
		P_S	<- matrix(x[13:16],nrow=2)
		P_A	<- matrix(x[17:20],nrow=2)
		R	<- matrix(x[21:24],nrow=2)
		
		# Define partner of index case
		i <- ii%%2+1
		
		notification <- matrix(0,nrow=4,ncol=8)
		index <- c(I_S[ii,],I_A[ii,])
		# Rates at wich index cases ii will show up for treatment
		treatment1 <- c(sigma[ii]+tau[ii],sigma[ii]+tau[ii],gamma[ii]+tau[ii],gamma[ii]+tau[ii])
		treatment2 <- c(sigma[ii]+tau[ii],sigma[ii]+tau[ii],tau[ii],tau[ii])
		# Probability of clearance in partners i
		clear <- function(interval) {
			c(exp(-(sigma[i]+tau[i])*interval), exp(-(sigma[i]+tau[i])*interval), exp(-(gamma[i]+tau[i])*interval), exp(-(gamma[i]+tau[i])*interval),
			  exp(-(sigma[i]+tau[i]+treat[i])*interval), exp(-(sigma[i]+tau[i]+treat[i])*interval), exp(-(gamma[i]+tau[i]+treat[i])*interval), exp(-(gamma[i]+tau[i]+treat[i])*interval))
		}
		
		# Loop over all index cases
		for(j in 1:4) {
			# Set activity class of index case
			jj <- (j+1)%%2+1
			# Determine probability that most recent contact was with infector
			p <- exp(-c[ii,jj]/(treatment1[j]))
			# Determine type of infector
			v <- beta[i,,jj]*rho[jj,]*c(I_S[i,],I_A[i,],P_S[i,],P_A[i,])/N[i,]
			v <- v/sum(v)
			notification[j,] <- v*p*clear(1/treatment1[j])
			# Determine type of infectee
			v <- c(f[i]*beta[ii,jj,]*rho[jj,]*S[i,]/N[i,],(1-f[i])*beta[ii,jj,]*rho[jj,]*S[i,]/N[i,],0,0,0,0)
			# Determine type of (infected) random
			v <- v + rho[jj,]*c(I_S[i,],I_A[i,],P_S[i,],P_A[i,])/N[i,]
			notification[j,] <- notification[j,] + v*(1-p)*clear(1/c[ii,jj])
		}
		notification <- cbind(notification,rowSums(notification))
		notification <- cbind(notification,index*treatment2*f_P[ii])
		return(notification)
	})
}

# Transmission model
model <- function(t, x, parms) {
    with(as.list(c(parms,x)),{
        # Assign variables
        S	<- matrix(x[1:4],nrow=2)
        I_S <- matrix(x[5:8],nrow=2)
        I_A <- matrix(x[9:12],nrow=2)
        P_S	<- matrix(x[13:16],nrow=2)
        P_A	<- matrix(x[17:20],nrow=2)
        R	<- matrix(x[21:24],nrow=2)

        # Define derivatives
        dS 	    <- array(0,c(2,2))
        dI_S 	<- array(0,c(2,2))
        dI_A 	<- array(0,c(2,2))
        dP_S 	<- array(0,c(2,2))
        dP_A 	<- array(0,c(2,2))
        dR 	    <- array(0,c(2,2))
        
        # Loop over sex (1: female; 2: male)
        for(i in 1:2) {
            ii <- i%%2+1
            # # Calculate PN matrix (PN only starts from I, but not from P)
            # #			1		2		3		4		5		6		7		8
            # #			I_S1	I_S2	I_A1	I_A3	P_S1	P_S2	P_A1	P_A2
            # # 1	I_S1
            # # 2	I_S2
            # # 3	I_A1
            # # 4	I_A2
            notification <- matrix(0,nrow=4,ncol=8)
            index <- c(I_S[ii,],I_A[ii,])
            # Rates at wich index cases ii will show up for treatment
            treatment1 <- c(sigma[ii]+tau[ii],sigma[ii]+tau[ii],gamma[ii]+tau[ii],gamma[ii]+tau[ii])
            treatment2 <- c(sigma[ii]+tau[ii],sigma[ii]+tau[ii],tau[ii],tau[ii])
            # Probability of clearance in partners i
            clear <- function(interval) {
                c(exp(-(sigma[i]+tau[i])*interval), exp(-(sigma[i]+tau[i])*interval), exp(-(gamma[i]+tau[i])*interval), exp(-(gamma[i]+tau[i])*interval),
                  exp(-(sigma[i]+tau[i]+treat[i])*interval), exp(-(sigma[i]+tau[i]+treat[i])*interval), exp(-(gamma[i]+tau[i]+treat[i])*interval), exp(-(gamma[i]+tau[i]+treat[i])*interval))
            }

            # Loop over all index cases
            for(j in 1:4) {
                # Set activity class of index case
                jj <- (j+1)%%2+1
                # Determine probability that most recent contact was with infector
                p <- exp(-c[ii,jj]/(treatment1[j]))
                # Determine type of infector
                v <- beta[i,,jj]*rho[jj,]*c(I_S[i,],I_A[i,],P_S[i,],P_A[i,])/N[i,]
                v <- v/sum(v)
                notification[j,] <- v*p*clear(1/treatment1[j])
                # Determine type of infectee
                v <- c(f[i]*beta[ii,jj,]*rho[jj,]*S[i,]/N[i,],(1-f[i])*beta[ii,jj,]*rho[jj,]*S[i,]/N[i,],0,0,0,0)
                # Determine type of (infected) random
                v <- v + rho[jj,]*c(I_S[i,],I_A[i,],P_S[i,],P_A[i,])/N[i,]
                notification[j,] <- notification[j,] + v*(1-p)*clear(1/c[ii,jj])
                notification[j,] <- notification[j,]*index[j]*treatment2[j]*f_P[ii]
            }
            notification <- colSums(notification)
            
            # Loop activity class (1: low; 2: high)
            for(j in 1:2) {
                # Calculate transmission term
                trans  <- c[i,j]*S[i,j]*sum(beta[ii,,j]*rho[j,]*(I_S[ii,]+I_A[ii,]+P_S[ii,]+P_A[ii,])/N[ii,])
                # Calculate derivatives
                dS[i,j]   <- delta*N[i,j] - trans + sigma[i]*(I_S[i,j]+P_S[i,j]) + tau[i]*(I_S[i,j]+I_A[i,j]+P_S[i,j]+P_A[i,j]) + treat[i]*(P_S[i,j]+P_A[i,j]) + omega[i]*R[i,j] - delta*S[i,j]
                dI_S[i,j] <- f[i]*trans - sigma[i]*I_S[i,j] - tau[i]*I_S[i,j] - notification[j] - delta*I_S[i,j]
                dI_A[i,j] <- (1-f[i])*trans - gamma[i]*I_A[i,j] - tau[i]*I_A[i,j] - notification[j+2] - delta*I_A[i,j]
                dP_S[i,j] <- notification[j] - treat[i]*P_S[i,j] - sigma[i]*P_S[i,j] - tau[i]*P_S[i,j] - delta*P_S[i,j]
                dP_A[i,j] <- notification[j+2] - treat[i]*P_A[i,j] - gamma[i]*P_A[i,j] - tau[i]*P_A[i,j] - delta*P_A[i,j]
                dR[i,j]   <- gamma[i]*(I_A[i,j]+P_A[i,j]) - omega[i]*R[i,j] - delta*R[i,j]
            }
        }
        
        # Return results
        return(list(c(dS,dI_S,dI_A,dP_S,dP_A,dR)))
    })
}	
