# Fit piecewise exponential models using Poisson method - fixed effect model
# Cut points at 6 and 12 months

##################################################################################

# Start with a blank environment
rm(list=ls())

# Load libraries
library(coda)
library(survival)
library(R2WinBUGS)

# Set the root directory 
path <- 'R:/SCFreeman/scf20/Melanoma_NMA/Piecewise/Fixed_6_12'
setwd(path)

# Load data
data <- read.csv("../melanoma_aggregate_data.csv")

# Order trials by trialid
data <- data[order(data$trialid, data$treatment, data$spgrp),]

# Set the location for WinBUGS
bugs.directory <- "Z:\\winbugs143_unrestricted\\winbugs14_full_patched\\WinBUGS14"

# WinBUGS burn-in & simulation size
num.sims <- 30000
burn.in <- 30000

#-----------------------------------------------------------------------------
# Data
#-----------------------------------------------------------------------------

# Create treatment indicator variables - caution needed to esnure consistency equations hold
data$trt2 <- 0
data$trt2[data$treatment=="TRL"] <- 1

data$trt3 <- 0
data$trt3[data$treatment=="IPI+DTIC"] <- 1

data$trt7 <- 0
data$trt7[data$treatment=="VM+COB"] <- 1

data$trt10 <- 0
data$trt10[data$treatment=="NIV+IPI"] <- 1

data$trt11 <- 0
data$trt11[data$treatment=="PEM"] <- 1

data$trt12 <- 0
data$trt12[data$treatment=="IPI+SRG"] <- 1

data$trt13 <- 0
data$trt13[data$treatment=="SEL+DTIC"] <- 1

data$trt4 <- 0
data$trt4[data$treatment=="DB" & data$trialid==1] <- 1
data$trt4[data$treatment=="DB" & data$trialid==8] <- -1

data$trt5 <- 0
data$trt5[data$treatment=="DB+TR" & data$trialid==7] <- 1
data$trt5[data$treatment=="DB+TR" & data$trialid==8] <- -1

data$trt6 <- 0
data$trt6[data$treatment=="VM" & (data$trialid==2 | data$trialid==8)] <- 1

data$trt8 <- 0
data$trt8[data$treatment=="NIV" & data$trialid==4] <- 1

data$trt9 <- 0
data$trt9[data$treatment=="NIV" & data$trialid==3] <- 1


# Calculate indicator variables for before/after cut point
# Cut point at 6 months (so 6 months included as before)
data$z1[data$time<=6] <- 0
data$z1[data$time>6] <- 1

# Cut point at 12 months (so 12 months included as before)
data$z4[data$time<=12] <- 0
data$z4[data$time>12] <- 1

# ints=no. of time intervals
ints <- max(data$spgrp)

# J=no. of studies
J <- length(unique(data$trialid))

# ntrt=no. of trts (only needed in WinBUGS data if calculating survival)
ntrt <- 13

# N=rows of data
N <- nrow(data)

nma_data <- list(spgrp=data$spgrp, y=data$y, d=data$nevents, trial=data$trialid, trt2=data$trt2,
                 trt3=data$trt3, trt4=data$trt4, trt5=data$trt5, trt6=data$trt6, trt7=data$trt7,
                 trt8=data$trt8, trt9=data$trt9, trt10=data$trt10, trt11=data$trt11, 
                 trt12=data$trt12, trt13=data$trt13, z1=data$z1, 
                 z4=data$z4,
                 N=N, ints=ints, J=J, ntrt=ntrt)


#-----------------------------------------------------------------------------
# Model
#-----------------------------------------------------------------------------

# Create initial values for model
# alpha is a matrix (no. of trials by no. of treatment contrast parameters)
alpha1 <- c(rep(0.1, ntrt-1))
alpha2 <- c(rep(0.2, ntrt-1))
alpha3 <- c(rep(-0.1, ntrt-1))

# beta is a matrix (J trials by no. of time intervals)
beta1 <- array(rep(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),J*ints), dim=c(J, ints))
beta2 <- array(rep(c(0.8, 0.9, 1, 0.1, 0.2, 0.3, 0.4, 0.5),J*ints), dim=c(J, ints))
beta3 <- array(rep(c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 0.1, 0.2),J*ints), dim=c(J, ints))

phi1 <- c(rep(0.1, ntrt-1))
phi2 <- c(rep(0.2, ntrt-1))
phi3 <- c(rep(0.3, ntrt-1))

phic1 <- c(rep(-0.1, ntrt-1))
phic2 <- c(rep(0.1, ntrt-1))
phic3 <- c(rep(0.2, ntrt-1))

inits <- list(list(alpha=alpha1, beta=beta1, phi=phi1, phic=phic1), 
              list(alpha=alpha2, beta=beta2, phi=phi2, phic=phic2),
              list(alpha=alpha3, beta=beta3, phi=phi3, phic=phic3))

bugs.object<-bugs(data=nma_data, inits=inits, 
                  parameters.to.save=c("alpha", "phi", "phic", "S"), 
                  model.file="FE_model.txt", clearWD=F, 
                  summary.only=FALSE, n.iter=(num.sims+burn.in), 
                  n.sims=num.sims, n.burnin=burn.in, n.chains=3, 
                  bugs.seed=318659, bugs.directory=bugs.directory, 
                  debug=F, working.directory=path)


#-----------------------------------------------------------------------------
# Model Results
#-----------------------------------------------------------------------------

results <- bugs.object$summary

# Save results in csv file
write.csv(results,file="results.csv")

# Check results
results2 <- bugs.object$sims.matrix[,grep("alpha",rownames(bugs.object$summary))]
results2 <- cbind(rep(0,dim(results2)[1]),results2)
summary(results2)

results_mcmc<-mcmc(results2)
par(mfrow=c(3,2))

#Check autocorrelation
autocorr.plot(results_mcmc[,2:13])

# Check trace for convergence
traceplot(results_mcmc[,2])
traceplot(results_mcmc[,3])
traceplot(results_mcmc[,4])
traceplot(results_mcmc[,5])
traceplot(results_mcmc[,6])
traceplot(results_mcmc[,7])
traceplot(results_mcmc[,8])
traceplot(results_mcmc[,9])
traceplot(results_mcmc[,10])
traceplot(results_mcmc[,11])
traceplot(results_mcmc[,12])
traceplot(results_mcmc[,13])

# Histograms of posterior distributions
densplot(results_mcmc[,2])
densplot(results_mcmc[,3])
densplot(results_mcmc[,4])
densplot(results_mcmc[,5])
densplot(results_mcmc[,6])
densplot(results_mcmc[,7])
densplot(results_mcmc[,8])
densplot(results_mcmc[,9])
densplot(results_mcmc[,10])
densplot(results_mcmc[,11])
densplot(results_mcmc[,12])
densplot(results_mcmc[,13])

