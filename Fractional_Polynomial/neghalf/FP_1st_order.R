# Fractional Polynomial model - 1st Order fixed effect

##################################################################################

# Start with an empty environment
rm(list = ls())

# Load libraries
library(coda)
library(survival)
library(R2WinBUGS)

# Set the working directory
path <- "R:/SCFreeman/scf20/Melanoma_NMA/Fractional_Polynomial/neghalf"
setwd(path)

# Load the data
data <- read.csv("../melanoma_aggregate_data.csv")

# Create a treatment code variable as an integer
data$txCode[data$treatment=="DTIC"] <- 1
data$txCode[data$treatment=="TRL"] <- 2
data$txCode[data$treatment=="IPI+DTIC"] <- 3
data$txCode[data$treatment=="DB"] <- 4
data$txCode[data$treatment=="DB+TR"] <- 5
data$txCode[data$treatment=="VM"] <- 6
data$txCode[data$treatment=="VM+COB"] <- 7
data$txCode[data$treatment=="IPI"] <- 8
data$txCode[data$treatment=="NIV"] <- 9
data$txCode[data$treatment=="NIV+IPI"] <- 10
data$txCode[data$treatment=="PEM"] <- 11
data$txCode[data$treatment=="IPI+SRG"] <- 12
data$txCode[data$treatment=="SEL+DTIC"] <- 13

# Order data
data <- data[order(data$trialid, data$txCode, data$spgrp),]

# Set the location for WinBUGS
bugs.directory <- "Z:\\winbugs143_unrestricted\\winbugs14_full_patched\\WinBUGS14"

# WinBUGS burn-in & simulation size
num.sims <- 90000
burn.in <- 90000

# Fractional polynomial powers
P1 <- -0.5

#-----------------------------------------------------------------------------
# Data formatting
#-----------------------------------------------------------------------------

# Need to number the treatment arms within each trial
data$arm[data$treatment=="DTIC"] <- 1
data$arm[data$treatment=="TRL" | data$treatment=="IPI+DTIC" | data$treatment=="SEL+DTIC" | data$treatment=="IPI+SRG" 
         | data$treatment=="PEM" | data$treatment=="VM+COB"] <- 2
data$arm[data$trialid==1 & data$treatment=="DTIC"] <- 1
data$arm[data$trialid==1 & data$treatment=="DB"] <- 2
data$arm[data$trialid==7 & data$treatment=="DB"] <- 1
data$arm[data$trialid==7 & data$treatment=="DB+TR"] <- 2
data$arm[data$trialid==8 & data$treatment=="DB+TR"] <- 1
data$arm[data$trialid==8 & data$treatment=="VM"] <- 2
data$arm[data$trialid==6 & data$treatment=="VM"] <- 1
data$arm[data$trialid==3 & data$treatment=="NIV"] <- 2
data$arm[data$trialid==10 & data$treatment=="IPI"] <- 1
data$arm[data$trialid==9 & data$treatment=="IPI"] <- 1
data$arm[data$trialid==5 & data$treatment=="IPI"] <- 1
data$arm[data$trialid==5 & data$treatment=="NIV+IPI"] <- 2
data$arm[data$trialid==4 & data$treatment=="IPI"] <- 1
data$arm[data$trialid==4 & data$treatment=="NIV"] <- 2
data$arm[data$trialid==4 & data$treatment=="NIV+IPI"] <- 3
data$arm[data$trialid==2 & data$treatment=="VM"] <- 2

# Check all arms coded
data$arm==1 | data$arm==2 | data$arm==3

# Length of time intervals
data$length <- data$time-data$start

# Number of treatments
nt <- length(unique(data$treatment))

# Number of studies
ns <- length(unique(data$trialid)) 

# Number of rows in dataset
N <- nrow(data)

# Maximum time
maxt <- 120
  
# Mean & precision
mean <- c(0,0)
prec <- array(c(0.0001, 0, 0, 0.0001), dim=c(2,2))

# Number of treatment arms for each trial
na <- c(2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2)

# Treatment in each trial arm - This fills the down the columns first
t <- array(data=c(1, 1, 1, 8, 8, 6, 4, 5, 8, 8, 1, 1, 1,
                  4, 6, 9, 9, 10, 7, 5, 6, 12, 11, 2, 3, 13,
                  NA, NA, NA, 10, rep(NA, 9)), dim=c(13, 3))

bugs_data <- list(s=data$trialid, r=data$nevents, z=data$natrisk, a=data$arm, time=data$time,
                 dt=data$length, P1=P1, N=N, nt=nt, ns=ns, maxt=maxt, mean=mean, prec=prec,
                 t=t,  na=na)


#-----------------------------------------------------------------------------
# Model
#-----------------------------------------------------------------------------

# Create initial values for model
d1 <- array(c(NA, rep(0.1, 12), NA, rep(0.2, 12)), dim=c(nt,2))
d2 <- array(c(NA, rep(0.2, 12), NA, rep(-0.1, 12)), dim=c(nt,2))
d3 <- array(c(NA, rep(-0.1, 12), NA, rep(0.1, 12)), dim=c(nt,2))

mu1 <- array(rep(c(0.4, 0.5, 0.6, 0.2, 0.3, 0.1, 0.1),4), dim=c(ns,2))
mu2 <- array(rep(c(0.3, 0.4, 0.5, 0.6, 0.7, -0.1, -0.2),4), dim=c(ns,2))
mu3 <- array(rep(c(0.5, 0.6, 0.7, 0.1, -0.1,  0.2, 0.2),4), dim=c(ns,2))

inits <- list(list(d=d1, mu=mu1), 
              list(d=d2, mu=mu2),
              list(d=d3, mu=mu3))

bugs.object <- bugs(data=bugs_data, inits=inits, 
                  parameters.to.save=c("d", "S", "rk60"), 
                  model.file="FE_1st_order_model.txt", clearWD=F, 
                  summary.only=FALSE, n.iter=(num.sims+burn.in), 
                  n.sims=num.sims, n.burnin=burn.in, n.chains=3, 
                  bugs.seed=385916, bugs.directory=bugs.directory, 
                  debug=F, working.directory=path)

results <- bugs.object$summary

# Save results in csv file
write.csv(results,file="results.csv")


#-----------------------------------------------------------------------------
# Model Results
#-----------------------------------------------------------------------------

results2 <- bugs.object$sims.matrix[,grep("d",rownames(bugs.object$summary))]
results2 <- cbind(rep(0,dim(results2)[1]),results2)
summary(results2)

results_mcmc<-mcmc(results2)
par(mfrow=c(3,3))

#Check autocorrelation
autocorr.plot(results_mcmc[,2:24])

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
traceplot(results_mcmc[,14])
traceplot(results_mcmc[,15])
traceplot(results_mcmc[,16])
traceplot(results_mcmc[,17])
traceplot(results_mcmc[,18])
traceplot(results_mcmc[,19])
traceplot(results_mcmc[,20])
traceplot(results_mcmc[,21])
traceplot(results_mcmc[,22])
traceplot(results_mcmc[,23])
traceplot(results_mcmc[,24])

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
densplot(results_mcmc[,14])
densplot(results_mcmc[,15])
densplot(results_mcmc[,16])
densplot(results_mcmc[,17])
densplot(results_mcmc[,18])
densplot(results_mcmc[,19])
densplot(results_mcmc[,20])
densplot(results_mcmc[,21])
densplot(results_mcmc[,22])
densplot(results_mcmc[,23])
densplot(results_mcmc[,24])

