# Fit NMA of Cox HR in WinBUGS

###############################################################

# Load library
library(R2WinBUGS)

# Empty environment
rm(list=ls())

# Set the working directory
path <- "R:/SCFreeman/scf20/Melanoma_NMA/HR/Fixed"
setwd(path)

# load HR data
data <- read.csv("hr_data.csv")

# Now create a version which has one row per study
# At the moment the three-arm trial has two rows of data and need to change this to one row
data_wide <- data.frame(study=data$study)
data_wide$t1 <- data$treat1
data_wide$t2 <- data$treat2
data_wide$t3 <- NA
data_wide$y2 <- data$lhr
data_wide$y3 <- NA
data_wide$se2 <- data$se
data_wide$se3 <- NA
data_wide$V <- data$V

data_wide$t3[13] <- 10
data_wide$y3[13] <- data$lhr[14]
data_wide$se3[13] <- data$se[14]
data_wide <- data_wide[1:13,]

data <- data_wide
data$y1 <- 0
data$se1 <- 0

# Set the location for WinBUGS
bugs.directory <- "Z:\\winbugs143_unrestricted\\winbugs14_full_patched\\WinBUGS14"

# WinBUGS burn-in & simulation size
num.sims <- 60000
burn.in <- 60000


#-----------------------------------------------------------------------------
# Data
#-----------------------------------------------------------------------------

# No. of studies
ns <- nrow(data)

# No. of treatments
nt <- max(data$t2)

# No. of arms in each trial
data$na <- 2
data$na[data$study=="CheckMate 067"] <- 3

y <- array(c(data$y1, data$y2, data$y3), dim=c(ns,3))
se <- array(c(data$se1, data$se2, data$se3), dim=c(ns,3))
t <- array(c(data$t1, data$t2, data$t3), dim=c(ns,3))

bugs_data <- list(ns2=12, ns3=1, nt=nt, t=t, y=y, se=se, na=data$na, V=data$V)

#-----------------------------------------------------------------------------
# Initial values
#-----------------------------------------------------------------------------

d1 <- c(NA, rep(0,12))
d2 <- c(NA, rep(0.1, 12))
d3 <- c(NA, rep(-0.1, 12))

fe_inits <- list(list(d=d1), 
                 list(d=d2),
                 list(d=d3))

#-----------------------------------------------------------------------------
# Fit FE model in WinBUGS
#-----------------------------------------------------------------------------

bugs.fe <- bugs(data=bugs_data, inits=fe_inits, 
                parameters.to.save=c("d", "hrd", "best", "prob"), 
                model.file="FE_model.txt", clearWD=F, 
                summary.only=FALSE, n.iter=(num.sims+burn.in), 
                n.sims=num.sims, n.burnin=burn.in, n.chains=3, 
                bugs.seed=385916, bugs.directory=bugs.directory, 
                debug=F, working.directory=path)

fe_results <- bugs.fe$summary

# Save results in csv file
write.csv(fe_results,file="fe_results.csv")

# Check results
results2 <- bugs.fe$sims.matrix[,grep("d",rownames(bugs.fe$summary))]
results2 <- cbind(rep(0,dim(results2)[1]),results2)
summary(results2)

results_mcmc<-mcmc(results2)
par(mfrow=c(3,2))

# Check autocorrelation
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

