# Calculate area under the survival curve for each iteration from WinBUGS

# Load libraries
library(foreach)
library(pracma)
library(Rmisc)
library(survival)

# Run HR NMA model to get BUGS results which are stored as an object bugs.fe
df2 <- bugs.fe$sims.list$hrd

# Calculate Kaplan-Meier for DTIC arm from CheckMate 066
ipd_data <- read.csv("../Data/melanoma_ipd.csv")
ipd_data <- ipd_data[ipd_data$studyCode==3 & ipd_data$txCode==1,]
KM.est <- survfit(Surv(time,event)~1, data=ipd_data, type="kaplan-meier", conf.int=FALSE)

# Break df2 into a series of vectors
d2 <- df2[,1]
d3 <- df2[,2]
d4 <- df2[,3]
d5 <- df2[,4]
d6 <- df2[,5]
d7 <- df2[,6]
d8 <- df2[,7]
d9 <- df2[,8]
d10 <- df2[,9]
d11 <- df2[,10]
d12 <- df2[,11]
d13 <- df2[,12]

# No. of iterations
n <- length(d2)

# Calculate cumulative hazard
cumhaz <- replicate(n, data.frame(S=NA))
names(cumhaz) <- paste0("df", 1:n)

foreach(i=1:n) %do% {
  cumhaz[[i]]$time <- c(0, KM.est$time)
  cumhaz[[i]]$trt1 <- c(0, KM.est$cumhaz)
  cumhaz[[i]]$trt2 <- cumhaz[[i]]$trt1*d2[i]
  cumhaz[[i]]$trt3 <- cumhaz[[i]]$trt1*d3[i]
  cumhaz[[i]]$trt4 <- cumhaz[[i]]$trt1*d4[i]
  cumhaz[[i]]$trt5 <- cumhaz[[i]]$trt1*d5[i]
  cumhaz[[i]]$trt6 <- cumhaz[[i]]$trt1*d6[i]
  cumhaz[[i]]$trt7 <- cumhaz[[i]]$trt1*d7[i]
  cumhaz[[i]]$trt8 <- cumhaz[[i]]$trt1*d8[i]
  cumhaz[[i]]$trt9 <- cumhaz[[i]]$trt1*d9[i]
  cumhaz[[i]]$trt10 <- cumhaz[[i]]$trt1*d10[i]
  cumhaz[[i]]$trt11 <- cumhaz[[i]]$trt1*d11[i]
  cumhaz[[i]]$trt12 <- cumhaz[[i]]$trt1*d12[i]
  cumhaz[[i]]$trt13 <- cumhaz[[i]]$trt1*d13[i]
}

# Calculate survival
surv <- replicate(n, data.frame(S=NA))
names(surv) <- paste0("df", 1:n)

foreach(i=1:n) %do% {
  surv[[i]]$time <- c(0, KM.est$time)
  surv[[i]]$trt1 <- c(1, KM.est$surv)
  surv[[i]]$trt2 <- exp(-cumhaz[[i]]$trt2)
  surv[[i]]$trt3 <- exp(-cumhaz[[i]]$trt3)
  surv[[i]]$trt4 <- exp(-cumhaz[[i]]$trt4)
  surv[[i]]$trt5 <- exp(-cumhaz[[i]]$trt5)
  surv[[i]]$trt6 <- exp(-cumhaz[[i]]$trt6)
  surv[[i]]$trt7 <- exp(-cumhaz[[i]]$trt7)
  surv[[i]]$trt8 <- exp(-cumhaz[[i]]$trt8)
  surv[[i]]$trt9 <- exp(-cumhaz[[i]]$trt9)
  surv[[i]]$trt10 <- exp(-cumhaz[[i]]$trt10)
  surv[[i]]$trt11 <- exp(-cumhaz[[i]]$trt11)
  surv[[i]]$trt12 <- exp(-cumhaz[[i]]$trt12)
  surv[[i]]$trt13 <- exp(-cumhaz[[i]]$trt13)
}


# Now calculate area under the curve for each treatment for each iteration

# Calculate AUC at max follow-up 

auc.trt1 <- c(rep(NA, n))
auc.trt2 <- c(rep(NA, n))
auc.trt3 <- c(rep(NA, n))
auc.trt4 <- c(rep(NA, n))
auc.trt5 <- c(rep(NA, n))
auc.trt6 <- c(rep(NA, n))
auc.trt7 <- c(rep(NA, n))
auc.trt8 <- c(rep(NA, n))
auc.trt9 <- c(rep(NA, n))
auc.trt10 <- c(rep(NA, n))
auc.trt11 <- c(rep(NA, n))
auc.trt12 <- c(rep(NA, n))
auc.trt13 <- c(rep(NA, n))

foreach(i=1:n) %do% {
  auc.trt1[i] <- trapz(surv[[i]]$time, surv[[i]]$trt1)
  auc.trt2[i] <- trapz(surv[[i]]$time, surv[[i]]$trt2)
  auc.trt3[i] <- trapz(surv[[i]]$time, surv[[i]]$trt3)
  auc.trt4[i] <- trapz(surv[[i]]$time, surv[[i]]$trt4)
  auc.trt5[i] <- trapz(surv[[i]]$time, surv[[i]]$trt5)
  auc.trt6[i] <- trapz(surv[[i]]$time, surv[[i]]$trt6)
  auc.trt7[i] <- trapz(surv[[i]]$time, surv[[i]]$trt7)
  auc.trt8[i] <- trapz(surv[[i]]$time, surv[[i]]$trt8)
  auc.trt9[i] <- trapz(surv[[i]]$time, surv[[i]]$trt9)
  auc.trt10[i] <- trapz(surv[[i]]$time, surv[[i]]$trt10)
  auc.trt11[i] <- trapz(surv[[i]]$time, surv[[i]]$trt11)
  auc.trt12[i] <- trapz(surv[[i]]$time, surv[[i]]$trt12)
  auc.trt13[i] <- trapz(surv[[i]]$time, surv[[i]]$trt13)
}


# Calculate mean AUC and 95% confidence interval
a1 <- CI(auc.trt1, ci=0.95)
a2 <- CI(auc.trt2, ci=0.95)
a3 <- CI(auc.trt3, ci=0.95)
a4 <- CI(auc.trt4, ci=0.95)
a5 <- CI(auc.trt5, ci=0.95)
a6 <- CI(auc.trt6, ci=0.95)
a7 <- CI(auc.trt7, ci=0.95)
a8 <- CI(auc.trt8, ci=0.95)
a9 <- CI(auc.trt9, ci=0.95)
a10 <- CI(auc.trt10, ci=0.95)
a11 <- CI(auc.trt11, ci=0.95)
a12 <- CI(auc.trt12, ci=0.95)
a13 <- CI(auc.trt13, ci=0.95)

auc60 <- data.frame(mean=c(a1[2], a2[2], a3[2], a4[2], a5[2], a6[2], a7[2], a8[2], a9[2], a10[2], 
                           a11[2], a12[2], a13[2]),
                    lower=c(a1[3], a2[3], a3[3], a4[3], a5[3], a6[3], a7[3], a8[3], a9[3], a10[3],
                            a11[3], a12[3], a13[3]),
                    upper=c(a1[1], a2[1], a3[1], a4[1], a5[1], a6[1], a7[1], a8[1], a9[1], a10[1],
                            a11[1], a12[1], a13[1]))

# Save as excel file
write.csv(auc60, "60month_AUC.csv")
