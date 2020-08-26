# Area under the survival curve for each iteration from WinBUGS

# Load libraries
library(foreach)
library(pracma)
library(Rmisc)

# Run piecewise model to get BUGS results which are stored as an object bugs.object
# Extract estimates of alpha
df2 <- bugs.object$sims.list$alpha

# beta for CheckMate 066
beta <- bugs.object$sims.list$beta[,3,] 

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

# Calculate log hazard for each iteration
loghaz <- replicate(n, data.frame(S=NA))
names(loghaz) <- paste0("df", 1:n)

foreach(i=1:n) %do% {
  loghaz[[i]]$time <- c(6, 12, 18, 24, 120)
  
  loghaz[[i]]$trt1 <- c(beta[i,1], beta[i,2], beta[i,3], beta[i,4], beta[i,5])
  
  loghaz[[i]]$trt2 <- c(loghaz[[i]]$trt1[1] + d2[i], loghaz[[i]]$trt1[2] + d2[i], loghaz[[i]]$trt1[3] + d2[i],
                        loghaz[[i]]$trt1[4] + d2[i], loghaz[[i]]$trt1[5] + d2[i])
  
  loghaz[[i]]$trt3 <- c(loghaz[[i]]$trt1[1] + d3[i], loghaz[[i]]$trt1[2] + d3[i], loghaz[[i]]$trt1[3] + d3[i],
                        loghaz[[i]]$trt1[4] + d3[i], loghaz[[i]]$trt1[5] + d3[i])
  
  loghaz[[i]]$trt4 <- c(loghaz[[i]]$trt1[1] + d4[i], loghaz[[i]]$trt1[2] + d4[i], loghaz[[i]]$trt1[3] + d4[i],
                        loghaz[[i]]$trt1[4] + d4[i], loghaz[[i]]$trt1[5] + d4[i])                     
  
  loghaz[[i]]$trt5 <- c(loghaz[[i]]$trt1[1] + d4[i] + d5[i], loghaz[[i]]$trt1[2] + d4[i] + d5[i], 
                        loghaz[[i]]$trt1[3] + d4[i] + d5[i], loghaz[[i]]$trt1[4] + d4[i] + d5[i], 
                        loghaz[[i]]$trt1[5] + d4[i] + d5[i]) 
  
  loghaz[[i]]$trt6 <- c(loghaz[[i]]$trt1[1] + d6[i], loghaz[[i]]$trt1[2] + d6[i], loghaz[[i]]$trt1[3] + d6[i],
                        loghaz[[i]]$trt1[4] + d6[i], loghaz[[i]]$trt1[5] + d6[i])  
  
  loghaz[[i]]$trt7 <- c(loghaz[[i]]$trt1[1] + d6[i] + d7[i], loghaz[[i]]$trt1[2] + d6[i] + d7[i], 
                        loghaz[[i]]$trt1[3] + d6[i] + d7[i], loghaz[[i]]$trt1[4] + d6[i] + d7[i], 
                        loghaz[[i]]$trt1[5] + d6[i] + d7[i])
  
  loghaz[[i]]$trt8 <- c(loghaz[[i]]$trt1[1] + d9[i] - d8[i], loghaz[[i]]$trt1[2] + d9[i] - d8[i], 
                        loghaz[[i]]$trt1[3] + d9[i] - d8[i], loghaz[[i]]$trt1[4] + d9[i] - d8[i], 
                        loghaz[[i]]$trt1[5] + d9[i] - d8[i])
  
  loghaz[[i]]$trt9 <- c(loghaz[[i]]$trt1[1] + d9[i], loghaz[[i]]$trt1[2] + d9[i], loghaz[[i]]$trt1[3] + d9[i],
                        loghaz[[i]]$trt1[4] + d9[i], loghaz[[i]]$trt1[5] + d9[i])  
  
  loghaz[[i]]$trt10 <- c(loghaz[[i]]$trt1[1] + d9[i] - d8[i] + d10[i], loghaz[[i]]$trt1[2] + d9[i] - d8[i] + d10[i], 
                         loghaz[[i]]$trt1[3] + d9[i] - d8[i] + d10[i], loghaz[[i]]$trt1[4] + d9[i] - d8[i] + d10[i], 
                         loghaz[[i]]$trt1[5] + d9[i] - d8[i] + d10[i])
  
  loghaz[[i]]$trt11 <- c(loghaz[[i]]$trt1[1] + d9[i] - d8[i] + d11[i], loghaz[[i]]$trt1[2] + d9[i] - d8[i] + d11[i], 
                         loghaz[[i]]$trt1[3] + d9[i] - d8[i] + d11[i], loghaz[[i]]$trt1[4] + d9[i] - d8[i] + d11[i], 
                         loghaz[[i]]$trt1[5] + d9[i] - d8[i] + d11[i])
  
  loghaz[[i]]$trt12 <- c(loghaz[[i]]$trt1[1] + d9[i] - d8[i] + d12[i], loghaz[[i]]$trt1[2] + d9[i] - d8[i] + d12[i], 
                         loghaz[[i]]$trt1[3] + d9[i] - d8[i] + d12[i], loghaz[[i]]$trt1[4] + d9[i] - d8[i] + d12[i], 
                         loghaz[[i]]$trt1[5] + d9[i] - d8[i] + d12[i])
  
  loghaz[[i]]$trt13 <- c(loghaz[[i]]$trt1[1] + d13[i], loghaz[[i]]$trt1[2] + d13[i], loghaz[[i]]$trt1[3] + d13[i],
                         loghaz[[i]]$trt1[4] + d13[i], loghaz[[i]]$trt1[5] + d13[i])
}

# Calculate hazard
haz <- replicate(n, data.frame(S=NA))
names(haz) <- paste0("df", 1:n)

foreach(i=1:n) %do% {
  haz[[i]]$time <- c(6, 12, 18, 24, 120)
  haz[[i]]$trt1 <- exp(loghaz[[i]]$trt1)
  haz[[i]]$trt2 <- exp(loghaz[[i]]$trt2)
  haz[[i]]$trt3 <- exp(loghaz[[i]]$trt3)
  haz[[i]]$trt4 <- exp(loghaz[[i]]$trt4)
  haz[[i]]$trt5 <- exp(loghaz[[i]]$trt5)
  haz[[i]]$trt6 <- exp(loghaz[[i]]$trt6)
  haz[[i]]$trt7 <- exp(loghaz[[i]]$trt7)
  haz[[i]]$trt8 <- exp(loghaz[[i]]$trt8)
  haz[[i]]$trt9 <- exp(loghaz[[i]]$trt9)
  haz[[i]]$trt10 <- exp(loghaz[[i]]$trt10)
  haz[[i]]$trt11 <- exp(loghaz[[i]]$trt11)
  haz[[i]]$trt12 <- exp(loghaz[[i]]$trt12)
  haz[[i]]$trt13 <- exp(loghaz[[i]]$trt13)
}

# Calculate cumulative hazard up to 60 months
cumhaz <- replicate(n, data.frame(S=NA))
names(cumhaz) <- paste0("df", 1:n)

foreach(i=1:n) %do% {
  cumhaz[[i]]$time <- c(6, 12, 18, 24, 60)
  
  cumhaz[[i]]$trt1 <- c(6*haz[[i]]$trt1[1], 
                        6*(haz[[i]]$trt1[1] + haz[[i]]$trt1[2]), 
                        6*(haz[[i]]$trt1[1] + haz[[i]]$trt1[2] + haz[[i]]$trt1[3]),
                        6*(haz[[i]]$trt1[1] + haz[[i]]$trt1[2] + haz[[i]]$trt1[3] + haz[[i]]$trt1[4]),
                        6*(haz[[i]]$trt1[1] + haz[[i]]$trt1[2] + haz[[i]]$trt1[3] + haz[[i]]$trt1[4]) + 36*haz[[i]]$trt1[5])
                   
  
  cumhaz[[i]]$trt2 <- c(6*haz[[i]]$trt2[1], 
                        6*(haz[[i]]$trt2[1] + haz[[i]]$trt2[2]), 
                        6*(haz[[i]]$trt2[1] + haz[[i]]$trt2[2] + haz[[i]]$trt2[3]),
                        6*(haz[[i]]$trt2[1] + haz[[i]]$trt2[2] + haz[[i]]$trt2[3] + haz[[i]]$trt2[4]),
                        6*(haz[[i]]$trt2[1] + haz[[i]]$trt2[2] + haz[[i]]$trt2[3] + haz[[i]]$trt2[4]) + 36*haz[[i]]$trt2[5])
                        
  cumhaz[[i]]$trt3 <- c(6*haz[[i]]$trt3[1], 
                        6*(haz[[i]]$trt3[1] + haz[[i]]$trt3[2]), 
                        6*(haz[[i]]$trt3[1] + haz[[i]]$trt3[2] + haz[[i]]$trt3[3]),
                        6*(haz[[i]]$trt3[1] + haz[[i]]$trt3[2] + haz[[i]]$trt3[3] + haz[[i]]$trt3[4]),
                        6*(haz[[i]]$trt3[1] + haz[[i]]$trt3[2] + haz[[i]]$trt3[3] + haz[[i]]$trt3[4]) + 36*haz[[i]]$trt3[5])
  
  cumhaz[[i]]$trt4 <- c(6*haz[[i]]$trt4[1], 
                        6*(haz[[i]]$trt4[1] + haz[[i]]$trt4[2]), 
                        6*(haz[[i]]$trt4[1] + haz[[i]]$trt4[2] + haz[[i]]$trt4[3]),
                        6*(haz[[i]]$trt4[1] + haz[[i]]$trt4[2] + haz[[i]]$trt4[3] + haz[[i]]$trt4[4]),
                        6*(haz[[i]]$trt4[1] + haz[[i]]$trt4[2] + haz[[i]]$trt4[3] + haz[[i]]$trt4[4]) + 36*haz[[i]]$trt4[5])
                        
  cumhaz[[i]]$trt5 <- c(6*haz[[i]]$trt5[1], 
                        6*(haz[[i]]$trt5[1] + haz[[i]]$trt5[2]), 
                        6*(haz[[i]]$trt5[1] + haz[[i]]$trt5[2] + haz[[i]]$trt5[3]),
                        6*(haz[[i]]$trt5[1] + haz[[i]]$trt5[2] + haz[[i]]$trt5[3] + haz[[i]]$trt5[4]),
                        6*(haz[[i]]$trt5[1] + haz[[i]]$trt5[2] + haz[[i]]$trt5[3] + haz[[i]]$trt5[4]) + 36*haz[[i]]$trt5[5])
                        
  cumhaz[[i]]$trt6 <- c(6*haz[[i]]$trt6[1], 
                        6*(haz[[i]]$trt6[1] + haz[[i]]$trt6[2]), 
                        6*(haz[[i]]$trt6[1] + haz[[i]]$trt6[2] + haz[[i]]$trt6[3]),
                        6*(haz[[i]]$trt6[1] + haz[[i]]$trt6[2] + haz[[i]]$trt6[3] + haz[[i]]$trt6[4]),
                        6*(haz[[i]]$trt6[1] + haz[[i]]$trt6[2] + haz[[i]]$trt6[3] + haz[[i]]$trt6[4]) + 36*haz[[i]]$trt6[5])
                       
  cumhaz[[i]]$trt7 <- c(6*haz[[i]]$trt7[1], 
                        6*(haz[[i]]$trt7[1] + haz[[i]]$trt7[2]), 
                        6*(haz[[i]]$trt7[1] + haz[[i]]$trt7[2] + haz[[i]]$trt7[3]),
                        6*(haz[[i]]$trt7[1] + haz[[i]]$trt7[2] + haz[[i]]$trt7[3] + haz[[i]]$trt7[4]),
                        6*(haz[[i]]$trt7[1] + haz[[i]]$trt7[2] + haz[[i]]$trt7[3] + haz[[i]]$trt7[4]) + 36*haz[[i]]$trt7[5])
                        
  cumhaz[[i]]$trt8 <- c(6*haz[[i]]$trt8[1], 
                        6*(haz[[i]]$trt8[1] + haz[[i]]$trt8[2]), 
                        6*(haz[[i]]$trt8[1] + haz[[i]]$trt8[2] + haz[[i]]$trt8[3]),
                        6*(haz[[i]]$trt8[1] + haz[[i]]$trt8[2] + haz[[i]]$trt8[3] + haz[[i]]$trt8[4]),
                        6*(haz[[i]]$trt8[1] + haz[[i]]$trt8[2] + haz[[i]]$trt8[3] + haz[[i]]$trt8[4]) + 36*haz[[i]]$trt8[5])
                        
  cumhaz[[i]]$trt9 <- c(6*haz[[i]]$trt9[1], 
                        6*(haz[[i]]$trt9[1] + haz[[i]]$trt9[2]), 
                        6*(haz[[i]]$trt9[1] + haz[[i]]$trt9[2] + haz[[i]]$trt9[3]),
                        6*(haz[[i]]$trt9[1] + haz[[i]]$trt9[2] + haz[[i]]$trt9[3] + haz[[i]]$trt9[4]),
                        6*(haz[[i]]$trt9[1] + haz[[i]]$trt9[2] + haz[[i]]$trt9[3] + haz[[i]]$trt9[4]) + 36*haz[[i]]$trt9[5])
                        
  cumhaz[[i]]$trt10 <- c(6*haz[[i]]$trt10[1], 
                         6*(haz[[i]]$trt10[1] + haz[[i]]$trt10[2]), 
                         6*(haz[[i]]$trt10[1] + haz[[i]]$trt10[2] + haz[[i]]$trt10[3]),
                         6*(haz[[i]]$trt10[1] + haz[[i]]$trt10[2] + haz[[i]]$trt10[3] + haz[[i]]$trt10[4]),
                         6*(haz[[i]]$trt10[1] + haz[[i]]$trt10[2] + haz[[i]]$trt10[3] + haz[[i]]$trt10[4]) + 36*haz[[i]]$trt10[5])
                         
  cumhaz[[i]]$trt11 <- c(6*haz[[i]]$trt11[1], 
                         6*(haz[[i]]$trt11[1] + haz[[i]]$trt11[2]), 
                         6*(haz[[i]]$trt11[1] + haz[[i]]$trt11[2] + haz[[i]]$trt11[3]),
                         6*(haz[[i]]$trt11[1] + haz[[i]]$trt11[2] + haz[[i]]$trt11[3] + haz[[i]]$trt11[4]),
                         6*(haz[[i]]$trt11[1] + haz[[i]]$trt11[2] + haz[[i]]$trt11[3] + haz[[i]]$trt11[4]) + 36*haz[[i]]$trt11[5])
                        
  cumhaz[[i]]$trt12 <- c(6*haz[[i]]$trt12[1], 
                         6*(haz[[i]]$trt12[1] + haz[[i]]$trt12[2]), 
                         6*(haz[[i]]$trt12[1] + haz[[i]]$trt12[2] + haz[[i]]$trt12[3]),
                         6*(haz[[i]]$trt12[1] + haz[[i]]$trt12[2] + haz[[i]]$trt12[3] + haz[[i]]$trt12[4]),
                         6*(haz[[i]]$trt12[1] + haz[[i]]$trt12[2] + haz[[i]]$trt12[3] + haz[[i]]$trt12[4]) + 36*haz[[i]]$trt12[5])
                         
  cumhaz[[i]]$trt13 <- c(6*haz[[i]]$trt13[1], 
                         6*(haz[[i]]$trt13[1] + haz[[i]]$trt13[2]), 
                         6*(haz[[i]]$trt13[1] + haz[[i]]$trt13[2] + haz[[i]]$trt13[3]),
                         6*(haz[[i]]$trt13[1] + haz[[i]]$trt13[2] + haz[[i]]$trt13[3] + haz[[i]]$trt13[4]),
                         6*(haz[[i]]$trt13[1] + haz[[i]]$trt13[2] + haz[[i]]$trt13[3] + haz[[i]]$trt13[4]) + 36*haz[[i]]$trt13[5])
                         
  
}

# Calculate survival
surv <- replicate(n, data.frame(S=NA))
names(surv) <- paste0("df", 1:n)

foreach(i=1:n) %do% {
  surv[[i]]$time <- c(0, 6, 12, 18, 24, 60)
  surv[[i]]$trt1 <- c(1, exp(-cumhaz[[i]]$trt1))
  surv[[i]]$trt2 <- c(1, exp(-cumhaz[[i]]$trt2))
  surv[[i]]$trt3 <- c(1, exp(-cumhaz[[i]]$trt3))
  surv[[i]]$trt4 <- c(1, exp(-cumhaz[[i]]$trt4))
  surv[[i]]$trt5 <- c(1, exp(-cumhaz[[i]]$trt5))
  surv[[i]]$trt6 <- c(1, exp(-cumhaz[[i]]$trt6))
  surv[[i]]$trt7 <- c(1, exp(-cumhaz[[i]]$trt7))
  surv[[i]]$trt8 <- c(1, exp(-cumhaz[[i]]$trt8))
  surv[[i]]$trt9 <- c(1, exp(-cumhaz[[i]]$trt9))
  surv[[i]]$trt10 <- c(1, exp(-cumhaz[[i]]$trt10))
  surv[[i]]$trt11 <- c(1, exp(-cumhaz[[i]]$trt11))
  surv[[i]]$trt12 <- c(1, exp(-cumhaz[[i]]$trt12))
  surv[[i]]$trt13 <- c(1, exp(-cumhaz[[i]]$trt13))
}


# Now calculate area under the curve for each treatment for each iteration

# Calculate AUC at 60 months 

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
