# Calculate area under the survival curve for each WinBUGS iteration

# Load libraries
library(foreach)
library(pracma)
library(Rmisc)

# Run GG model to get BUGS results which are stored as an object called nmaAggrFE

# Fit generalised gamma model just for the CheckMate 066 trial
df <- data[data$studyCode==3 & data$txCode==1,]

ma <- flexsurvreg(formula=Surv(time,event)~1,
                  data=df, dist="gengamma", method="BFGS")

# Identify coefficicents needed for predicting survival
mu <- ma$coefficients["mu"]
sigma <- exp(ma$coefficients["sigma"])
q <- ma$coefficients["Q"]


# Store estimates of beta from each iteration
df2 <- nmaAggrFE$sims.list$beta

# Break df2 into a series of vectors
beta2 <- df2[,1]
beta3 <- df2[,2]
beta4 <- df2[,3]
beta5 <- df2[,4]
beta6 <- df2[,5]
beta7 <- df2[,6]
beta8 <- df2[,7]
beta9 <- df2[,8]
beta10 <- df2[,9]
beta11 <- df2[,10]
beta12 <- df2[,11]
beta13 <- df2[,12]

# No. of iterations
n <- length(beta2)

# time values to calculate survival for
x <- seq(0,60,0.1)

# No. of time points needed for data frames
p <- length(x)

# Generate empty data frames ready for calculation of survival time
temp2 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp2) <- paste0("df", 1:n)

temp3 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp3) <- paste0("df", 1:n)

temp4 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp4) <- paste0("df", 1:n)

temp5 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp5) <- paste0("df", 1:n)

temp6 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp6) <- paste0("df", 1:n)

temp7 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp7) <- paste0("df", 1:n)

temp8 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp8) <- paste0("df", 1:n)

temp9 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp9) <- paste0("df", 1:n)

temp10 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp10) <- paste0("df", 1:n)

temp11 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp11) <- paste0("df", 1:n)

temp12 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp12) <- paste0("df", 1:n)

temp13 <- replicate(n, data.frame(S=c(rep(NA, p))))
names(temp13) <- paste0("df", 1:n)

# Fill the data frames
foreach(i=1:n) %do% {
  temp2[[i]] <- 1-pgengamma(x, mu = mu+beta2[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp3[[i]] <- 1-pgengamma(x, mu = mu+beta3[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp4[[i]] <- 1-pgengamma(x, mu = mu+beta4[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp5[[i]] <- 1-pgengamma(x, mu = mu+beta5[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp6[[i]] <- 1-pgengamma(x, mu = mu+beta6[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp7[[i]] <- 1-pgengamma(x, mu = mu+beta7[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp8[[i]] <- 1-pgengamma(x, mu = mu+beta8[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp9[[i]] <- 1-pgengamma(x, mu = mu+beta9[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp10[[i]] <- 1-pgengamma(x, mu = mu+beta10[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp11[[i]] <- 1-pgengamma(x, mu = mu+beta11[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp12[[i]] <- 1-pgengamma(x, mu = mu+beta12[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
  temp13[[i]] <- 1-pgengamma(x, mu = mu+beta13[i], sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
}

# Survival in DTIC arm
S.trt1 <- 1-pgengamma(x, mu = mu, sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)

# Now calculate area under the curve for each treatment for each iteration

tempg <- replicate(n, data.frame(S=c(rep(NA, p))))
names(tempg) <- paste0("df", 1:n)

foreach(i=1:n) %do% {
tempg[[i]] <- data.frame(time=x, trt1=S.trt1, trt2=temp2[[i]], trt3=temp3[[i]], trt4=temp4[[i]], trt5=temp5[[i]], 
                         trt6=temp6[[i]], trt7=temp7[[i]], trt8=temp8[[i]], trt9=temp9[[i]], trt10=temp10[[i]], 
                         trt11=temp11[[i]], trt12 = temp12[[i]], trt13=temp13[[i]])
}

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
  auc.trt1[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt1)
  auc.trt2[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt2)
  auc.trt3[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt3)
  auc.trt4[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt4)
  auc.trt5[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt5)
  auc.trt6[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt6)
  auc.trt7[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt7)
  auc.trt8[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt8)
  auc.trt9[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt9)
  auc.trt10[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt10)
  auc.trt11[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt11)
  auc.trt12[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt12)
  auc.trt13[i] <- trapz(tempg[[i]]$time, tempg[[i]]$trt13)
}

rm(temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13)

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
