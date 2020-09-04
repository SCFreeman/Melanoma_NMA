# Use Roche code to run fractional polynomial NMA for melanoma network

# 1st September 2020
# Suzanne Freeman

rm(list = ls()) 

library(survival)
library(doBy)

setwd("R:/SCFreeman/scf20/Survival/Fractional_Polynomial/FP3")

# Start with IPD data and convert to ANOVA format by choosing time points for aggregating over

data <- read.csv("../../Data/melanoma_ipd.csv")

# Add treatment as a factor variable
data$treatment[data$txCode==1] <- "DTIC"
data$treatment[data$txCode==2] <- "TRL"
data$treatment[data$txCode==3] <- "IPI+DTIC"
data$treatment[data$txCode==4] <- "DB"
data$treatment[data$txCode==5] <- "DB+TR"
data$treatment[data$txCode==6] <- "VM"
data$treatment[data$txCode==7] <- "VM+COB"
data$treatment[data$txCode==8] <- "IPI"
data$treatment[data$txCode==9] <- "NIV"
data$treatment[data$txCode==10] <- "NIV+IPI"
data$treatment[data$txCode==11] <- "PEM"
data$treatment[data$txCode==12] <- "IPI+SRG"
data$treatment[data$txCode==13] <- "SEL+DTIC"
data$treatment <- as.factor(data$treatment)

# Need to create treatment as a factor variable before applying this function
source("anova_data.R")

# Select time points for aggregating data
timepoints=c(6, 12, 120)

# Time points including zero
timepoints2=c(0, 6, 12, 120)


# Apply function
anova <- data.frame(spgrp=NA, treatment=NA, trialid=NA, y=NA, nevents=NA,
                    natrisk=NA, y.max=NA, start=NA, time=NA, trt=NA, treatnumf=NA,
                    studynumf=NA)

anova <- anova_data(timepoints=timepoints, timepoints2=timepoints2, ref.study=3, 
                      df=data)

# Save as csv file
write.csv(anova, "melanoma_aggregate_data.csv")

# Save data as km for use in Roche code below
km <- anova

#--------------------------------------------------------------------------------

####  PART 2: R-code to easily fit several FP models and to compare the AIC values
#### ------------

#list of models to be fitted - 1st order FP models
models <- list(
  "Exponential" = list(b1=function(x){0},b2=function(x){0}),
  "Weibull,p1=0" = list(b1=function(x){log(x)},b2=function(x){0}),
  "p1=0.5" = list(b1=function(x){x^0.5},b2=function(x){0}),
  "Gompertz,p1=1" = list(b1=function(x){x},b2=function(x){0}),
  "p1=2" = list(b1=function(x){x^2},b2=function(x){0}),
  "p1=3" = list(b1=function(x){x^3},b2=function(x){0}),
  "p1=-0.5" = list(b1=function(x){x^-0.5},b2=function(x){0}),
  "p1=-1" = list(b1=function(x){x^-1},b2=function(x){0}),
  "p1=-2" = list(b1=function(x){x^-2},b2=function(x){0}),
  "Second order, p1=3, p2=3" = list(b1=function(x){x^3},b2=function(x){x^3*log(x)}),
  "Second order, p1=3, p2=2" = list(b1=function(x){x^3},b2=function(x){x^2}),
  "Second order, p1=3, p2=1" = list(b1=function(x){x^3},b2=function(x){x}),
  "Second order, p1=3, p2=0.5" = list(b1=function(x){x^3},b2=function(x){x^-0.5}),
  "Second order, p1=3, p2=0" = list(b1=function(x){x^3},b2=function(x){log(x)}),
  "Second order, p1=3, p2=-0.5" = list(b1=function(x){x^3},b2=function(x){x^-0.5}),
  "Second order, p1=3, p2=-1" = list(b1=function(x){x^3},b2=function(x){x^-1}),
  "Second order, p1=3, p2=-2" = list(b1=function(x){x^3},b2=function(x){x^-2}),
  "Second order, p1=2, p2=2" = list(b1=function(x){x^2},b2=function(x){x^2*log(x)}),
  "Second order, p1=2, p2=1" = list(b1=function(x){x^2},b2=function(x){x}),
  "Second order, p1=2, p2=0.5" = list(b1=function(x){x^2},b2=function(x){x^0.5}),
  "Second order, p1=2, p2=0" = list(b1=function(x){x^2},b2=function(x){log(x)}),
  "Second order, p1=2, p2=-0.5" = list(b1=function(x){x^2},b2=function(x){x^-0.5}),
  "Second order, p1=2, p2=-1" = list(b1=function(x){x^2},b2=function(x){x^-1}),
  "Second order, p1=2, p2=-2" = list(b1=function(x){x^2},b2=function(x){x^-2}),
  "Second order, p1=1, p2=1" = list(b1=function(x){x},b2=function(x){x*log(x)}),
  "Second order, p1=1, p2=0.5" = list(b1=function(x){x},b2=function(x){x^0.5}),
  "Second order, p1=1, p2=0" = list(b1=function(x){x},b2=function(x){log(x)}),
  "Second order, p1=1, p2=-0.5" = list(b1=function(x){x},b2=function(x){x^-0.5}),
  "Second order, p1=1, p2=-1" = list(b1=function(x){x},b2=function(x){x^-1}),
  "Second order, p1=1, p2=-2" = list(b1=function(x){x},b2=function(x){x^-2}),
  "Second order, p1=0.5, p2=0.5" = list(b1=function(x){x^0.5},b2=function(x){x^0.5*log(x)}),
  "Second order, p1=0.5, p2=0" = list(b1=function(x){x^0.5},b2=function(x){log(x)}),
  "Second order, p1=0.5, p2=-0.5" = list(b1=function(x){x^0.5},b2=function(x){x^-0.5}),
  "Second order, p1=0.5, p2=-1" = list(b1=function(x){x^0.5},b2=function(x){x^-1}),
  "Second order, p1=0.5, p2=-2" = list(b1=function(x){x^0.5},b2=function(x){x^-2}),
  "Second order, p1=0, p2=0" = list(b1=function(x){log(x)},b2=function(x){log(x)*log(x)}),
  "Second order, p1=0, p2=-0.5" = list(b1=function(x){log(x)},b2=function(x){x^-0.5}),
  "Second order, p1=0, p2=-1" = list(b1=function(x){log(x)},b2=function(x){x^-1}),
  "Second order, p1=0, p2=-2" = list(b1=function(x){log(x)},b2=function(x){x^-2}),
  "Second order, p1=-0.5, p2=-0.5" = list(b1=function(x){x^-0.5},b2=function(x){x^-0.5*log(x)}),
  "Second order, p1=-0.5, p2=-1" = list(b1=function(x){x^-0.5},b2=function(x){x^-1}),
  "Second order, p1=-0.5, p2=-2" = list(b1=function(x){x^-0.5},b2=function(x){x^-2}),
  "Second order, p1=-1, p2=-1" = list(b1=function(x){x^-1},b2=function(x){x^-1*log(x)}),
  "Second order, p1=-1, p2=-2" = list(b1=function(x){x^-1},b2=function(x){x^-2}),
  "Second order, p1=-2, p2=-2" = list(b1=function(x){x^-2},b2=function(x){x^-2*log(x)})
)

#Fit all models
fit.KM.NMA<-function(bf){
  km.new=km
  km.new$beta1=bf[[1]](km.new$time)
  km.new$beta2=bf[[2]](km.new$time)
  #model formula
  f=cbind(nevents,natrisk-nevents)~treatnumf+studynumf+treatnumf*beta1+treatnumf*beta2+studynumf*beta1+studynumf*beta2
  glm(f,family=binomial(link=cloglog),data=km.new)
}
fits=lapply(models,fit.KM.NMA)
#Get AIC from each model
aics=lapply(fits,AIC)
#Print the AICs
data.frame(AIC=round(unlist(aics),2))

# Sort AIC into ascending order
a <- data.frame(AIC=round(unlist(aics),2))
sort(a[,1])
#--------------------------------------------------------------------------------

####  PART 3: Fit different fractional polynomial models
#### ------------

#--------------------------------------------------------------------------------
#First order model with p1=0.5 
km.new       <- km
km.new$beta1 <- km.new$time^0.5 #beta1 is the time variable as first fp parameter


#model formula
f <-  cbind(nevents,natrisk-nevents) ~
  treatnumf + studynumf +              #beta0 terms for treatment and study
  treatnumf*beta1 + studynumf*beta1   #beta1 terms for treatment and study

#fit the model
fit <- glm( formula = f,
            family  = binomial(link=cloglog),
            data    = km.new)
#Data processing to extract the contrast coefficient

coeff.data            <- as.data.frame(coef(summary(fit)))[,1:2]
names(coeff.data)     <- c("est","std.err")
attach(coeff.data)
coeff.data$conf.int.lower   <- est-std.err*qnorm(0.975)
coeff.data$conf.int.upper   <- est+std.err*qnorm(0.975)
coeff.data[,1:4]      <- round(coeff.data[,1:4],3)
coeff.data$pn         <- rownames(coeff.data)
coeff.data            <- coeff.data[grepl("treatnumf",coeff.data$pn),]
coeff.data$contrast   <- paste(as.vector(gsub("treatnumf|beta1|beta2|:","",coeff.data$pn)),"vs",ref.trt)
coeff.data$parameter  <- ifelse(grepl("beta1",coeff.data$pn),"D1",ifelse(grepl("beta2",coeff.data$pn),"D2","D0"))
detach(coeff.data)

#Print the parameter estimates and AIC
coeff.data[,c("parameter","contrast","est","std.err","conf.int.lower","conf.int.upper")]
AIC(fit)


#--------------------------------------------------------------------------------
#First order model with p1=0 
km.new       <- km
km.new$beta1 <- log(km.new$time) #beta1 is the time variable as first fp parameter


#model formula
f <-  cbind(nevents,natrisk-nevents) ~
  treatnumf + studynumf +              #beta0 terms for treatment and study
  treatnumf*beta1 + studynumf*beta1   #beta1 terms for treatment and study

#fit the model
fit <- glm( formula = f,
            family  = binomial(link=cloglog),
            data    = km.new)
#Data processing to extract the contrast coefficient

coeff.data            <- as.data.frame(coef(summary(fit)))[,1:2]
names(coeff.data)     <- c("est","std.err")
attach(coeff.data)
coeff.data$conf.int.lower   <- est-std.err*qnorm(0.975)
coeff.data$conf.int.upper   <- est+std.err*qnorm(0.975)
coeff.data[,1:4]      <- round(coeff.data[,1:4],3)
coeff.data$pn         <- rownames(coeff.data)
coeff.data            <- coeff.data[grepl("treatnumf",coeff.data$pn),]
coeff.data$contrast   <- paste(as.vector(gsub("treatnumf|beta1|beta2|:","",coeff.data$pn)),"vs",ref.trt)
coeff.data$parameter  <- ifelse(grepl("beta1",coeff.data$pn),"D1",ifelse(grepl("beta2",coeff.data$pn),"D2","D0"))
detach(coeff.data)

#Print the parameter estimates and AIC
coeff.data[,c("parameter","contrast","est","std.err","conf.int.lower","conf.int.upper")]
AIC(fit)


#--------------------------------------------------------------------------------
#Second order model with p1=-0.5 and p2=-0.5
km.new       <- km
km.new$beta1 <- km.new$time^-0.5 #beta1 is the time variable as first fp parameter
km.new$beta2 <- km.new$time^-0.5*log(km.new$time)    #beta2 is the time variable as second fp parameter

#model formula
f <- cbind(nevents,natrisk-nevents)~
  treatnumf + studynumf +               #beta0 terms for treatment and study
  treatnumf*beta1 + studynumf*beta1 +   #beta1 terms for treatment and study
  treatnumf*beta2 + studynumf*beta2     #beta2 terms for treatment and study

#fit the model
fit <- glm( formula = f,
            family  = binomial(link=cloglog),
            data    = km.new)

#Data processing to extract the contrast coefficient

coeff.data            <- as.data.frame(coef(summary(fit)))[,1:2]
names(coeff.data)     <- c("est","std.err")
attach(coeff.data)
coeff.data$conf.int.lower   <- est-std.err*qnorm(0.975)
coeff.data$conf.int.upper   <- est+std.err*qnorm(0.975)
coeff.data[,1:4]      <- round(coeff.data[,1:4],3)
coeff.data$pn         <- rownames(coeff.data)
coeff.data            <- coeff.data[grepl("treatnumf",coeff.data$pn),]
coeff.data$contrast   <- paste(as.vector(gsub("treatnumf|beta1|beta2|:","",coeff.data$pn)),"vs",ref.trt)
coeff.data$parameter  <- ifelse(grepl("beta1",coeff.data$pn),"D1",ifelse(grepl("beta2",coeff.data$pn),"D2","D0"))
detach(coeff.data)

#Print the parameter estimates and AIC
coeff.data[,c("parameter","contrast","est","std.err","conf.int.lower","conf.int.upper")]
AIC(fit)


#--------------------------------------------------------------------------------
#Second order model with p1=0 and p2=-1
km.new       <- km
km.new$beta1 <- log(km.new$time) #beta1 is the time variable as first fp parameter
km.new$beta2 <- km.new$time^-1    #beta2 is the time variable as second fp parameter

#model formula
f <- cbind(nevents,natrisk-nevents)~
  treatnumf + studynumf +               #beta0 terms for treatment and study
  treatnumf*beta1 + studynumf*beta1 +   #beta1 terms for treatment and study
  treatnumf*beta2 + studynumf*beta2     #beta2 terms for treatment and study

#fit the model
fit <- glm( formula = f,
            family  = binomial(link=cloglog),
            data    = km.new)

#Data processing to extract the contrast coefficient

coeff.data            <- as.data.frame(coef(summary(fit)))[,1:2]
names(coeff.data)     <- c("est","std.err")
attach(coeff.data)
coeff.data$conf.int.lower   <- est-std.err*qnorm(0.975)
coeff.data$conf.int.upper   <- est+std.err*qnorm(0.975)
coeff.data[,1:4]      <- round(coeff.data[,1:4],3)
coeff.data$pn         <- rownames(coeff.data)
coeff.data            <- coeff.data[grepl("treatnumf",coeff.data$pn),]
coeff.data$contrast   <- paste(as.vector(gsub("treatnumf|beta1|beta2|:","",coeff.data$pn)),"vs",ref.trt)
coeff.data$parameter  <- ifelse(grepl("beta1",coeff.data$pn),"D1",ifelse(grepl("beta2",coeff.data$pn),"D2","D0"))
detach(coeff.data)

#Print the parameter estimates and AIC
coeff.data[,c("parameter","contrast","est","std.err","conf.int.lower","conf.int.upper")]
AIC(fit)


#--------------------------------------------------------------------------------
#Second order model with p1=-0.5 and p2=-1
km.new       <- km
km.new$beta1 <- km.new$time^-0.5 #beta1 is the time variable as first fp parameter
km.new$beta2 <- km.new$time^-1    #beta2 is the time variable as second fp parameter

#model formula
f <- cbind(nevents,natrisk-nevents)~
  treatnumf + studynumf +               #beta0 terms for treatment and study
  treatnumf*beta1 + studynumf*beta1 +   #beta1 terms for treatment and study
  treatnumf*beta2 + studynumf*beta2     #beta2 terms for treatment and study

#fit the model
fit <- glm( formula = f,
            family  = binomial(link=cloglog),
            data    = km.new)

#Data processing to extract the contrast coefficient

coeff.data            <- as.data.frame(coef(summary(fit)))[,1:2]
names(coeff.data)     <- c("est","std.err")
attach(coeff.data)
coeff.data$conf.int.lower   <- est-std.err*qnorm(0.975)
coeff.data$conf.int.upper   <- est+std.err*qnorm(0.975)
coeff.data[,1:4]      <- round(coeff.data[,1:4],3)
coeff.data$pn         <- rownames(coeff.data)
coeff.data            <- coeff.data[grepl("treatnumf",coeff.data$pn),]
coeff.data$contrast   <- paste(as.vector(gsub("treatnumf|beta1|beta2|:","",coeff.data$pn)),"vs",ref.trt)
coeff.data$parameter  <- ifelse(grepl("beta1",coeff.data$pn),"D1",ifelse(grepl("beta2",coeff.data$pn),"D2","D0"))
detach(coeff.data)

#Print the parameter estimates and AIC
coeff.data[,c("parameter","contrast","est","std.err","conf.int.lower","conf.int.upper")]
AIC(fit)


#--------------------------------------------------------------------------------

sessionInfo()
date()

####  END OF PROGRAM
#### ------------