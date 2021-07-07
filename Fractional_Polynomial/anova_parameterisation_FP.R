# Format data for use with fractional polynomial models

rm(list = ls()) 

library(survival)
library(doBy)

setwd("R:/SCFreeman/scf20/Melanoma_NMA/Fractional_Polynomial")

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
