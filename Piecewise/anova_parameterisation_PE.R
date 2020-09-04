# Aggregate IPD for fitting piecewise models

# Load libraries
library(survival)
library(doBy)

# Start with empty environment
rm(list = ls()) 

# Set the workign directory
setwd("R:/SCFreeman/scf20/Melanoma_NMA/Piecewise")

# Start with IPD data
data <- read.csv("../Data/melanoma_ipd.csv")

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

# Import function for aggregating data
source("anova_data.R")

# Select time points for aggregating data
timepoints=c(6, 12, 120)

# Time points including zero
timepoints2=c(0, 6, 12, 120)

# Empty data frame for aggregated data
anova <- data.frame(spgrp=NA, treatment=NA, trialid=NA, y=NA, nevents=NA,
                    natrisk=NA, y.max=NA, start=NA, time=NA, trt=NA, treatnumf=NA,
                    studynumf=NA)

# Apply function
anova <- anova_data(timepoints=timepoints, timepoints2=timepoints2, ref.study=3, 
                    df=data)

# Save as csv file
write.csv(anova, "melanoma_aggregate_data.csv")
