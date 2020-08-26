# KM plot of Dacarbazine arms

# Load libraries
library(survival)
library(survminer)
library(dplyr)
library(RColorBrewer)

# Start with empty environment
rm(list=ls())

# Set the working directory
path <- "R:/SCFreeman/scf20/Melanoma_NMA/Data/"
setwd(path)

# Import data
data <- read.csv("melanoma_ipd.csv")

# Reduce to Dacarbazine arms only
data <- data[data$txCode==1,]
df1 <- data[data$studyCode==1,]
df2 <- data[data$studyCode==2,]
df3 <- data[data$studyCode==3,]
df11 <- data[data$studyCode==11,]
df12 <- data[data$studyCode==12,]
df13 <- data[data$studyCode==13,]

fit1 <- survfit(Surv(time,event)~1, data=df1)
fit2 <- survfit(Surv(time,event)~1, data=df2)
fit3 <- survfit(Surv(time,event)~1, data=df3)
fit4 <- survfit(Surv(time,event)~1, data=df11)
fit5 <- survfit(Surv(time,event)~1, data=df12)
fit6 <- survfit(Surv(time,event)~1, data=df13)

# Taken from Set1 palette from RcolorBrewer
color=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#f781bf")

# List of survfit objects needed for ggsurvplot
fit <- list(fit1=fit1, fit2=fit2, fit3=fit3, fit4=fit4, fit5=fit5, fit6=fit6)

# use combine=T to get multiple lines on one graph
p <- ggsurvplot(fit=fit, data=data, combine=T, censor=F, palette=color,
                xlab="Time (months)", ylab="Overall Survival", 
                xlim=c(0, 72), ylim=c(0, 1),
                legend=c(0.7,0.8), legend.title="", 
                legend.labs=c("BREAK3", "BRIM3", "CheckMate 066", "Ribas 2013", "Robert 2011",
                              "Robert 2013"),
                font.legend=16)
p$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

# save plot
setwd("R:/SCFreeman/scf20/Melanoma_NMA/Exploratory analyses")
dev.copy(pdf, "DTIC_arms.pdf")
dev.off()
