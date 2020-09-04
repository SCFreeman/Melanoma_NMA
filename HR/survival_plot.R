# Plot survival curves for each treatment

# Load library
library(survival)

# Start with empty environment
rm(list=ls())

# Load results
setwd("R:/SCFreeman/scf20/Melanoma_NMA/HR")
data <- read.csv("fe_results.csv")
data <- data[1:24,]

# Calculate Kaplan-Meier for DTIC arm from CheckMate 066
ipd_data <- read.csv("../Data/melanoma_ipd.csv")
ipd_data <- ipd_data[ipd_data$studyCode==3 & ipd_data$txCode==1,]
KM.est <- survfit(Surv(time,event)~1, data=ipd_data, type="kaplan-meier", conf.int=FALSE)


# Calculate cumulative hazard
cumhaz <- data.frame(time=c(0, KM.est$time), trt1=c(0, KM.est$cumhaz))
cumhaz$trt2 <- cumhaz$trt1*data$mean[13]
cumhaz$trt3 <- cumhaz$trt1*data$mean[14]
cumhaz$trt4 <- cumhaz$trt1*data$mean[15]
cumhaz$trt5 <- cumhaz$trt1*data$mean[16]
cumhaz$trt6 <- cumhaz$trt1*data$mean[17]
cumhaz$trt7 <- cumhaz$trt1*data$mean[18]
cumhaz$trt8 <- cumhaz$trt1*data$mean[19]
cumhaz$trt9 <- cumhaz$trt1*data$mean[20]
cumhaz$trt10 <- cumhaz$trt1*data$mean[21]
cumhaz$trt11 <- cumhaz$trt1*data$mean[22]
cumhaz$trt12 <- cumhaz$trt1*data$mean[23]
cumhaz$trt13 <- cumhaz$trt1*data$mean[24]

# Calcuate survival
surv <- data.frame(time=cumhaz$time, trt0=c(1, KM.est$surv))
surv$trt1 <- exp(-cumhaz$trt1)
surv$trt2 <- exp(-cumhaz$trt2)
surv$trt3 <- exp(-cumhaz$trt3)
surv$trt4 <- exp(-cumhaz$trt4)
surv$trt5 <- exp(-cumhaz$trt5)
surv$trt6 <- exp(-cumhaz$trt6)
surv$trt7 <- exp(-cumhaz$trt7)
surv$trt8 <- exp(-cumhaz$trt8)
surv$trt9 <- exp(-cumhaz$trt9)
surv$trt10 <- exp(-cumhaz$trt10)
surv$trt11 <- exp(-cumhaz$trt11)
surv$trt12 <- exp(-cumhaz$trt12)
surv$trt13 <- exp(-cumhaz$trt13)

graph_data <- surv

# Colors for the plot
color=c("#070d43", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
        "#cab2d6", "#6a3d9a", "#581845", "#b15928", "#184558")

par(xpd = T, mar = par()$mar + c(6,0,0,0))

# Start by plotting the Kaplan_Meier DTIC curve
plot(KM.est,xlab="Time (months)",ylab="Overall Survival",xaxt="n",yaxt="n",main=" ",xlim=c(0,60),ylim=c(0,1),
     mark.time=FALSE, col=color[14], conf.int=F)
#Add y axis (2 specifies that axis goes on the left of the plot)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
#Add x axis (1 specified that axis goes at the bottom of the plot)
axis(1, at=c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60))

# Add prediction lines
lines(graph_data$time, graph_data$trt1, col=color[1])
lines(graph_data$time, graph_data$trt2, col=color[2])
lines(graph_data$time, graph_data$trt3, col=color[3])
lines(graph_data$time, graph_data$trt4, col=color[4])
lines(graph_data$time, graph_data$trt5, col=color[5])
lines(graph_data$time, graph_data$trt6, col=color[6])
lines(graph_data$time, graph_data$trt7, col=color[7])
lines(graph_data$time, graph_data$trt8, col=color[8])
lines(graph_data$time, graph_data$trt9, col=color[9])
lines(graph_data$time, graph_data$trt10, col=color[10])
lines(graph_data$time, graph_data$trt11, col=color[11])
lines(graph_data$time, graph_data$trt12, col=color[12])
lines(graph_data$time, graph_data$trt13, col=color[13])

legend(0, -0.3,
       c("KM", "DTIC", "TRL", "IPI+DTIC", "DB", "DB+TR", "VM", "VM+COB", "IPI", "NIV", "NIV+IPI",
         "PEM", "IPI+SRG", "SEL+DTIC"),
       col=c(color[14], color[1], color[2], color[3], color[4], color[5], color[6], color[7], color[8], color[9],
             color[10], color[11], color[12], color[13]),
       lty=c(1,1,1,1), ncol=3, text.width=12, box.lty=0)

# save plot
dev.copy(pdf, "survival_plot.pdf")
dev.off()

library(pracma)

# Calculate AUC at 51.5 months

auc <- data.frame(trt=c(1:13), auc=NA)

auc$auc[1] <- trapz(graph_data$time, graph_data$trt1)
auc$auc[2] <- trapz(graph_data$time, graph_data$trt2)
auc$auc[3] <- trapz(graph_data$time, graph_data$trt3)
auc$auc[4] <- trapz(graph_data$time, graph_data$trt4)
auc$auc[5] <- trapz(graph_data$time, graph_data$trt5)
auc$auc[6] <- trapz(graph_data$time, graph_data$trt6)
auc$auc[7] <- trapz(graph_data$time, graph_data$trt7)
auc$auc[8] <- trapz(graph_data$time, graph_data$trt8)
auc$auc[9] <- trapz(graph_data$time, graph_data$trt9)
auc$auc[10] <- trapz(graph_data$time, graph_data$trt10)
auc$auc[11] <- trapz(graph_data$time, graph_data$trt11)
auc$auc[12] <- trapz(graph_data$time, graph_data$trt12)
auc$auc[13] <- trapz(graph_data$time, graph_data$trt13)

write.csv(auc, file="auc.csv")

