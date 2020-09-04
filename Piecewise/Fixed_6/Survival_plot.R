# Survival plot from FE Pieceiwse model

# Start with a blank environment
rm(list=ls())

# Load libraries
library(coda)
library(survival)
library(R2WinBUGS)

# Set the root directory 
path <- 'R:/SCFreeman/scf20/Melanoma_NMA/Piecewise/Fixed_6'
setwd(path)

# Load data
pred_results <- read.csv("results.csv")

graph_data <- data.frame(time=c(0:120),
                         trt1=c(1, pred_results$mean[25:144]),
                         trt2=c(1, pred_results$mean[145:264]),
                         trt3=c(1, pred_results$mean[265:384]),
                         trt4=c(1, pred_results$mean[385:504]),
                         trt5=c(1, pred_results$mean[505:624]),
                         trt6=c(1, pred_results$mean[625:744]),
                         trt7=c(1, pred_results$mean[745:864]),
                         trt8=c(1, pred_results$mean[865:984]),
                         trt9=c(1, pred_results$mean[985:1104]),
                         trt10=c(1, pred_results$mean[1105:1224]),
                         trt11=c(1, pred_results$mean[1225:1344]),
                         trt12=c(1, pred_results$mean[1345:1464]),
                         trt13=c(1, pred_results$mean[1465:1584]))

# Colors for the plot
color=c("#070d43", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
        "#cab2d6", "#6a3d9a", "#581845", "#b15928", "#184558")

# Kaplan-Meier data for DTIC arm from CheckMate 066
ipd_data <- read.csv("../../Data/melanoma_ipd.csv")
ipd_data <- ipd_data[ipd_data$studyCode==3 & ipd_data$txCode==1,]
KM.est<-survfit(Surv(time,event)~1, data=ipd_data, type="kaplan-meier", conf.int=FALSE)

par(xpd = T, mar = par()$mar + c(6,0,0,0))


# Start by plotting the Kaplan_Meier DTIC curve
plot(KM.est,xlab="Time (months)",ylab="Overall Survival",xaxt="n",yaxt="n",main=" ",xlim=c(0,120),ylim=c(0,1),
     mark.time=FALSE, col=color[14], conf.int=F)
#Add y axis (2 specifies that axis goes on the left of the plot)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))
#Add x axis (1 specified that axis goes at the bottom of the plot)
axis(1, at=c(0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120))


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
       lty=c(1,1,1,1), ncol=3, text.width=20, box.lty=0)

# save plot
dev.copy(pdf, "survival_plot.pdf")
dev.off()

library(pracma)

# Calculate AUC at 60 months
graph60 <- graph_data[1:61, ]

auc60 <- data.frame(trt=c(1:13), auc=NA)

auc60$auc[1] <- trapz(graph60$time, graph60$trt1)
auc60$auc[2] <- trapz(graph60$time, graph60$trt2)
auc60$auc[3] <- trapz(graph60$time, graph60$trt3)
auc60$auc[4] <- trapz(graph60$time, graph60$trt4)
auc60$auc[5] <- trapz(graph60$time, graph60$trt5)
auc60$auc[6] <- trapz(graph60$time, graph60$trt6)
auc60$auc[7] <- trapz(graph60$time, graph60$trt7)
auc60$auc[8] <- trapz(graph60$time, graph60$trt8)
auc60$auc[9] <- trapz(graph60$time, graph60$trt9)
auc60$auc[10] <- trapz(graph60$time, graph60$trt10)
auc60$auc[11] <- trapz(graph60$time, graph60$trt11)
auc60$auc[12] <- trapz(graph60$time, graph60$trt12)
auc60$auc[13] <- trapz(graph60$time, graph60$trt13)

write.csv(auc60, file="auc60.csv")
