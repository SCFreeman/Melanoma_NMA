# Plot results from RP fixed effect model

# Empty environment
rm(list=ls())

# Load library
library(survival)

# Working directory
setwd("R:/SCFreeman/scf20/Melanoma_NMA/Royston_Parmar/FE_nonPH")

# Load results
pred_results <- read.csv("results.csv")

graph_data <- data.frame(time=c(0:120),
                         trt1=c(1, pred_results$mean[246:365]),
                         trt2=c(1, pred_results$mean[366:485]),
                         trt3=c(1, pred_results$mean[486:605]),
                         trt4=c(1, pred_results$mean[606:725]),
                         trt5=c(1, pred_results$mean[726:845]),
                         trt6=c(1, pred_results$mean[846:965]),
                         trt7=c(1, pred_results$mean[966:1085]),
                         trt8=c(1, pred_results$mean[1086:1205]),
                         trt9=c(1, pred_results$mean[1206:1325]),
                         trt10=c(1, pred_results$mean[1326:1445]),
                         trt11=c(1, pred_results$mean[1446:1565]),
                         trt12=c(1, pred_results$mean[1566:1685]),
                         trt13=c(1, pred_results$mean[1686:1805]))

# Colors for the plot
color=c("#070d43", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
        "#cab2d6", "#6a3d9a", "#581845", "#b15928", "#184558")

# Kaplan-Meier for DTIC arm from CheckMate 066
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
dev.copy(pdf, "survival_plot_10yrs.pdf")
dev.off()
