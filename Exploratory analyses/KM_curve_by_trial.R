# Plot KM curves for each trial

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
df1 <- data[data$studyCode==1,]
df2 <- data[data$studyCode==2,]
df3 <- data[data$studyCode==3,]
df4 <- data[data$studyCode==4,]
df5 <- data[data$studyCode==5,]
df6 <- data[data$studyCode==6,]
df7 <- data[data$studyCode==7,]
df8 <- data[data$studyCode==8,]
df9 <- data[data$studyCode==9,]
df10 <- data[data$studyCode==10,]
df11 <- data[data$studyCode==11,]
df12 <- data[data$studyCode==12,]
df13 <- data[data$studyCode==13,]

# Get median survival by printing fit
fit1 <- survfit(Surv(time, event) ~ arm, data = df1)
fit2 <- survfit(Surv(time, event) ~ arm, data = df2)
fit3 <- survfit(Surv(time, event) ~ arm, data = df3)
fit4 <- survfit(Surv(time, event) ~ arm, data = df4)
fit5 <- survfit(Surv(time, event) ~ arm, data = df5)
fit6 <- survfit(Surv(time, event) ~ arm, data = df6)
fit7 <- survfit(Surv(time, event) ~ arm, data = df7)
fit8 <- survfit(Surv(time, event) ~ arm, data = df8)
fit9 <- survfit(Surv(time, event) ~ arm, data = df9)
fit10 <- survfit(Surv(time, event) ~ arm, data = df10)
fit11 <- survfit(Surv(time, event) ~ arm, data = df11)
fit12 <- survfit(Surv(time, event) ~ arm, data = df12)
fit13 <- survfit(Surv(time, event) ~ arm, data = df13)

fit <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12, fit13)
df_list <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13)

ggsurvplot_list(fit=fit, data=df_list, censor=F)

# Using the Paired colour palette from RColorBrewer
colors=c("#070d43", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
         "#cab2d6", "#6a3d9a", "#581845", "#b15928")

p1 <- ggsurvplot(fit=fit1, data = df1, censor=F, palette=c(colors[1], colors[4]),
           xlab="Time (months)", ylab="Overall Survival", 
           xlim=c(0, 72), ylim=c(0, 1), title="BREAK-3",
           legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Dabrafenib"),
           font.legend=18)
p1$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

# save plot
setwd("R:/SCFreeman/scf20/Melanoma_NMA/Exploratory analyses")
dev.copy(pdf, "Break3.pdf")
dev.off()

p2 <- ggsurvplot(fit=fit2, data = df2, censor=F, palette=c(colors[1], colors[6]),
           xlab="Time (months)", ylab="Overall Survival", 
           xlim=c(0, 72), ylim=c(0, 1), title="BRIM-3",
           legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Vemurafenib"),
           font.legend=18)
p2$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

# save plot
dev.copy(pdf, "Brim3.pdf")
dev.off()

p3 <- ggsurvplot(fit=fit3, data = df3, censor=F, palette=c(colors[1], colors[9]),
           xlab="Time (months)", ylab="Overall Survival", 
           xlim=c(0, 72), ylim=c(0, 1), title="CheckMate 066",
           legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Nivolumab"),
           font.legend=18)
p3$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "CheckMate_066.pdf")
dev.off()

p4 <- ggsurvplot(fit=fit4, data = df4, censor=F, palette=c(colors[8], colors[9], colors[10]),
           xlab="Time (months)", ylab="Overall Survival", 
           xlim=c(0, 72), ylim=c(0, 1), title="CheckMate 067",
           legend=c(0.5,0.9), legend.title="", legend.labs=c("Ipilimumab", "Nivolumab", "Nivolumab + Ipilimumab"),
           font.legend=18)
p4$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "CheckMate_067.pdf")
dev.off()

p5 <- ggsurvplot(fit=fit5, data = df5, censor=F, palette=c(colors[8], colors[10]),
           xlab="Time (months)", ylab="Overall Survival", 
           xlim=c(0, 72), ylim=c(0, 1), title="CheckMate 069",
           legend=c(0.5,0.9), legend.title="", legend.labs=c("Ipilimumab", "Nivolumab + Ipilimumab"),
           font.legend=18)
p5$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "CheckMate_069.pdf")
dev.off()

p6 <- ggsurvplot(fit=fit6, data = df6, censor=F, palette=c(colors[6], colors[7]),
                 xlab="Time (months)", ylab="Overall Survival", 
                 xlim=c(0, 72), ylim=c(0, 1), title="COBRIM",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Vemurafenib", "Vemurafenib + Cobimetinib"),
                 font.legend=18)
p6$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "cobrim.pdf")
dev.off()

p7 <- ggsurvplot(fit=fit7, data = df7, censor=F, palette=c(colors[4], colors[5]),
                 xlab="Time (months)", ylab="Overall Survival", 
                 xlim=c(0, 72), ylim=c(0, 1), title="COMBI-d",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Dabrafenib", "Dabrafenib + Trametinib"),
                 font.legend=18)
p7$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "combi_d.pdf")
dev.off()

p8 <- ggsurvplot(fit=fit8, data = df8, censor=F, palette=c(colors[5], colors[6]),
                 xlab="Time (months)", ylab="Overall Survival", 
                 xlim=c(0, 72), ylim=c(0, 1), title="COMBI-v",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Dabrafenib + Trametinib", "Vemurafenib"),
                 font.legend=18)
p8$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "combi_v.pdf")
dev.off()

p9 <- ggsurvplot(fit=fit9, data = df9, censor=F, palette=c(colors[8], colors[12]),
                 xlab="Time (months)", ylab="Overall Survival", 
                 xlim=c(0, 72), ylim=c(0, 1), title="Hodi 2014",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Ipilimumab", "Ipilimumab + Sargramostin"),
                 font.legend=18)
p9$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "Hodi2014.pdf")
dev.off()

p10 <- ggsurvplot(fit=fit10, data = df10, censor=F, palette=c(colors[8], colors[11]),
                 xlab="Time (months)", ylab="Overall Survival", 
                 xlim=c(0, 72), ylim=c(0, 1), title="Keynote 006",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Ipilimumab", "Pembrolizumab"),
                 font.legend=18)
p10$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "Keynote006.pdf")
dev.off()

p11 <- ggsurvplot(fit=fit11, data = df11, censor=F, palette=c(colors[1], colors[2]),
                 xlab="Time (months)", ylab="Overall Survival", 
                 xlim=c(0, 72), ylim=c(0, 1), title="Ribas 2013",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Tremelimumab"),
                 font.legend=18)
p11$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "Ribas2013.pdf")
dev.off()

p12 <- ggsurvplot(fit=fit12, data = df12, censor=F, palette=c(colors[1], colors[3]),
                 xlab="Time (months)", ylab="Overall Survival", 
                 xlim=c(0, 72), ylim=c(0, 1), title="Robert 2011",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Ipilimumab + Dacarbazine"),
                 font.legend=18)
p12$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "Robert2011.pdf")
dev.off()

p13 <- ggsurvplot(fit=fit13, data = df13, censor=F, palette=c(colors[1], colors[13]),
                 xlab="Time (months)", ylab="Overall Survival", 
                 xlim=c(0, 72), ylim=c(0, 1), title="Robert 2013",
                 legend=c(0.5,0.9), legend.title="", legend.labs=c("Dacarbazine", "Selumetinib + Dacarbazine"),
                 font.legend=18)
p13$plot + scale_x_continuous(breaks=sort(c(0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72))) + scale_y_continuous(breaks=sort(c(0, 0.2, 0.4, 0.6, 0.8, 1)))

dev.copy(pdf, "Robert2013.pdf")
dev.off()

