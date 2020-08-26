# Test for non-PH

# Load libraries
library(survival)
library(broom)
library(metafor)
library(survminer)

# Start with an empty environment
rm(list=ls())

# Set the working directory
path <- "R:/SCFreeman/scf20/Melanoma_NMA/Data/"
setwd(path)

# Import the data and split into data frames for each trial
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

# Fit the Cox PH model
cox1 <- coxph(formula = Surv(time, event) ~ arm, data=df1)
cox2 <- coxph(formula = Surv(time, event) ~ arm, data=df2)
cox3 <- coxph(formula = Surv(time, event) ~ arm, data=df3)
cox4 <- coxph(formula = Surv(time, event) ~ arm, data=df4)
cox5 <- coxph(formula = Surv(time, event) ~ arm, data=df5)
cox6 <- coxph(formula = Surv(time, event) ~ arm, data=df6)
cox7 <- coxph(formula = Surv(time, event) ~ arm, data=df7)
cox8 <- coxph(formula = Surv(time, event) ~ arm, data=df8)
cox9 <- coxph(formula = Surv(time, event) ~ arm, data=df9)
cox10 <- coxph(formula = Surv(time, event) ~ arm, data=df10)
cox11 <- coxph(formula = Surv(time, event) ~ arm, data=df11)
cox12 <- coxph(formula = Surv(time, event) ~ arm, data=df12)
cox13 <- coxph(formula = Surv(time, event) ~ arm, data=df13)

# Test the PH assumption
test1 <- cox.zph(cox1)
test2 <- cox.zph(cox2)
test3 <- cox.zph(cox3)
test4 <- cox.zph(cox4)
test5 <- cox.zph(cox5)
test6 <- cox.zph(cox6)
test7 <- cox.zph(cox7)
test8 <- cox.zph(cox8)
test9 <- cox.zph(cox9)
test10 <- cox.zph(cox10)
test11 <- cox.zph(cox11)
test12 <- cox.zph(cox12)
test13 <- cox.zph(cox13)

# Graphical test of PH
ggcoxzph(test1)
ggcoxzph(test2)
ggcoxzph(test3)
ggcoxzph(test4)
ggcoxzph(test5)
ggcoxzph(test6)
ggcoxzph(test7)
ggcoxzph(test8)
ggcoxzph(test9)
ggcoxzph(test10)
ggcoxzph(test11)
ggcoxzph(test12)
ggcoxzph(test13)

# Log-log plots
s1 <- Surv(df1$time, df1$event)
s2 <- Surv(df2$time, df2$event)
s3 <- Surv(df3$time, df3$event)
s4 <- Surv(df4$time, df4$event)
s5 <- Surv(df5$time, df5$event)
s6 <- Surv(df6$time, df6$event)
s7 <- Surv(df7$time, df7$event)
s8 <- Surv(df8$time, df8$event)
s9 <- Surv(df9$time, df9$event)
s10 <- Surv(df10$time, df10$event)
s11 <- Surv(df11$time, df11$event)
s12 <- Surv(df12$time, df12$event)
s13 <- Surv(df13$time, df13$event)

plot(survfit(s1 ~ df1$txCode), col=c("blue", "red"), fun="cloglog", 
     xlab="Log(time)", ylab="log(-log(survival probability))")
plot(survfit(s2 ~ df2$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s3 ~ df3$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s4 ~ df4$txCode), col=c("blue", "red", "green"), fun="cloglog")
plot(survfit(s5 ~ df5$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s6 ~ df6$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s7 ~ df7$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s8 ~ df8$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s9 ~ df9$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s10 ~ df10$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s11 ~ df11$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s12 ~ df12$txCode), col=c("blue", "red"), fun="cloglog")
plot(survfit(s13 ~ df13$txCode), col=c("blue", "red"), fun="cloglog")

