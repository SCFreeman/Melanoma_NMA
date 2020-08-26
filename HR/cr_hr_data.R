# HR for each trial

###############################################################

library(survival)
library(broom)
library(metafor)

rm(list=ls())

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

# 2 arm trials
cox1 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df1), exp=T)
cox2 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df2), exp=T)
cox3 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df3), exp=T)
cox5 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df5), exp=T)
cox6 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df6), exp=T)
cox7 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df7), exp=T)
cox8 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df8), exp=T)
cox9 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df9), exp=T)
cox10 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df10), exp=T)
cox11 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df11), exp=T)
cox12 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df12), exp=T)
cox13 <- broom::tidy(coxph(formula = Surv(time, event) ~ arm, data=df13), exp=T)


# 3 arm trial - need to convert arm to a factor variable
df4$treat <- factor(df4$arm, labels=c("A", "B", "C"))
coxph(formula = Surv(time, event) ~ treat, data=df4)
cox4 <- broom::tidy(coxph(formula = Surv(time, event) ~ treat, data=df4), exp=T)

# Vector of HR
hr <- c(cox1$estimate, cox2$estimate, cox3$estimate, cox5$estimate, cox6$estimate,
        cox7$estimate, cox8$estimate, cox9$estimate, cox10$estimate, cox11$estimate, cox12$estimate,
        cox13$estimate, cox4$estimate[1], cox4$estimate[2])

# Vector of SE
se <- c(cox1$std.error, cox2$std.error, cox3$std.error, cox5$std.error, 
        cox6$std.error, cox7$std.error, cox8$std.error, cox9$std.error, cox10$std.error, cox11$std.error, 
        cox12$std.error, cox13$std.error, cox4$std.error[1], cox4$std.error[2])

# Vector of LCI
lci <- c(cox1$conf.low, cox2$conf.low, cox3$conf.low, cox5$conf.low, 
         cox6$conf.low, cox7$conf.low, cox8$conf.low, cox9$conf.low, cox10$conf.low, cox11$conf.low, 
         cox12$conf.low, cox13$conf.low, cox4$conf.low[1], cox4$conf.low[2])

# Vector of UCI
uci <- c(cox1$conf.high, cox2$conf.high, cox3$conf.high, cox5$conf.high, 
         cox6$conf.high, cox7$conf.high, cox8$conf.high, cox9$conf.high, cox10$conf.high, cox11$conf.high, 
         cox12$conf.high, cox13$conf.high, cox4$conf.high[1], cox4$conf.high[2])

# Calculate LHR
lhr <- log(hr)

# Create a data frame
hr_data <- data.frame(study=c("BREAK-3", "BRIM-3", "CheckMate 066", "CheckMate 069", 
                              "COBRIM", "COMBI-d", "COMBI-v", "Hodi 2014", "Keynote 006", "Ribas 2013", "Robert 2011",
                              "Robert 2013", "CheckMate 067", "CheckMate 067"),
                      treat1=c(1, 1, 1, 8, 6, 4, 5, 8, 8, 1, 1, 1, 8, 8),
                      treat2=c(4, 6, 9, 10, 7, 5, 6, 12, 11, 2, 3, 13, 9, 10),
                      hr=hr,
                      se=se, lci=lci, uci=uci, lhr=lhr)

# Calculate variance of baseline treatment
sd1 <- sd(df4$time[df4$txCode==8])
sd2 <- sd(df4$time[df4$txCode==9])
n1 <- 315
n2 <- 316

swithin <- (((n1-1)*sd1^2)+((n2-1)*sd2^2))/((n1+n2-2)^0.5)

V <- (sd1^2/n1) / (swithin^2)

hr_data$V <- NA
hr_data$V[13:14] <- V

# save as csv file
path <- "R:/SCFreeman/scf20/Melanoma_NMA/HR/"
setwd(path)
write.csv(hr_data, "hr_data.csv")
