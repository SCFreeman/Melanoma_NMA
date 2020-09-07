# Generalised gamma model with treatment modelled as a location and scale parameter

##################################################################################

# Start with a blank environment
rm(list=ls())

# Load packages
library(survival)
library(flexsurv)
library(R2WinBUGS)

# Set the working directory
path <- 'R:/SCFreeman/scf20/Melanoma_NMA/Generalised_Gamma/Fixed_scale'
setwd(path)

# load data
data <- read.csv("../../Data/melanoma_ipd.csv")
data <- data[, 2:8]

# Set the location for WinBUGS
bugs.directory <- "Z:\\winbugs143_unrestricted\\winbugs14_full_patched\\WinBUGS14"

# WinBUGS burn-in & simulation size (this will be divided by the number of chains)
num.sims <- 60000
burn.in <- 30000

#-----------------------------------------------------------------------------
# Initial data formatting
#-----------------------------------------------------------------------------

# Add treatment as a character variable
data$trt[data$txCode==1] <- "DTIC"
data$trt[data$txCode==2] <- "TRL"
data$trt[data$txCode==3] <- "IPI+DTIC"
data$trt[data$txCode==4] <- "DB"
data$trt[data$txCode==5] <- "DB+TR"
data$trt[data$txCode==6] <- "VM"
data$trt[data$txCode==7] <- "VM+COB"
data$trt[data$txCode==8] <- "IPI"
data$trt[data$txCode==9] <- "NIV"
data$trt[data$txCode==10] <- "NIV+IPI"
data$trt[data$txCode==11] <- "PEM"
data$trt[data$txCode==12] <- "IPI+SRG"
data$trt[data$txCode==13] <- "SEL+DTIC"

#select data set
selected <- data.frame(trt=data$trt, study=data$study, time=data$time, event=data$event, 
                       studyCode=data$studyCode, txCode=data$txCode)

# Convert trt and study to factor variables
selected$trt <- as.factor(selected$trt)
selected$study <- as.factor(selected$study)

# Re-order trt so that the treatments are in the order we want (rather than alphabetical)
selected$trt=factor(selected$trt, levels(selected$trt)[c(3, 11, 5, 1, 2, 12, 13, 4, 7, 8, 9, 6, 10)])

# Vector of study names
Studies <- levels(selected$study)

# Vector of treatment names
Treat <- levels(selected$trt)

# Number of studies
nStudies<-length(Studies)

# Number of treatments
nTreat<-length(Treat)

# Select the greatest f/up time and round it to the nearest 10
maxFU<-round(max(selected$time),-1)


#-----------------------------------------------------------------------------
# Two-stage BUGS Analysis: First stage - Analyse each trial individually and 
# create a data frame containing results from each trial
#-----------------------------------------------------------------------------

# Create data set for winbugs NMA 
study_GG<-rep(NA, nStudies)

# Loop over the number of studies
for (i in 1:nStudies){
  # Select one study
  data_i<-selected[selected$study==Studies[i],]
  # Vector of treatments used in study 
  data_i$trt<-factor(data_i$trt)
  #Fit generalised gamma model with trt modelled on location and scale parameters
  nmastud<-flexsurvreg(formula=Surv(time,event)~trt + sigma(trt),
                       data=data_i, dist="gengamma", method="BFGS")
  # Identify treatments other than the reference
  labels<-paste("trt",levels(data_i$trt)[2:length(levels(data_i$trt))],
                sep="")
  labels2<-paste("sigma(trt",levels(data_i$trt)[2:length(levels(data_i$trt))], ")",
                sep="")
  #Number of treatment contrasts in study
  contrasts<-length(labels)
  # Save results of beta coefficients and standard errors
  beta<-nmastud$coefficients[labels]
  se<-sqrt(diag(nmastud$cov))[labels]
  sigma <- nmastud$coefficients[labels2]
  sigma.se <- sqrt(diag(nmastud$cov))[labels2]
  #Set up covariance matrix but leave empty if only two treatments in study
  cov<-NA
  cov_s <- NA
  if (contrasts>1) {cov<-sqrt(nmastud$cov[labels,labels][1,2])}
  if (contrasts>1) {cov_s<-sqrt(nmastud$cov[labels2,labels2][1,2])}
  # Add results for study to a data frame 
  study_i<-data.frame(STUDY=rep(Studies[i],contrasts),
                      COMP=gsub("trt","",labels),
                      REF=rep(levels(data_i$trt)[1],contrasts),
                      MEAN=beta,
                      MEANSE=se,
                      sigma=sigma,
                      sigma.se = sigma.se,
                      COV=cov,
                      COV_s=cov_s)
  # Add results for most recent study to a data frame which contains the results of previous studies
  study_GG<-rbind(study_GG,study_i)
}
study_GG<-study_GG[-1,]


#---------------------------------
# Correction for multi-arm trials
#--------------------------------

study_GG$multi<-0

# Create row vector which contains FALSE if a two-arm trial and TRUE if a multi-arm trial
multi_index<-study_GG$STUDY%in%rownames(table(study_GG$STUDY))[table(study_GG$STUDY)>1]

# Use this vector to replace 0's with 1's in the column multi if multi_index is TRUE
study_GG[multi_index,]$multi<-1
#Adjust the variable MEANSE to take into account multi-arm trials
study_GG[multi_index,"MEANSE"]<-sqrt(study_GG[multi_index,"MEANSE"]^2-study_GG[multi_index,"COV"]^2)
study_GG[multi_index,"sigma.se"]<-sqrt(study_GG[multi_index,"sigma.se"]^2-study_GG[multi_index,"COV_s"]^2)

#Identify multi-arm trials
multi_study<-unique(study_GG$STUDY[study_GG$multi==1])
#Identify first row in study_GG which belongs to the multi-arm trial
multi_first<-match(multi_study,study_GG$STUDY)

# Add a row to the data frame for each multi-arm trial
study_GG<-rbind(study_GG,
                data.frame(STUDY=multi_study,
                            COMP=study_GG$REF[multi_first],
                            REF=study_GG$REF[multi_first],
                            MEAN=rep(0,length(multi_study)),
                            MEANSE=study_GG$COV[multi_first],
                           sigma=rep(0,length(multi_study)),
                           sigma.se=study_GG$COV_s[multi_first],
                            COV=rep(NA,length(multi_study)),
                           COV_s=rep(NA,length(multi_study)),
                            multi=rep(1,length(multi_study))
                           )   
                )

#Remove the row names (turns them back to numbers)
row.names(study_GG)<-NULL

# Make study a factor variable
study_GG$studyCode <- as.factor(study_GG$STUDY)

#-----------------------------------------------------------------------------
# Prepare data for analysis in WinBUGS
#-----------------------------------------------------------------------------

#Use the list command to group together data needed for WinBUGS
BUGS_data<-list(LnObs=dim(study_GG)[1], nTx=nTreat, nStudies=nStudies, Lstudy=as.numeric(study_GG$studyCode),
                Ltx=match(study_GG$COMP,Treat), Lbase=match(study_GG$REF,Treat),
                Lmean=study_GG[,"MEAN"], Lse=study_GG[,"MEANSE"], multi=study_GG$multi,
                Lmean2=study_GG[,"sigma"], Lse2=study_GG[,"sigma.se"])

# Set up initial values for FE model
initsFE1<-list(alpha = rep(-0.5,nStudies), beta = c(NA,rep(-0.5,nTreat-1)), gamma = c(NA,rep(-0.5,nTreat-1)))
initsFE2<-list(alpha = rep(0.5,nStudies), beta = c(NA,rep(0.5,nTreat-1)), gamma = c(NA,rep(0.5,nTreat-1)))
initsFE3<-list(alpha = rep(0.1,nStudies), beta = c(NA,rep(0.1,nTreat-1)), gamma = c(NA,rep(0.1,nTreat-1)))
initsFE<-list(initsFE1, initsFE2, initsFE3)

#-----------------------------------------------------------------------------
# Fixed Effect Model for treatment & scale effect
#-----------------------------------------------------------------------------

# Run WinBUGS model 

nmaAggrFE<-bugs(data=BUGS_data, inits=initsFE, parameters.to.save=c("beta", "aft", "rk.beta", "gamma", 
                                                                    "rk.gamma", "scale.aft"),
                model.file="FE_model.txt", clearWD=F, summary.only=FALSE, n.iter=(num.sims+burn.in), 
                n.sims=num.sims, n.burnin=burn.in, n.chains=3, 
                debug=FALSE,  bugs.seed=385916, bugs.directory = bugs.directory, working.directory=path)

# Save results in csv file
FE_results <- nmaAggrFE$summary
write.csv(FE_results,file="FE_results.csv")

# Check results
nmaAggrFEbetas<-nmaAggrFE$sims.matrix[,grep("beta",rownames(nmaAggrFE$summary))]
nmaAggrFEbetas<-cbind(rep(0,dim(nmaAggrFEbetas)[1]),nmaAggrFEbetas)
colnames(nmaAggrFEbetas)<-Treat
summary(nmaAggrFEbetas)

coda.options(combine.plots=FALSE)
nmaAggrFE_mcmc<-mcmc(nmaAggrFEbetas)

all(nmaAggrFE$summary[,"Rhat"] < 1.1) 
print(nmaAggrFE,digits=3)

# Check trace, density and autocorrelation plots
par(mfrow=c(3,2))

autocorr.plot(nmaAggrFE_mcmc)

traceplot(nmaAggrFE_mcmc[,2])
traceplot(nmaAggrFE_mcmc[,3])
traceplot(nmaAggrFE_mcmc[,4])
traceplot(nmaAggrFE_mcmc[,5])
traceplot(nmaAggrFE_mcmc[,6])
traceplot(nmaAggrFE_mcmc[,7])
traceplot(nmaAggrFE_mcmc[,8])
traceplot(nmaAggrFE_mcmc[,9])
traceplot(nmaAggrFE_mcmc[,10])
traceplot(nmaAggrFE_mcmc[,11])
traceplot(nmaAggrFE_mcmc[,12])
traceplot(nmaAggrFE_mcmc[,13])

densplot(nmaAggrFE_mcmc[,2])
densplot(nmaAggrFE_mcmc[,3])
densplot(nmaAggrFE_mcmc[,4])
densplot(nmaAggrFE_mcmc[,5])
densplot(nmaAggrFE_mcmc[,6])
densplot(nmaAggrFE_mcmc[,7])
densplot(nmaAggrFE_mcmc[,8])
densplot(nmaAggrFE_mcmc[,9])
densplot(nmaAggrFE_mcmc[,10])
densplot(nmaAggrFE_mcmc[,11])
densplot(nmaAggrFE_mcmc[,12])
densplot(nmaAggrFE_mcmc[,13])

#-----------------------------------------------------------------------------
# Survival prediction - CheckMate 066
#-----------------------------------------------------------------------------

# Fit generalised gamma model just for the DTIC arm of the CheckMate 066 trial
df <- data[data$studyCode==3 & data$txCode==1,]

ma <- flexsurvreg(formula=Surv(time,event)~1,
                  data=df, dist="gengamma", method="BFGS")

# Identify coefficicents needed for predicting survival
mu <- ma$coefficients["mu"]
sigma <- exp(ma$coefficients["sigma"])
q <- ma$coefficients["Q"]

# Store trt effects (beta) for each treatment
trt2 <- nmaAggrFE$summary[1,1]
trt3 <- nmaAggrFE$summary[2,1]
trt4 <- nmaAggrFE$summary[3,1]
trt5 <- nmaAggrFE$summary[4,1]
trt6 <- nmaAggrFE$summary[5,1]
trt7 <- nmaAggrFE$summary[6,1]
trt8 <- nmaAggrFE$summary[7,1]
trt9 <- nmaAggrFE$summary[8,1]
trt10 <- nmaAggrFE$summary[9,1]
trt11 <- nmaAggrFE$summary[10,1]
trt12 <- nmaAggrFE$summary[11,1]
trt13 <- nmaAggrFE$summary[12,1]

# Store sigma effects (gamma) for each treatment
sigma.trt2 <- nmaAggrFE$summary[194,1]
sigma.trt3 <- nmaAggrFE$summary[195,1]
sigma.trt4 <- nmaAggrFE$summary[196,1]
sigma.trt5 <- nmaAggrFE$summary[197,1]
sigma.trt6 <- nmaAggrFE$summary[198,1]
sigma.trt7 <- nmaAggrFE$summary[199,1]
sigma.trt8 <- nmaAggrFE$summary[200,1]
sigma.trt9 <- nmaAggrFE$summary[201,1]
sigma.trt10 <- nmaAggrFE$summary[202,1]
sigma.trt11 <- nmaAggrFE$summary[203,1]
sigma.trt12 <- nmaAggrFE$summary[204,1]
sigma.trt13 <- nmaAggrFE$summary[205,1]


#Calculate survival across 120 months
x <- seq(0,120,0.1)
S.trt1 <- 1-pgengamma(x, mu = mu, sigma = sigma, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt2 <- 1-pgengamma(x, mu = mu+trt2, sigma = sigma+sigma.trt2, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt3 <- 1-pgengamma(x, mu = mu+trt3, sigma = sigma+sigma.trt3, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt4 <- 1-pgengamma(x, mu = mu+trt4, sigma = sigma+sigma.trt4, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt5 <- 1-pgengamma(x, mu = mu+trt5, sigma = sigma+sigma.trt5, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt6 <- 1-pgengamma(x, mu = mu+trt6, sigma = sigma+sigma.trt6, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt7 <- 1-pgengamma(x, mu = mu+trt7, sigma = sigma+sigma.trt7, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt8 <- 1-pgengamma(x, mu = mu+trt8, sigma = sigma+sigma.trt8, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt9 <- 1-pgengamma(x, mu = mu+trt9, sigma = sigma+sigma.trt9, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt10 <- 1-pgengamma(x, mu = mu+trt10, sigma = sigma+sigma.trt10, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt11 <- 1-pgengamma(x, mu = mu+trt11, sigma = sigma+sigma.trt11, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt12 <- 1-pgengamma(x, mu = mu+trt12, sigma = sigma+sigma.trt12, Q=q, lower.tail = TRUE, log.p = FALSE)
S.trt13 <- 1-pgengamma(x, mu = mu+trt13, sigma = sigma+sigma.trt13, Q=q, lower.tail = TRUE, log.p = FALSE)

graph_data <- data.frame(time=x, trt1=S.trt1, trt2=S.trt2, trt3=S.trt3, trt4=S.trt4, trt5=S.trt5, trt6=S.trt6,
                         trt7=S.trt7, trt8=S.trt8, trt9=S.trt9, trt10=S.trt10, trt11=S.trt11, 
                         trt12 = S.trt12, trt13=S.trt13)
write.csv(graph_data, "graph_data.csv")

# Calculate KM estimate for DTIC from CheckMate 066
ipd_data <- df[df$txCode==1,]
KM.est<-survfit(Surv(time,event)~1, data=ipd_data, type="kaplan-meier", conf.int=FALSE)

# Using the Paired colour palette from RColorBrewer
color=c("#070d43", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
        "#cab2d6", "#6a3d9a", "#581845", "#b15928", "#184558")

# Graph region parameters
par(xpd = T, mar = par()$mar + c(6,0,0,0))

# Start by plotting the Kaplan-Meier DTIC curve for the CheckMate 066 trial
plot(KM.est,xlab="Time (months)",ylab="Overall Survival",xaxt="n",yaxt="n",main=" ",xlim=c(0,120),ylim=c(0,1),
     mark.time=FALSE, col=color[14], conf.int=F)

#Add y axis (2 specifies that axis goes on the left of the plot)
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1))

#Add x axis (1 specified that axis goes at the bottom of the plot)
axis(1, at=c(0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120))

# Add survival curves for each treatment
lines(x,S.trt1, col=color[1])
lines(x,S.trt2, col=color[2])
lines(x,S.trt3, col=color[3])
lines(x,S.trt4, col=color[4])
lines(x,S.trt5, col=color[5])
lines(x,S.trt6, col=color[6])
lines(x,S.trt7, col=color[7])
lines(x,S.trt8, col=color[8])
lines(x,S.trt9, col=color[9])
lines(x,S.trt10, col=color[10])
lines(x,S.trt11, col=color[11])
lines(x,S.trt12, col=color[12])
lines(x,S.trt13, col=color[13])

# Add legend
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
graph60 <- graph_data[1:601, ]
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