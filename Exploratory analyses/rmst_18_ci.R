# Restricted mean survival time at 18 months with confidence intervals

# Load libraries
library(survRM2) 
library(foreach)

# Start with empty environment
rm(list=ls())

# Set the working directory & load data
setwd("R:/SCFreeman/scf20/Melanoma_NMA/Data")
data <- read.csv("melanoma_ipd.csv")

# Create an empty data frame to store results
rmst_18 <- data.frame(rmst0=NA, rmst0_se=NA, rmst0_lci=NA, rmst0_uci=NA,
                      rmst1=NA, rmst1_se=NA, rmst1_lci=NA, rmst1_uci=NA,
                      rmstD=NA, rmstD_lci=NA, rmstD_uci=NA)

# Two-arm studies with 1 as the reference category
id1 <- c(1, 2, 3, 11, 12, 13)

# Loop over studies
for(i in id1) {
  temp1 <- data[data$studyCode==i,]
  
  # Re-code treatments as 0 and 1
  temp1$newTx[temp1$txCode==1] <- 0
  temp1$newTx[temp1$txCode!=1] <- 1
  
  # Calculate rmst
  a <- rmst2(time=temp1$time, status=temp1$event, arm=temp1$newTx, tau=18)
  print(a)
  
  adf <- data.frame(rmst0=a$RMST.arm0$result[1,1], rmst0_se=a$RMST.arm0$result[1,2], 
                    rmst0_lci=a$RMST.arm0$result[1,3], rmst0_uci=a$RMST.arm0$result[1,4],
                    rmst1=a$RMST.arm1$result[1,1], rmst1_se=a$RMST.arm1$result[1,2],
                    rmst1_lci=a$RMST.arm1$result[1,3], rmst1_uci=a$RMST.arm1$result[1,4],
                    rmstD=a$unadjusted.result[1,1], rmstD_lci=a$unadjusted.result[1,2],
                    rmstD_uci=a$unadjusted.result[1,3])
  
  rmst_18 <- rbind(rmst_18, adf)
}
rmst_18 <- rmst_18[-1,]


# Two arm trials with 8 as the reference category
id2 <- c(5, 9, 10)

for(i in id2) {
  temp1 <- data[data$studyCode==i,]
  
  # Re-code treatments as 0 and 1
  temp1$newTx[temp1$txCode==8] <- 0
  temp1$newTx[temp1$txCode!=8] <- 1
  
  # Calculate rmst
  a <- rmst2(time=temp1$time, status=temp1$event, arm=temp1$newTx, tau=18)
  print(a)
  
  adf <- data.frame(rmst0=a$RMST.arm0$result[1,1], rmst0_se=a$RMST.arm0$result[1,2], 
                    rmst0_lci=a$RMST.arm0$result[1,3], rmst0_uci=a$RMST.arm0$result[1,4],
                    rmst1=a$RMST.arm1$result[1,1], rmst1_se=a$RMST.arm1$result[1,2],
                    rmst1_lci=a$RMST.arm1$result[1,3], rmst1_uci=a$RMST.arm1$result[1,4],
                    rmstD=a$unadjusted.result[1,1], rmstD_lci=a$unadjusted.result[1,2],
                    rmstD_uci=a$unadjusted.result[1,3])
  
  rmst_18 <- rbind(rmst_18, adf)
}

# Two arm trial with 4 as the reference category

temp1 <- data[data$studyCode==7,]

# Re-code treatments as 0 and 1
temp1$newTx[temp1$txCode==4] <- 0
temp1$newTx[temp1$txCode!=4] <- 1

# Calculate rmst
a <- rmst2(time=temp1$time, status=temp1$event, arm=temp1$newTx, tau=18)
print(a)

adf <- data.frame(rmst0=a$RMST.arm0$result[1,1], rmst0_se=a$RMST.arm0$result[1,2], 
                  rmst0_lci=a$RMST.arm0$result[1,3], rmst0_uci=a$RMST.arm0$result[1,4],
                  rmst1=a$RMST.arm1$result[1,1], rmst1_se=a$RMST.arm1$result[1,2],
                  rmst1_lci=a$RMST.arm1$result[1,3], rmst1_uci=a$RMST.arm1$result[1,4],
                  rmstD=a$unadjusted.result[1,1], rmstD_lci=a$unadjusted.result[1,2],
                  rmstD_uci=a$unadjusted.result[1,3])

rmst_18 <- rbind(rmst_18, adf)


# Two arm trial with 5 as the reference category

temp1 <- data[data$studyCode==8,]

# Re-code treatments as 0 and 1
temp1$newTx[temp1$txCode==5] <- 0
temp1$newTx[temp1$txCode!=5] <- 1

# Calculate rmst
a <- rmst2(time=temp1$time, status=temp1$event, arm=temp1$newTx, tau=18)
print(a)

adf <- data.frame(rmst0=a$RMST.arm0$result[1,1], rmst0_se=a$RMST.arm0$result[1,2], 
                  rmst0_lci=a$RMST.arm0$result[1,3], rmst0_uci=a$RMST.arm0$result[1,4],
                  rmst1=a$RMST.arm1$result[1,1], rmst1_se=a$RMST.arm1$result[1,2],
                  rmst1_lci=a$RMST.arm1$result[1,3], rmst1_uci=a$RMST.arm1$result[1,4],
                  rmstD=a$unadjusted.result[1,1], rmstD_lci=a$unadjusted.result[1,2],
                  rmstD_uci=a$unadjusted.result[1,3])

rmst_18 <- rbind(rmst_18, adf)


# Two arm trial with 6 as the reference category

temp1 <- data[data$studyCode==6,]

# Re-code treatments as 0 and 1
temp1$newTx[temp1$txCode==6] <- 0
temp1$newTx[temp1$txCode!=6] <- 1

# Calculate rmst
a <- rmst2(time=temp1$time, status=temp1$event, arm=temp1$newTx, tau=18)
print(a)

adf <- data.frame(rmst0=a$RMST.arm0$result[1,1], rmst0_se=a$RMST.arm0$result[1,2], 
                  rmst0_lci=a$RMST.arm0$result[1,3], rmst0_uci=a$RMST.arm0$result[1,4],
                  rmst1=a$RMST.arm1$result[1,1], rmst1_se=a$RMST.arm1$result[1,2],
                  rmst1_lci=a$RMST.arm1$result[1,3], rmst1_uci=a$RMST.arm1$result[1,4],
                  rmstD=a$unadjusted.result[1,1], rmstD_lci=a$unadjusted.result[1,2],
                  rmstD_uci=a$unadjusted.result[1,3])

rmst_18 <- rbind(rmst_18, adf)


# Three arm trial with 8 as the reference category
id3 <- 4

for(i in id3) {
  for(j in 9:10) {
    temp1 <- data[data$studyCode==i,]
    temp1 <- temp1[temp1$txCode==8 | temp1$txCode==j,]
    
    # Re-code treatments as 0 and 1
    temp1$newTx[temp1$txCode==8] <- 0
    temp1$newTx[temp1$txCode!=8] <- 1
    
    # Calculate rmst
    a <- rmst2(time=temp1$time, status=temp1$event, arm=temp1$newTx, tau=12)
    print(a)
    
    adf <- data.frame(rmst0=a$RMST.arm0$result[1,1], rmst0_se=a$RMST.arm0$result[1,2], 
                      rmst0_lci=a$RMST.arm0$result[1,3], rmst0_uci=a$RMST.arm0$result[1,4],
                      rmst1=a$RMST.arm1$result[1,1], rmst1_se=a$RMST.arm1$result[1,2],
                      rmst1_lci=a$RMST.arm1$result[1,3], rmst1_uci=a$RMST.arm1$result[1,4],
                      rmstD=a$unadjusted.result[1,1], rmstD_lci=a$unadjusted.result[1,2],
                      rmstD_uci=a$unadjusted.result[1,3])
    
    rmst_18 <- rbind(rmst_18, adf)
  }
}

# Look at 9 v 10 in three arm trial

temp1 <- data[data$studyCode==4,]

# Re-code treatments as 0 and 1
temp1$newTx[temp1$txCode==9] <- 0
temp1$newTx[temp1$txCode!=9] <- 1

# Calculate rmst
a <- rmst2(time=temp1$time, status=temp1$event, arm=temp1$newTx, tau=18)
print(a)

adf <- data.frame(rmst0=a$RMST.arm0$result[1,1], rmst0_se=a$RMST.arm0$result[1,2], 
                  rmst0_lci=a$RMST.arm0$result[1,3], rmst0_uci=a$RMST.arm0$result[1,4],
                  rmst1=a$RMST.arm1$result[1,1], rmst1_se=a$RMST.arm1$result[1,2],
                  rmst1_lci=a$RMST.arm1$result[1,3], rmst1_uci=a$RMST.arm1$result[1,4],
                  rmstD=a$unadjusted.result[1,1], rmstD_lci=a$unadjusted.result[1,2],
                  rmstD_uci=a$unadjusted.result[1,3])

rmst_18 <- rbind(rmst_18, adf)


# Create data frame with study details
studies <- data.frame(studyCode=c(1, 2, 3, 11, 12, 13, 5, 9, 10, 7, 8, 6, 4, 4, 4),
                      study=c("BREAK-3", "BRIM-3", "CheckMate 066", "Ribas 2013", "Robert 2011",
                              "Robert 2013", "CheckMate 069", "Hodi 2014", "Keynote 006",
                              "COMBI-d", "COMBI-v", "COBRIM", "CheckMate 067", "CheckMate 067", "CheckMate 067"),
                      t1=c(1, 1, 1, 1, 1, 1, 8, 8, 8, 4, 5, 6, 8, 8, 9),
                      t2=c(4, 6, 9, 2, 3, 13, 10, 12, 11, 5, 6, 7, 9, 10, 10))

rmst_18 <- cbind(studies, rmst_18)

rmst_18

# Store results as a csv file
setwd("../Exploratory analyses")
write.csv(rmst_18, "rmst_18_plus_ci.csv", row.names=F)
