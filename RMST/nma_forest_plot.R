#### Create forest plot of RMST results

# Load library
library(metafor)

# Start with empty environment
rm(list=ls())

# Set working directory and load data
setwd("R:/SCFreeman/scf20/Melanoma_NMA/RMST")
data <- read.csv("fe_results.csv") 
data <- data[1:12,]

data$trt <- c("TRL", "IPI + DTIC", "DB", "DB+TR", "VM", "VM+COB", "IPI", "NIV", "NIV+IPI", "PEM", 
              "IPI+SRG", "SEL+DTIC")

forest(x=data$mean, ci.lb=data$X2.5, ci.ub=data$X97.5,
       slab = data$trt,
       refline=0, top=1, xlab="RMST (months)", cex=0.8)

dev.copy(pdf, "nma_forest_plot.pdf")
dev.off()
