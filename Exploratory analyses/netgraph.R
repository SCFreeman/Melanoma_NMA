# Create a network plot

# Load library
library(netmeta)

# Start with empty environment
rm(list=ls())

# Set the workign directory
setwd("R:/SCFreeman/scf20/Melanoma_NMA/Exploratory analyses")

# Load the data
data <- read.csv("netgraph_data.csv")

# Create a netmeta object
a <- netmeta(data$lhr, data$se, treat1=data$t1, treat2=data$t2, studlab=data$Study, reference=1)
netgraph(a)

# Treatment labels
lab <- c("DTIC", "SEL+DTIC", "IPI+DTIC", "DB", "DB+TR", "VM", "VM+COB", "NIV", "NIV+IPI", "IPI",
         "PEM", "IPI+SRG", "TRL")

netgraph(a, labels=lab, offset=0.02, plastic=F, col="#434343", multiarm=T, col.multiarm="purple", points=T,
         col.points="blue", number.of.studies = T, cex=1.5,
         cex.points=c(6,1,1,2,2,3,1,3,3,5,1,1,1))

dev.copy(pdf, "netgraph.pdf")
dev.off()
