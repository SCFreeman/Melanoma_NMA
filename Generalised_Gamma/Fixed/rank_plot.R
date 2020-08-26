# Plot of probability of each treatment obtaining each rank

##################################################################################

# Load libraries
library(reshape2)
library(ggplot2)

# Emty environment
rm(list=ls())

# Set working directory
path <- "R:/SCFreeman/scf20/Melanoma_NMA/Generalised_Gamma/Fixed"
setwd(path)

# Read in results
data <- read.csv("FE_results.csv")

# Keep rows for ranking only
data <- data[25:193,]

# Variable for rank
data$rank_code <- rep(13:1, 13)

# Restrict probability to 2 decimal places
data$prob <- round(data$mean, 2)

# Add a treatment label
data$Treatment <- c(rep("DTIC", 13), rep("TRL", 13), rep("IPI+DTIC", 13), rep("DB", 13), rep("DB+TR", 13),
                    rep("VM", 13), rep("VM+COB", 13), rep("IPI", 13), rep("NIV", 13), rep("NIV+IPI", 13),
                    rep("PEM", 13), rep("IPI+SRG", 13), rep("SEL+DTIC", 13))

# Rename mean column
names(data)[names(data)=="mean"] <- "Probability"

# Plot with text
q <- ggplot(data, aes(x=rank_code, y=Treatment)) +
  geom_point(aes(size=Probability), shape=21, colour="skyblue", fill="skyblue") +
  theme(panel.background=element_blank(), panel.border=element_rect(colour="black", fill=NA, size=1),
        legend.position="bottom") +
  scale_size_area(max_size=10) +
  scale_x_continuous(name="Rank", limits=c(1, 13), breaks=seq(1,13,1)) +
  scale_y_discrete(name="Treatment") +  
  geom_text(aes(label=prob))
q

# save plot
dev.copy(pdf, "rank_plot.pdf")
dev.off()
