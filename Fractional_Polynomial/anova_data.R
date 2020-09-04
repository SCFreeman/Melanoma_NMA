# Create function that takes a generated dataset and formats the data ready to apply the anova
# parameterisation

anova_data <- function(timepoints, timepoints2, ref.study=1, df){

  # Split the data at timepoints
  df2 <- survSplit(Surv(time, event) ~., data=df,
                   cut=timepoints, episode ="timegroup")
  
  # Calculate offset
  df2$y <- df2$time - df2$tstart
  
  # Add a variable that equals one for all patients - this is so the number at risk
  # can be calculated when we collapse the data
  df2$n <- 1
  
  # Collapse data
  df3 <- summaryBy(y + event + n ~ timegroup + treatment + studyCode, FUN=c(sum, max), data=df2)
  df3 <- subset(df3, select=-c(event.max, n.max))
  names(df3) <- c("spgrp", "treatment", "trialid", "y", "nevents", "natrisk", "y.max")
  
  
  # Add in a start time variable
  df3$start <- NA
  for(i in unique(df3$spgrp)){
    df3$start[df3$spgrp==i] <- timepoints2[i]
  }
  
  # Add in a time variable (i.e. how long since time 0 to max value of y for each row)
  df3$time <- df3$start + df3$y.max
  
  # Return the formatted dataset
  return(df3)

}
