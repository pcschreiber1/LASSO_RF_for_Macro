library(ggplot2)
library(dplyr)

setwd("C:/Projects/Comp_Stat_Project")

sim1 <- read.csv("sim1.csv", header=TRUE)
sim1 = sim1[1:40,]

df <- sim1 %>% group_by(SNR, Method) %>%  summarize(Mean_Ret = mean(Retention, na.rm=TRUE),
                                                    Mean_Zero = mean(Nonzero, na.rm=TRUE),
                                                    Mean_Pred = mean(Prediction, na.rm=TRUE))
snr.breaks = round(exp(seq(from=min(log(sim1$SNR)),
                           to=max(log(sim1$SNR)),length=4)),2)

ggplot(data=df, aes(x=SNR, y=Mean_Ret, color=Method)) +
  geom_line(lwd=1) +
  geom_point(pch=19) +
  theme_bw() +
  #facet_grid(rows = vars(Method)) +
  #facet_grid(formula(paste(1,"~",2))) +
  xlab("Signal-to-noise ratio") +
  ylab("Retention") +
  geom_line(aes(x=SNR, y=5), lwd=0.5, linetype=3, color="black") +
  ggtitle("Simulation 1") + 
  scale_x_continuous(trans="log", breaks=snr.breaks)
  
  