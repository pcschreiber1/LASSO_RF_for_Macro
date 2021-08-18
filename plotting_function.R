library(ggplot2)
library(dplyr)
library(cowplot)

setwd("C:/Projects/Comp_Stat_Project")

sim1 <- read.csv("sim1.csv", header=TRUE)

#Prepare the data
df <- sim1 %>% #create new data frame
  na_if(Inf) %>% #Change INF values produced by RF
  group_by(SNR, Method) %>%  
  summarize(Mean_Ret = mean(Retention, na.rm=TRUE),
            Mean_Zero = mean(Nonzero, na.rm=TRUE),
            Mean_Pred = mean(Prediction, na.rm=TRUE))

# Line breaks for SNR (logarithmic scale)
snr.breaks = round(exp(seq(from=min(log(sim1$SNR)),
                           to=max(log(sim1$SNR)),length=4)),2)

p1 <- ggplot(data=df, aes(x=SNR, y=Mean_Ret, color=Method)) +
  geom_line(lwd=1) +
  geom_point(pch=19) +
  theme_bw() +
  #facet_grid(rows = vars(Method)) +
  #facet_grid(formula(paste(1,"~",2))) +
  xlab("Signal-to-noise ratio") +
  ylab("Retention") +
  geom_line(aes(x=SNR, y=5), lwd=0.5, linetype=3, color="black") +
  #ggtitle("Simulation 1") + 
  scale_x_continuous(trans="log", breaks=snr.breaks)

p2 <- ggplot(data=df, aes(x=SNR, y=Mean_Zero, color=Method)) +
  geom_line(lwd=1) +
  geom_point(pch=19) +
  theme_bw() +
  #facet_grid(rows = vars(Method)) +
  #facet_grid(formula(paste(1,"~",2))) +
  xlab("Signal-to-noise ratio") +
  ylab("Number of Nonzero Coeff") +
  geom_line(aes(x=SNR, y=5), lwd=0.5, linetype=3, color="black") +
  #ggtitle("Simulation 1") + 
  scale_x_continuous(trans="log", breaks=snr.breaks)

p3 <- ggplot(data=df, aes(x=SNR, y=Mean_Pred, color=Method)) +
  geom_line(lwd=1) +
  geom_point(pch=19) +
  theme_bw() +
  #facet_grid(rows = vars(Method)) +
  #facet_grid(formula(paste(1,"~",2))) +
  xlab("Signal-to-noise ratio") +
  ylab("Mean-squared Prediction Error") +
  geom_line(aes(x=SNR, y=5), lwd=0.5, linetype=3, color="black") +
  #ggtitle("Simulation 1") + 
  scale_x_continuous(trans="log", breaks=snr.breaks)

# get legend manually
legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
plot_grid(p1 +   theme(legend.position="none"),
          p2 +   theme(legend.position="none"),
          p3 +   theme(legend.position="none"),
          legend,
          ncol=2,
          nrow=2)
