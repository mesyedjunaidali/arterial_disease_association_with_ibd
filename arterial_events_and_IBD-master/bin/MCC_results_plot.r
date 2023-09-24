setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

list.of.packages <- c("easypackages", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("easypackages")
libraries(new.packages)

ylim_low <- -0.5
ylim_upp <- 0.8

# Arterial

data_arterial <- read.csv(file="arterial_MCC_results.csv", header=TRUE, stringsAsFactor=FALSE, sep=",")

data_arterial$method <- factor(data_arterial$method, levels = data_arterial$method)

p_arterial <- ggplot(data_arterial, aes(x=method, y=mean, fill=method)) + geom_bar(stat="identity", color="black",  position=position_dodge()) + ylim(ylim_low, ylim_upp) + ylab("mean MCC ± s.d.") + xlab("") + ggtitle("arterial event classification")  +
  theme(plot.title = element_text(hjust = 0.5)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,  position=position_dodge(.9)) 
  
p_arterial <- p_arterial + scale_y_continuous(limits=c(-1, +1), breaks=c(1:10))
p_arterial <- p_arterial + theme_classic()

ggsave("arterial_MCC_barplot.pdf")

data_event_type <- read.csv(file="event_type_MCC_results.csv", header=TRUE, stringsAsFactor=FALSE, sep=",")

data_event_type$method <- factor(data_event_type$method, levels = data_arterial$method)

p_event_type <- ggplot(data_event_type, aes(x=method, y=mean, fill=method)) + geom_bar(stat="identity", color="black",  position=position_dodge()) + ylim(ylim_low, ylim_upp) + ylab("mean MCC ± s.d.") + xlab("") + ggtitle("event type classification")  +
  theme(plot.title = element_text(hjust = 0.5)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,  position=position_dodge(.9)) 

ggsave("event_type_MCC_barplot.pdf")