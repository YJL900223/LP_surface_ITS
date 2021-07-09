library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(graphics)
library(reshape)
library(ggsci)
pro_data <- read.csv("pro_data.csv",header = T)
data <- melt(pro_data,id.vars='names')
# Sources with filtration ratio less than 1%
data <- data[data$value>0.01,]
data$value <- round(data$value*100, digits = 1)
library(plyr)
library(RColorBrewer)

data1 <- aggregate(value ~ names + variable, data = data, sum)
(p1 <- ggplot(data=data1,aes(variable,value,fill=names))+
    geom_bar(stat="identity", position="fill",color="black", width=0.8,size=0.25)+
    scale_fill_npg()+
    labs(x = '', y = '',fill='Source(%)') +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
      legend.position = "right"
    ))
sampFile <- read.table("mapping.txt",header = T)
data_all = merge(data1, sampFile, by.x="variable", by.y = "SampleID")
data_all$time3  = factor(data_all$time3, levels=c("D58","D90","D123","D156","D216","D316","D330"))
p2 = ggplot(data_all, aes(x=time3, y = value, fill = names )) + 
  geom_bar(stat = "identity",position="fill", width=1)+ 
  scale_y_continuous(labels = scales::percent) + scale_fill_npg()+
  
  facet_grid( ~ site2, scales = "free_x", switch = "x") +  theme(strip.background = element_blank())+
  
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Site")+ylab("Percentage (%)")+ theme_classic()+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))
p2

ggsave("Fig.6 FEAST_Multi_sinks.pdf",p2,width = 8,height = 6)
