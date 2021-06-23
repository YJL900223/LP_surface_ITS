rm(list = ls())
gc()


source("src//src.R")

#Set the arguments of your data
metadata_file = "Data_files//metadata_example_multi.txt"
count_matrix = "Data_files//otu_example_multi.txt"
EM_iterations = 1000 #default value
##if you use different sources for each sink, different_sources_flag = 1, otherwise = 0
different_sources_flag = 1


# Load sample metadata
metadata <- read.csv(metadata_file,h=T, sep = "\t", row.names = 1)

# Load OTU table
otus <- read.table(count_matrix, header = T, comment = '', check = F, sep = '\t',row.names = 1)
otus <- otus[,-(which(colSums(otus)<100))]

otus <- t(as.matrix(otus))


# Extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# Double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}


if(different_sources_flag == 0){
  
  metadata$id[metadata$SourceSink == 'Source'] = NA
  metadata$id[metadata$SourceSink == 'Sink'] = c(1:length(which(metadata$SourceSink == 'Sink')))
}


envs <- metadata$Env
Ids <- na.omit(unique(metadata$id))
Proportions_est <- list()


for(it in 1:length(Ids)){
  
  
  # Extract the source environments and source/sink indices
  if(different_sources_flag == 1){
    
    train.ix <- which(metadata$SourceSink=='Source' )
    test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])
    
  }
  
  else{
    
    train.ix <- which(metadata$SourceSink=='Source')
    test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])
  }
  
  num_sources <- length(train.ix)
  COVERAGE =  min(rowSums(otus[c(train.ix, test.ix),]))  #Can be adjusted by the user
  str(COVERAGE)
  # Define sources and sinks
  
  sources <- as.matrix(rarefy(otus[train.ix,], COVERAGE))
  sinks <- as.matrix(rarefy(t(as.matrix(otus[test.ix,])), COVERAGE))
  
  
  print(paste("Number of OTUs in the sink sample = ",length(which(sinks > 0))))
  print(paste("Seq depth in the sources and sink samples = ",COVERAGE))
  print(paste("The sink is:", envs[test.ix]))
  
  # Estimate source proportions for each sink
  
  FEAST_output<-FEAST(source=sources, sinks = t(sinks), env = envs[train.ix], em_itr = EM_iterations, COVERAGE = COVERAGE)
  Proportions_est[[it]] <- FEAST_output$data_prop[,1]
  
  
  names(Proportions_est[[it]]) <- c(as.character(envs[train.ix]), "unknown")
  
  if(length(Proportions_est[[it]]) < num_sources +1){
    
    tmp = Proportions_est[[it]]
    Proportions_est[[it]][num_sources] = NA
    Proportions_est[[it]][num_sources+1] = tmp[num_sources]
  }
  
  print("Source mixing proportions")
  print(Proportions_est[[it]])
  
  
}

print(Proportions_est)
#可视化过程，参考文涛脚本并修正若干bug
#输出计算结果
FEAST_output = as.data.frame(Proportions_est)
colnames(FEAST_output) = paste("repeat_",Ids,sep = "") #取Ids作为平行代号
head(FEAST_output)


write.csv(FEAST_output,"FEAST.csv",quote = F)


#简单出图(每个repeat一张)
library(RColorBrewer)
library(dplyr)
library(graphics)
head(FEAST_output)

pdf(file = "FEAST_repeat.pdf",width = 12,height = 12)
par(mfrow=c((length(unique(metadata$SampleType))%/%2 +2 ),2), mar=c(1,1,1,1))
# layouts = as.character(unique(metadata$SampleType))

for (i in 1:length(colnames(FEAST_output))) {
  
  labs <- paste0(row.names(FEAST_output)," (", round(FEAST_output[,i]/sum(FEAST_output[,i])*100,2), "%)")
  
  pie(FEAST_output[,i],labels=labs, init.angle=90,col =  brewer.pal(nrow(FEAST_output), "Paired"),#最多可用12种颜色梯度
      border="black",main =colnames(FEAST_output)[i] )
}

dev.off()



#简单出图（所有repeat求平均后出1张图）
head(FEAST_output)
asx = as.data.frame(rowMeans(FEAST_output))

asx  = as.matrix(asx)
asx_norm = t(t(asx)/colSums(asx)) #* 100 # normalization to total 100
head(asx_norm)


pdf(file = "FEAST_mean.pdf",width = 6,height = 6)
labs <- paste0(row.names(asx_norm)," (", round(asx_norm[,1]/sum(asx_norm[,1])*100,2), "%)")

pie(asx_norm[,1],labels=labs, init.angle=90,col =  brewer.pal(nrow(FEAST_output), "Paired"),#最多12个颜色梯度
    border="black",main = "mean of source tracker")
dev.off()
names <- c(as.character(envs[train.ix]), "unknown")
pro_data <- data.frame(names)
for (i in 1:70) {
  pro_data[paste0("F", i)] <- c(Proportions_est[[i]])
}
write.csv(pro_data,"pro_data.csv")
data <- melt(pro_data,id.vars='names')
# 过滤比例小于 1% 的来源
data <- data[data$value>0.01,]
# 保留一位小数
data$value <- round(data$value*100, digits = 1)
library(plyr)
library(RColorBrewer)
data1 <- aggregate(value ~ names + variable, data = data, sum)
(p2 <- ggplot(data=data1,aes(variable,value,fill=names))+
    geom_bar(stat="identity", position="fill",color="black", width=0.8,size=0.25)+
    scale_fill_manual(values=brewer.pal(12,"Set3")[c(1:5)])+
    labs(x = '', y = '',fill='Source(%)') +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
      legend.position = "right"
    ))

ggsave("FEAST_Multi_sinks.pdf",width = 12,height = 4)
################两组均值#####################
merge_tax <- pro_data
merge_tax <- aggregate(. ~ names,data=merge_tax,FUN="sum")

rownames(merge_tax) <- merge_tax$names
merge_tax <- merge_tax[,-1]
mat_t = t(merge_tax)
sampFile <- read.table("mapping.txt",header = T,row.names = 1)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,c(-1,-3,-4,-5,-6,-7)]

# 按组求均值，转置，再添加列名
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

# 保存变量备份，并输出至文件
mean_sort=as.data.frame(mat_mean_final)

# 数据转换长表格并绘图
mean_sort$source = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("source")))
data_all$value <- as.numeric(data_all$value)
# 设置分类学顺序，默认字母，可选丰度或手动
# data_all$tax  = factor(data_all$tax, levels=rownames(mean_sort))   

p = ggplot(data_all, aes(x=variable, y = value, fill = source )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+ theme_classic()
p
ggsave("FEAST_Multi_sinks_group_mean.pdf",width = 8,height = 4)
##################################饼图
merge_tax <- pro_data
merge_tax <- aggregate(. ~ names,data=merge_tax,FUN="sum")

rownames(merge_tax) <- merge_tax$names
merge_tax <- merge_tax[,-1]
mat_t = t(merge_tax)
sampFile <- read.table("mapping.txt",header = T,row.names = 1)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,c(-1,-3,-4,-5,-6,-7,-8)]

# 按组求均值，转置，再添加列名
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$env
colnames(mat_mean_final) = geno

# 保存变量备份，并输出至文件
mean_sort=as.data.frame(mat_mean_final)
colnames(mean_sort) <- "Sink"
mean_sort$source <- rownames(mean_sort)
mean_sort$Sink <- as.numeric(mean_sort$Sink)
color <- rev(brewer.pal(nrow(mean_sort),"Oranges"))
labs <- paste0(mean_sort$source,"\n(",round(mean_sort$Sink*100,2),"%)")
pie(mean_sort$Sink,labels = labs,init.angle = 90,col = color,border = "black")


######################################################################
site_data <- merge_tax
site_data <- as.data.frame(t(site_data))
site_data$time3 <- sampFile$time3
site_data$site2 <- sampFile$site2

site_data_all <- melt(site_data, id=c("time3","site2"), variable.name="names", value.name = "value")


##########绘制冲击图
library("reshape2", quietly=T, warn.conflicts=F)
library(ggalluvial)

main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7))

site_time_data <- aggregate(value ~ names + time3+site2, data = site_data_all, sum)

ggplot(data = site_time_data, aes(x = time3, y = value, alluvium = names)) +
  geom_alluvium(aes(fill = names,colour = names), alpha = .75) +
  main_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Source Tracker")+ylab("Percentage(%)")+facet_wrap(~site2)


