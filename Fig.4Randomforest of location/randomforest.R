#  Random Forest Regression

## regression analysis

# Read the experimental design and species classification files
tc_map =read.table("design.txt",header = T, row.names = 1)

otu_table =read.table("sum_g.txt",header = T, row.names = 1)
# Screening training sets
sub_map = tc_map[tc_map$group %in% c("G1"),] # ,"IR24"
# Screening OTUs
idx = rownames(sub_map) %in% colnames(otu_table)
sub_map = sub_map[idx,]
sub_otu = otu_table[, rownames(sub_map)]   

## Random forest regression
library(randomForest)
set.seed(315)
#sub_map$site2 <- as.factor(sub_map$site2)
rf = randomForest(t(sub_otu), sub_map$site2, importance=TRUE, proximity=TRUE, ntree = 1000)
print(rf)

## Cross validation select features
set.seed(315) # Random data ensures repeatable results
# Random Forest Cross Validation
result = rfcv(t(sub_otu), sub_map$site2, cv.fold=10)
# Looking at the error rate table, 23 is the best model with the lowest error rate
result$error.cv
# Draw validation results
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

# Export feature importance
imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp)
write.table(imp,file = "importance_class.txt",quote = F,sep = '\t', row.names = T, col.names = T)
# Simple visualization
varImpPlot(rf, main = "Top 23 - Feature OTU importance",n.var = 25, bg = par("bg"),
           color = par("fg"), gcolor = par("fg"), lcolor = "gray" )


# Read all feature contributions
imp = read.table("importance_class.txt", header=T, row.names= 1, sep="\t") 
# Top 23 was the best choice
imp = head(imp, n=23)
# Sort the x-axis in reverse so that the histogram is drawn from top to bottom
imp=imp[order(1:23,decreasing = T),]

# Draw a histogram of species type and importance
library(ggplot2)
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
imp$genus <- rownames(imp)
imp1 <- imp[order(imp$IncNodePurity,decreasing = F),]
imp1$genus=factor(imp1$genus,levels = imp1$genus)
p=ggplot(data = imp1, mapping = aes(x=genus,y=IncNodePurity,fill=genus)) + 
  geom_bar(stat="identity")+coord_flip()+main_theme
p
ggsave(paste("rf_imp_feature",".pdf", sep=""), p, width = 6, height =4)


# Load heat map drawing package
library(pheatmap)

# Data filtering and 23 features display
sub_abu = sub_otu[rownames(imp1),]

rownames(sub_abu)=imp1[rownames(sub_abu),"genus"]

pheatmap(sub_abu, scale = "row")
# save results
pheatmap(sub_abu, scale = "row", filename = "heatmap_samples.pdf", width = 5, height = 5)


# Combined and averaged by time
sampFile = as.data.frame(sub_map$site2,row.names = row.names(sub_map))
colnames(sampFile)[1] = "group"
mat_t = t(sub_abu)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
otu_norm_group = do.call(rbind, mat_mean)[-1,]
colnames(otu_norm_group) = mat_mean$group
pheatmap(otu_norm_group,scale="row",cluster_cols = F, cluster_rows = T)
pheatmap(otu_norm_group, scale="row",cluster_cols = F, cluster_rows = T, filename = "heatmap_groups.pdf", width = 5, height = 5)



## Find the maximum value of each group
bak=otu_norm_group

otu_norm_group = otu_norm_group[as.character(imp1$genus),] 

for (i in 1:length(rownames(otu_norm_group))) {
#  i=1
  x=as.data.frame(sort(otu_norm_group[i,],decreasing = T))
  imp1[i,"order"]=rownames(x)[1]
}
library(dplyr)
imp1$order2 =  as.numeric(gsub("012","",imp1$order,perl=TRUE) )
taxonomy = arrange(imp1, desc(order2), genus)

otu_norm_group1 = otu_norm_group[match(taxonomy$genus,rownames(otu_norm_group)),] # 按初始排序

pheatmap(otu_norm_group1,scale="row",cluster_cols = F, cluster_rows = F)
#0,1,2 correspond to CC, PC, SC respectively
pheatmap(otu_norm_group1,scale="row",cluster_cols = F, cluster_rows = F,filename ="pheatmap_order_all.pdf",width=8, height=4)


