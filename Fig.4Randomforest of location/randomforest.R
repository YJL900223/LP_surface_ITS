# 随机森林回归 Random Forest Regression

## 回归分析

# 读取实验设计、和物种分类文件
tc_map =read.table("design.txt",header = T, row.names = 1)
# 物种分类文件，由qiime summarize_taxa.py生成，详见扩增子分析流程系列
otu_table =read.table("sum_g.txt",header = T, row.names = 1)
# 筛选品种作为训练集
sub_map = tc_map[tc_map$group %in% c("G1"),] # ,"IR24"
# 筛选OTU
idx = rownames(sub_map) %in% colnames(otu_table)
sub_map = sub_map[idx,]
sub_otu = otu_table[, rownames(sub_map)]   

## 随机森林回归
library(randomForest)
set.seed(315)
#sub_map$site2 <- as.factor(sub_map$site2)
rf = randomForest(t(sub_otu), sub_map$site2, importance=TRUE, proximity=TRUE, ntree = 1000)
print(rf)

## 交叉验证选择Features
set.seed(315) # 随机数据保证结果可重复，必须
# rfcv是随机森林交叉验证函数：Random Forest Cross Validation
result = rfcv(t(sub_otu), sub_map$site2, cv.fold=10)
# 查看错误率表，23时错误率最低，为最佳模型
result$error.cv
# 绘制验证结果 
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

##########################导出图片

# 导出feature重要性
imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp)
write.table(imp,file = "importance_class.txt",quote = F,sep = '\t', row.names = T, col.names = T)
# 简单可视化
varImpPlot(rf, main = "Top 23 - Feature OTU importance",n.var = 25, bg = par("bg"),
           color = par("fg"), gcolor = par("fg"), lcolor = "gray" )


## ggplot2美华feature贡献度柱状图

# 软件内部的varImpPlot可以快速可视化贡献度，简单全面，但发表还是要美美哒，美是需要代码的，就是花时间
# 基本思路同绘制Top 23 feature柱状图，按门着色，简化纲水平名字

# 读取所有feature贡献度
imp = read.table("importance_class.txt", header=T, row.names= 1, sep="\t") 
# 分析选择top23分组效果最好
imp = head(imp, n=23)
# 反向排序X轴，让柱状图从上往下画
imp=imp[order(1:23,decreasing = T),]

# 绘制物种类型种重要性柱状图
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



# 绘制时间序列热图

# 加载热图绘制包
library(pheatmap)

# 数据筛选23个feature展示
sub_abu = sub_otu[rownames(imp1),]

# 简化名字
rownames(sub_abu)=imp1[rownames(sub_abu),"genus"]

# 直接自动聚类出图
pheatmap(sub_abu, scale = "row")
# 保存结果
pheatmap(sub_abu, scale = "row", filename = "heatmap_samples.pdf", width = 5, height = 5)


# 按时间为组合并均值
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



## 求每组最大值
bak=otu_norm_group

otu_norm_group = otu_norm_group[as.character(imp1$genus),] # 按初始排序

for (i in 1:length(rownames(otu_norm_group))) {
#  i=1
  x=as.data.frame(sort(otu_norm_group[i,],decreasing = T))
  imp1[i,"order"]=rownames(x)[1]
}
library(dplyr)
imp1$order2 =  as.numeric(gsub("012","",imp1$order,perl=TRUE) )
taxonomy = arrange(imp1, desc(order2), genus)

otu_norm_group1 = otu_norm_group[match(taxonomy$genus,rownames(otu_norm_group)),] # 按初始排序

# 按初始排序
pheatmap(otu_norm_group1,scale="row",cluster_cols = F, cluster_rows = F)

pheatmap(otu_norm_group1,scale="row",cluster_cols = F, cluster_rows = F,filename ="pheatmap_order_all.pdf",width=8, height=4)
###########0,1,2分别对应CC,PC,SC

