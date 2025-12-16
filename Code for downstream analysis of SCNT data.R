setwd("##/##/##")

file_paths <- list.files(pattern = "\\.counts")

# 打印文件路径
print(file_paths)


dfs <- list()

# 遍历每个文件
for (file_path in file_paths) {
  # 打印当前处理的文件路径
  print(paste("Processing file:", file_path))
  
  # 读取文件时只选择 Geneid 和最后一列（计数列）
  df <- read.table(file_path, header = TRUE, sep = "\t", comment.char = "#")
  
  # 获取最后一列的列名（计数列）
  count_col <- colnames(df)[ncol(df)]
  
  # 只保留 Geneid 和计数列
  df <- df[, c("Geneid", count_col)]
  
  # 重命名列
  colnames(df) <- c("Geneid", gsub("\\.counts", "", basename(file_path)))  # 文件名去除后缀作为列名
  
  # 将处理过的数据框加入列表
  dfs[[file_path]] <- df
}

# 合并所有数据框
final_df <- dfs[[1]]
for (i in 2:length(dfs)) {
  final_df <- merge(final_df, dfs[[i]], by = "Geneid", all = TRUE)
}

# 查看合并后的结果
head(final_df)

# 将最终结果保存到文件
write.csv(final_df, "###_Counts.csv", row.names = FALSE)


two_cell_count<-read.csv("###_Counts.csv")

colnames(two_cell_count)<-c("Geneid","CON_1","CON_2","ICSI_1","ICSI_2",
                            "OE_1","OE_2")
rownames(two_cell_count)<-two_cell_count[,1]
two_cell_count<-two_cell_count[,-1]


keep<-rowSums(two_cell_count>0) >= floor(0.5*ncol(two_cell_count))
table(keep)
final_df_1<- two_cell_count[keep,]

zero_var_cols <- apply(final_df_1, 1, function(x) var(x) == 0)
final_df_1_filtered <- final_df_1[!zero_var_cols,]
pca = prcomp(t(final_df_1_filtered), center = TRUE,scale. = TRUE)
pca.var <- pca$sdev^2  ## sdev是标准偏差，十个样本，就有十个标准偏差，平方是避免负数的干扰
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  ##求每个样本的variation
barplot(pca.var.per, main="Screen Plot", xlab="Principal Component", ylab="Percent Variation")  ##用柱状图可视化
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Z=pca$x[,3])

# 添加一列“group"将样品名称加入
pca.data$group <- c("CON","CON","ICSI","ICSI",
                    "OE","OE")
# 利用pdf()函数对作图结果进行保存
library(ggforce)

p<-ggplot(data=pca.data,aes(x=X,y=Y,color=group))+
  geom_point()+
  geom_mark_ellipse(aes(fill=group),alpha=0.1)+
  theme_bw()+
  coord_cartesian(clip = "off")+
  theme(plot.margin = margin(10,10,10,50),
        legend.background = element_blank())+
  xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
  ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))
print(p)
p + scale_x_continuous(limits = c(-200,200)) + scale_y_continuous(limits = c(-400,400))

library(DESeq2)

final_df_filter_1<-round(final_df_1)####round是四舍五入取整
condition <- factor(c(rep("CON",2),rep("ICSI",2),rep("Nr5a2_OE",2)))
colData <- data.frame(row.names=colnames(final_df_filter_1), condition)
dds_1 <- DESeqDataSetFromMatrix(final_df_filter_1, DataFrame(condition), design= ~ condition )
dds_1 <- DESeq(dds_1) 
#res <- results(dds)
res_con_ICSI <- results(dds_1, contrast=c("condition","CON","ICSI"))
res_con_OE <- results(dds_1, contrast=c("condition","Nr5a2_OE","CON"))
summary(res_con_ICSI)
table(res_con_ICSI$pvalue<0.05) #取P值小于0.05的结果
res_con_ICSI <- res_con_ICSI[order(res_con_ICSI$pvalue),]
#diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1.5 | log2FoldChange < -1.5))
#diff_gene_deseq2 <- row.names(diff_gene_deseq2)
FC<-2
resdata <-  merge(as.data.frame(res_con_ICSI),as.data.frame(counts(dds_1,normalize=TRUE)),by="row.names",sort=FALSE)
row.names(resdata)<-resdata[,1]
resdata$regulated <- "normal"
loc_up<-intersect(which(resdata$log2FoldChange>log2(FC)),which(resdata$pvalue<0.05))
loc_down<-intersect(which(resdata$log2FoldChange<(-log2(FC))),which(resdata$pvalue<0.05))
resdata$regulated[loc_up]<-"up"
resdata$regulated[loc_down]<-"down"
table(resdata$regulated)

summary(res_con_OE)
table(res_con_OE$pvalue<0.05) #取P值小于0.05的结果
res_con_OE <- res_con_OE[order(res_con_OE$pvalue),]
#diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1.5 | log2FoldChange < -1.5))
#diff_gene_deseq2 <- row.names(diff_gene_deseq2)
FC<-2
resdata_con_OE <-  merge(as.data.frame(res_con_OE),as.data.frame(counts(dds_1,normalize=TRUE)),by="row.names",sort=FALSE)
row.names(resdata_con_OE)<-resdata_con_OE[,1]
resdata_con_OE$regulated <- "normal"
loc_up<-intersect(which(resdata_con_OE$log2FoldChange>log2(FC)),which(resdata_con_OE$pvalue<0.05))
loc_down<-intersect(which(resdata_con_OE$log2FoldChange<(-log2(FC))),which(resdata_con_OE$pvalue<0.05))
resdata_con_OE$regulated[loc_up]<-"up"
resdata_con_OE$regulated[loc_down]<-"down"
table(resdata_con_OE$regulated)


write.csv(resdata,file= "###_con_vs_ICSI_DESeq2.csv",row.names = F)
write.csv(resdata_con_OE,file= "###_con_vs_OE_DESeq2.csv",row.names = F)

library(tidyverse)
resdata[,1]<-str_sub(resdata[,1],1,18)##保留前18位，去除小数点后
resdata_con_OE[,1]<-str_sub(resdata_con_OE[,1],1,18)##保留前18位，去除小数点后
colnames(resdata)[1]<-"gene_id"
colnames(resdata_con_OE)[1]<-"gene_id"
resdata_1 <- left_join(gff, resdata, by = c("gene_id" = "gene_id"))
resdata_1<-na.omit(resdata_1)
resdata_1<-resdata_1[,-2]
resdata_1<-resdata_1[!duplicated(resdata_1[[1]]), ]
table(resdata_1$regulated)



resdata_con_OE_1 <- left_join(gff, resdata_con_OE, by = c("gene_id" = "gene_id"))
resdata_con_OE_1<-na.omit(resdata_con_OE_1)
resdata_con_OE_1<-resdata_con_OE_1[,-2]
resdata_con_OE_1<-resdata_con_OE_1[!duplicated(resdata_con_OE_1[[1]]), ]
table(resdata_con_OE_1$regulated)




####绘制两组差异火山图
data$Group <- "no_sig" 
loc_down_nosig<-intersect(which(data$logFC_CON_ICSI <= (-1)),which(data$p1 < 0.05))
loc_up_nosig<-intersect(which(data$logFC_CON_ICSI >= 1),which(data$p1<0.05))
data$Group[loc_down_nosig]<-"down_nosig"
data$Group[loc_up_nosig]<-"up_nosig"
loc_down_up<-intersect(intersect(which(data$logFC_Nr5a2_con >= 1),which(data$logFC_CON_ICSI <= (-1))),which(data$p2<0.05))
loc_up_up<-intersect(intersect(which(data$logFC_Nr5a2_con >= 1),which(data$logFC_CON_ICSI >= 1)),which(data$p2<0.05))
loc_down_down<-intersect(intersect(which(data$logFC_Nr5a2_con <= (-1)),which(data$logFC_CON_ICSI <= (-1))),which(data$p2<0.05))
loc_up_down<-intersect(intersect(which(data$logFC_Nr5a2_con <= (-1)),which(data$logFC_CON_ICSI >= 1)),which(data$p2<0.05))
data$Group[loc_down_up]<-"down_up"
data$Group[loc_up_up]<-"up_up"
data$Group[loc_down_down]<-"down_down"
data$Group[loc_up_down]<-"up_down"
table(data$Group)
data$group2 <- "nosig"
loc_down_up_1<-intersect(intersect(intersect(which(data$logFC_Nr5a2_con >= 1),which(data$logFC_CON_ICSI <= (-1))),
                                   which(data$p1 < 0.05)),which(data$p2 < 0.05))
data$group2[loc_down_up_1]<-"down_up_1"
data_group2<-data[data$group2 == "down_up_1",]
write.csv(data_group2,"down_up_p0.05.csv")
color <- c(down_up_1 = "red",nosig = "gray")
P<-ggplot(data,aes(x=logFC_CON_ICSI,y=logFC_Nr5a2_con,color="gray",alpha = 1))+
  geom_point(size=0.5)+
  scale_color_manual(values = color)+
  labs(x="logFC_CON_ICSI",y="logFC_Nr5a2_con ")+
  theme_classic(
    base_line_size = 0.5 ,# 坐标轴的粗细
    
  )
data_group2<-data[data$Group != "no_sig",]
color <- c(down_up = "#E19279",down_down = "gray",up_down = "#364E7D",up_up = "gray")
P + 
  # 设置 'gray' 颜色的点为底层
  # 使得 'gray' 颜色的点有较低透明度
  geom_point(data = data_group2[data_group2$Group != "nosig", ],
             aes(x = logFC_CON_ICSI, y = logFC_Nr5a2_con, colour = Group), alpha = 1,size=0.5) +
  scale_color_manual(values = color) + 
  scale_size_manual(values = 0.1) + # 设置每个组的点大小
  geom_hline(yintercept = 1) + 
  geom_hline(yintercept = -1) + 
  geom_vline(xintercept = -1) + 
  geom_vline(xintercept = 1) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.8, linetype = "solid")) 


#####绘制down_up的阶段特异性热图
library(dplyr)
down_up_stage_fpkm<-inner_join(down_up,stage_fpkm,by=c("name"="gene"))
down_up_stage_fpkm<-down_up_stage_fpkm[,c(2,9:12)]
rownames(down_up_stage_fpkm)<-down_up_stage_fpkm[,1]

down_up_stage_fpkm<-down_up_stage_fpkm[,-1]
down_up_stage_fpkm = down_up_stage_fpkm[apply(down_up_stage_fpkm, 1, function(x) sd(x)!=0),] 
exp <- apply(down_up_stage_fpkm, 1, scale)
rownames(exp) <- colnames(down_up_stage_fpkm)
exp <- t(exp)

library(ComplexHeatmap)
Heatmap(exp, # 表达矩阵
        
        cluster_rows = TRUE, 
        
        cluster_columns = FALSE,
        
        show_row_names = FALSE,
        
        show_column_names = TRUE,
        
        column_title = "Expression Matrix",
        row_title = "Genes")


library(ComplexHeatmap)
library(circlize)  # 提供 colorRamp2

# 假设表达值在 -2 到 2 之间
col_fun <- colorRamp2(c(-2, 0, 2), c("#364E7D", "white","#E19279" ))

Heatmap(exp,
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_title = "Expression Matrix",
        row_title = "Genes")

####绘制cluster1在各阶段的桑基图
library(dplyr)
library(tidyr)
install.packages("networkD3")
library(networkD3)
state_change <- function(current, previous){
  diff <- current - previous
  
  if (diff >= 1) {
    return("Up")
  } else if (diff <= -1) {
    return("Down")
  } else {
    return("Stable")
  }
}
stages <- c("MII_oocyte", "zygote", "early_2cell", "X2cell", "X4cell", "X8cell", "ICM")


# 生成状态矩阵
df_state <- data.frame(gene = rownames(down_up_stage_fpkm))

df_state$MII_oocyte <- "Stable"   # 基线

for (i in 2:length(stages)) {
  prev <- down_up_stage_fpkm[[ stages[i-1] ]]
  curr <- down_up_stage_fpkm[[ stages[i] ]]
  
  df_state[[ stages[i] ]] <- mapply(state_change, curr, prev)
}
library(tidyr)
library(dplyr)

df_long <- df_state %>%
  pivot_longer(-gene,
               names_to = "stage",
               values_to = "state")
df_sankey <- df_state %>%
  pivot_longer(cols = MII_oocyte:ICM, 
               names_to = "stage", 
               values_to = "state")
library(ggplot2)
library(ggalluvial)
library(dplyr)

# 阶段顺序
# 确保每列都是字符
df_state[] <- lapply(df_state, as.character)

# 转为长格式
df_long <- df_state %>%
  pivot_longer(
    cols = -gene,
    names_to = "stage",
    values_to = "state"
  )

# 设置阶段顺序
df_long$stage <- factor(df_long$stage, levels = c(
  "MII_oocyte", "zygote", "early_2cell",
  "X2cell", "X4cell", "X8cell", "ICM"
))

ggplot(df_long,
       aes(x = stage, stratum = state, alluvium = gene, fill = state)) +
  geom_flow(stat = "flow", alpha = 0.7) +
  geom_stratum() +
  scale_x_discrete(limits = c(
    "MII_oocyte", "zygote", "early_2cell",
    "X2cell", "X4cell", "X8cell", "ICM"
  )) +
  scale_fill_manual(values = c("#364E7D", "grey", "#E19279")) +
  labs(title = "Gene Expression State Transitions across Stages",
       x = "Developmental Stage",
       y = "Number of genes") +
  theme_minimal()


#####火山图绘制
library(ggrepel)
p<-ggplot(
  # 数据、映射、颜色
  resdata, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = regulated), size=0.5) +
  scale_color_manual(values = c("#19927D","grey", "#364E7D")) +
  scale_x_continuous(limits = c(-12, 12)) +  # 横坐标范围设置，超出这个范围的点将被裁剪
  scale_y_continuous(limits = c(0, 40)) + 
  # 坐标轴
  labs(x="log2FoldChange",
       y="-log10(pvalue)")+
  # 图例
  theme(legend.position = "bottom")
p
p + theme_bw()+##去除背景色 
  theme(panel.grid=element_blank())+##去除背景格
  theme(panel.border = element_blank(), axis.line = element_line())###去除边框

p<-ggplot(
  # 数据、映射、颜色
  resdata_2, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = regulated), size=0.5) +
  scale_color_manual(values = c("#364E7D","grey", "#E19279")) +
  scale_x_continuous(limits = c(-10, 10)) +  # 横坐标范围设置，超出这个范围的点将被裁剪
  scale_y_continuous(limits = c(0, 15)) + 
  # 坐标轴
  labs(x="log2FoldChange",
       y="-log10(pvalue)")+
  # 图例
  theme(legend.position = "bottom")

p
p + theme_bw()+##去除背景色 
  theme(panel.grid=element_blank())+##去除背景格
  theme(panel.border = element_blank(), axis.line = element_line())###去除边框


####绘制启动子箱式图
rm(list = ls())  ## 魔幻操作，一键清空~
###加载包
library(tidyverse)
library(stringr)
library(reshape2)
# 输入, sample A 和 B 的 output_denisty 路径
data_path_1 <- "###/###_promoter.density"
data_path_2 <- "###/###other_promoter.density"


# 整合A B的信号值
data_1 <- read_delim(data_path_1,
                     delim = "\t",
                     skip = 1,
                     col_names = F) %>%
  na.omit()

data_2 <- read_delim(data_path_2,
                     delim = "\t",
                     skip = 1,
                     col_names = F) %>%
  na.omit()


############# 样本A
# 保留区间id和信号值列
data_agg_1 <- data_1[,-c(1:3,5,6)]

data_agg_2 <- data_2[,-c(1:3,5,6)]

#####绘制启动子箱式
# 计算每个基因的信号和
data_agg_1_sum <- rowSums(data_agg_1[,-1]) %>%
  as.data.frame()

data_agg_2_sum <- rowSums(data_agg_2[,-1]) %>%
  as.data.frame()


data_agg_1_sum<-scale(data_agg_1_sum,center = F,scale = T)%>%
  as.data.frame()
data_agg_2_sum<-scale(data_agg_2_sum,center = F,scale = T)%>%
  as.data.frame()
# 添加样本标签，用于画图
data_agg_1_sum$type <- "cluster_1"

data_agg_2_sum$type <- "cluster_other"

# 添加列名
colnames(data_agg_1_sum)<-c("number","type")

colnames(data_agg_2_sum)<-c("number","type")


data_agg_sum<-rbind(data_agg_1_sum,data_agg_2_sum)


######ggplot2绘制箱式图
library(ggplot2)
library(ggpubr)

library(patchwork)
data_agg_sum$type <- factor(data_agg_sum$type,levels = c("cluster_1", "cluster_other"))
p<-ggplot(data_agg_sum, aes(x=type, y=log2(number))) + 
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条先不会被箱体覆盖
  geom_boxplot(aes(fill=type), #绘制箱线图函数
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        legend.position="none",plot.title = element_text(size=14))+#图例位置
  labs(title="Plot of Nr5a2 promoter",x="cluster", y = "normalize RPKM")

p + scale_fill_manual(values = c(  "#F0C986","gray"))+
  stat_compare_means(comparisons = list(c("cluster_1", "cluster_other")),
                     method = "t.test")

#####显示点的散点图
#p <- ggplot(data_agg_sum, aes(x = type, y = log2(number))) + 
  #stat_boxplot(geom = "errorbar", width = 0.1, size = 0.8) + # 添加误差线
  #geom_boxplot(aes(fill = type), outlier.colour = NA, size = 0.8) + # 箱线图, 隐藏离群点
  #geom_jitter(aes(color = type), width = 0.2, size = 1.5, alpha = 0.7) + # 散点图，调小点的大小
  #theme(panel.background = element_blank(),
   #     axis.line = element_line(),
    #    legend.position = "none",
     #   plot.title = element_text(size = 14)) +
  #labs(title = "Plot of Nr5a2 promoter", x = "cluster", y = "normalize RPKM") +
  #scale_fill_manual(values = c("#F0C986", "gray")) +
  #scale_color_manual(values = c("#F0C986", "gray")) +
  #stat_compare_means(comparisons = list(c("cluster_1", "cluster_other")),
   #                  method = "t.test")

#p



######绘制远端
data_path_1_distal <- "###/###_distal.density"
data_path_2_distal <- "###/###other_distal.density"


data_1_distal <- read_delim(data_path_1_distal,
                            delim = "\t",
                            skip = 1,
                            col_names = F) %>%
  na.omit()

data_2_distal <- read_delim(data_path_2_distal,
                            delim = "\t",
                            skip = 1,
                            col_names = F) %>%
  na.omit()



# 保留区间id和信号值列
data_agg_1_distal <- data_1_distal[,-c(1:3,5,6)]


data_agg_2_distal <- data_2_distal[,-c(1:3,5,6)]





#####绘制启动子箱式
# 计算每个区间的信号平均值
data_agg_1_distal_mean <- colMeans(data_agg_1_distal[,-1]) %>%
  as.data.frame()

data_agg_2_distal_mean <- colMeans(data_agg_2_distal[,-1]) %>%
  as.data.frame()


data_agg_1_distal_mean<-scale(data_agg_1_distal_mean,center = F,scale = T)%>%
  as.data.frame()
data_agg_2_distal_mean<-scale(data_agg_2_distal_mean,center = F,scale = T)%>%
  as.data.frame()


# 添加bin的顺序，用于画图
data_agg_1_distal_mean$bin_number <- c(nrow(data_agg_1_distal_mean):1)

data_agg_2_distal_mean$bin_number <- c(nrow(data_agg_2_distal_mean):1)


#####将矩阵进行升序排列
data_agg_1_distal_mean<-data_agg_1_distal_mean[order(data_agg_1_distal_mean[,2]),]
data_agg_2_distal_mean<-data_agg_2_distal_mean[order(data_agg_2_distal_mean[,2]),]

data_agg_1_distal_mean<-data_agg_1_distal_mean[c(16:120),]
data_agg_2_distal_mean<-data_agg_2_distal_mean[c(16:120),]

data_agg_1_distal_mean$rowsums <- 0
for(i in 1:nrow(data_agg_1_distal_mean)){
  data_agg_1_distal_mean$rowsums[i] <- sum(data_agg_1_distal_mean[1:i,1])
}



data_agg_2_distal_mean$rowsums <- 0
for(i in 1:nrow(data_agg_2_distal_mean)){
  data_agg_2_distal_mean$rowsums[i] <- sum(data_agg_2_distal_mean[1:i,1])
}



# 添加样本标签，用于画图
data_agg_1_distal_mean$type <- "cluster_1"


data_agg_2_distal_mean$type <- "cluster_other"

# 添加列名
colnames(data_agg_1_distal_mean) <- c("signal","bin_number","rowSums_signal" ,"type")

colnames(data_agg_2_distal_mean) <- c("signal", "bin_number","rowSums_signal","type")









# 合并 A and B信号值
data_agg_mean <-  rbind(data_agg_1_distal_mean, data_agg_2_distal_mean)

data_plot <- data_agg_mean

fill_value <- c(  "#F0C986","gray")
ggplot(data = data_plot, aes(x = bin_number,
                             y = log10(rowSums_signal+1),
                             colour = type),linetype=solid) +
  geom_line(size=1.5)+
  # legend color and label
  scale_color_manual(values = fill_value,
  )+
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
  )+#图例位置
  labs(title="Plot of Nr5a2 distal",x="distance to TSS", y = "normalize RPKM")
p + ylim(0,20) 
