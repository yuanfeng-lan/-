####cibersort####
setwd("cibersort")   
install.packages('e1071')
install.packages('parallel')
#install.packages("BiocManager")
BiocManager::install("preprocessCore", version = "3.15")
library(e1071)
library(parallel)
library(preprocessCore)
source("CIBERSORT.R")   
library(tidyverse)
sig_matrix <- "LM22.txt"   
mixture_file = 'normalize.txt'   #肿瘤患者表达谱
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
#QN T 为芯片数据 F为测序数据
save(res_cibersort,file = "res_cibersort.Rdata")   #保存中间文件

load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22]   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞
#可视化（阿琛老师）
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框写
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板


library(ggpubr)
library(ggplot2)
# boxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Cell composition",
  main = "TME Cell composition"
) +
  theme_base() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ))




####方法2分析####
remove(list = ls()) #一键清空
#加载包
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

source("CIBERSORT.R")

# 设置分析依赖的基础表达文件
# 每类免疫细胞的标志性基因及其表达
# 基因名字为Gene symbol
LM22.file <-"LM22.txt" 
# 1. Cibersort
exp <- read.table("normalize.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- 2^exp
group <- read.table("re1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
TCGA_TME.results <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

write.table(exp, file = "exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
TCGA_exp.file <- "exp2.txt"

TCGA_TME.results <- CIBERSORT(LM22.file ,TCGA_exp.file, perm = 1000, QN = F)  
# perm置换次数=1000
# QN如果是芯片设置为T，如果是测序就设置为F

save(TCGA_TME.results,file = "res_cibersort.Rdata") 

## 2. 分组信息

# TCGA的数据还可以从名字获取
# group_list <- ifelse(as.numeric(substring(rownames(TCGA_TME.results),14,15)) < 10,
#                    "Tumor","Normal") %>% 
#  factor(.,levels = c("Normal","Tumor"))

phenotype = read.csv("group.txt",header = T,row.names = 1)
phenotype <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
phenotype <- group
group_list <- phenotype$group1 %>% 
  factor(.,levels = c("HC","HO","SS","NASH","HCC"))

table(group_list) # Normal 43   Tumor 43 

## group_list
## Nontumor    Tumor 
##       43       43
## 3. 绘图
# 3.1 数据粗处理
TME_data <- as.data.frame(res_cibersort[,1:22])
TME_data <- as.data.frame(TCGA_TME.results[,1:22])
TME_data <- TME_data[rownames(phenotype),]
TME_data$group <- group_list
TME_data$sample <- row.names(TME_data)
TME_data <- cbind(TME_data,phenotype)
write.table(TME_data, file = "TME_data.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#可视化（阿琛老师）
TME_data <- TME_data[,1:22]
mycol <- ggplot2::alpha(rainbow(ncol(TME_data)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(TME_data)),
        border = NA, # 柱子无边框写
        names.arg = rep("",nrow(TME_data)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-10, # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(TME_data), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板
# 2.2 融合数据#宽数据转化为长数据
TME_New = melt(TME_data)

## Using group, sample as id variables

colnames(TME_New)=c("Group","Sample","Celltype","Composition")  #设置行名
head(TME_New)

##      Group          Sample      Celltype Composition
## 1    Tumor TCGA.CV.6943.01 B cells naive 0.007651678
## 2    Tumor TCGA.CV.6959.01 B cells naive 0.019549031
## 3 Nontumor TCGA.CV.7438.11 B cells naive 0.025349204
## 4 Nontumor TCGA.CV.7242.11 B cells naive 0.032583659
## 5    Tumor TCGA.CV.7432.01 B cells naive 0.000000000
## 6 Nontumor TCGA.CV.6939.11 B cells naive 0.074282293

# 3.3 按免疫细胞占比中位数排序绘图（可选）
plot_order = TME_New[TME_New$Group=="HC",] %>% 
  group_by(Celltype) %>% 
  summarise(m = median(Composition)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltype)

## `summarise()` ungrouping output (override with `.groups` argument)

TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)

# 3.3 出图
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }

box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "TME Cell composition")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369", '#ffe119', '#f58231', '#911eb4'))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "anova",
                     hide.ns = T)

box_TME


box_TME;ggsave("ciber2.pdf",box_TME,height=15,width=25,unit="cm")

head(gapminder,n=4)
library(gapminder)
TME_New %>%
  ggplot(aes(x=Group,y=Composition,fill=Group)) +
  geom_boxplot() + geom_jitter(width = 0.1,alpha = 0.2) +
  facet_wrap(~Celltype,ncol = 8,nrow = 3) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+ 
  stat_compare_means(label.y = 0.4,label = "p.signif")

TME_New2 <- TME_New[TME_New$Celltype==c("Neutrophils"),]
TME_New2 %>%
  ggplot(aes(x=Group,y=Composition,fill=Group)) +
  geom_boxplot() + geom_jitter(width = 0.1,alpha = 0.2) +
  facet_wrap(~Celltype) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+ 
  stat_compare_means(label.y = 0.4,label = "p.signif")
