####GEO数据库数据处理####
library(tidyverse)
library(GEOquery)
gset = getGEO('GSE50760', destdir=".", AnnotGPL = F, getGPL = F)
gset[[1]]
#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
write.table(pdata, file = "phe.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

pdata <- read.table("phe_her2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp <- read.table("GSE50760_norm_counts_FPKM_GRCh38.p13_NCBI.tsv",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

index <- read.csv("Human.GRCh38.p13.annot.tsv",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
table(index$GeneType)


exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- read.table("KIRC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- read.table("BRCA_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- read.table("GBM_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- read.table("CHOL_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- read.table("KIRP_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


exp <- read.table("normalize.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gene <- read.table("PDZ list.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- read.table("PDZ list2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gene <- gene$gene
exp2 <- exp[gene,]


exp <- cbind(Name = rownames(exp), DESCRIPTION = "na", exp)

write.table(exp2, file = "exp2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
exp2 <- read.table("exp2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

write.table(file = "KIRC_fpkm.txt", exp, sep = "\t", col.names = T, row.names = F, quote = F)
#热图展示结果
library(pheatmap)


group$group1 <- as.factor(group$group1)
group$group1 = factor(group$group1,
                       levels = c("HC","HO",'SS','NASH','HCC'))

group_list=factor(ifelse(substr(colnames(exp),14,16) == "01A","Tumor","Normal"),levels = c("Normal","Tumor"))
group_list=ifelse(substr(colnames(exp),14,16) == "01A","Tumor","Normal")

sample <- colnames(exp)
group <- cbind(sample,group_list)
group <- as.data.frame(group)
rownames(group) <- group$sample
group <- group[,-1,drop=F]
write.table(group, file = "group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


exp2 <- exp2[,rownames(group)]

#热图色板
install.packages("RColorBrewer")
library(RColorBrewer)

install.packages("gplots")
library(gplots)




library(limma) 
exp3=normalizeBetweenArrays(exp2)
boxplot(exp3,outline=FALSE, notch=T,col=group_list, las=2)

####差异分析####
group_list
design <- model.matrix(~0+factor(group$group1))
row.names(design) <- rownames(group)

colnames(design)=c('HC1','NASH2')
contrast.matrix<-makeContrasts("HC1-NASH2",levels=design)

##step1
fit <- lmFit(exp2,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "1-2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

exp2 <- exp2[,rownames(gs_exp)]

group <- gs_exp[,c(1,3,8,9,10,11)]
  
barplot(1:25,col=colorRampPalette(brewer.pal(11, "PiYG"))(25))
coul <- colorRampPalette(brewer.pal(11, "PiYG"))(25)
coul <- colorRampPalette(brewer.pal(9, "OrRd"))(50)

coul <- colorRampPalette(brewer.pal(8, "Dark2"))(50)

exp3 <- exp3^2

pheatmap(exp2,
         annotation_col=group,
         scale = "row",
         show_rownames = T,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=3)
dev.off()
bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))

pheatmap(exp2,
         annotation_col=group,
         scale = "row",
         show_rownames = T,
         show_colnames =F,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=3,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         
         legend_breaks=seq(-5,5,2),
         
         breaks=bk)
dev.off()


GSE48452 <- read.table("GSE48452.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE48452 <- GSE48452[gene,]
GSE48452 <- GSE48452[,rownames(GSE48452group)]
GSE48452group <- read.table("GSE48452-group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE48452group$group = factor(GSE48452group$group,
                             levels = c("Control","Healthy_obese",'Steatosis','Nash'))

pheatmap(GSE48452,
         annotation_col=GSE48452group,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)


GSE89632 <- read.table("GSE89632.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE89632 <- GSE164760[gene,]

GSE89632group <- read.table("GSE89632-group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE89632 <- GSE89632[,rownames(GSE164760group)]
GSE89632group$group = factor(GSE89632group$group,
                             levels = c("HC","SS",'NASH'))
pheatmap(GSE89632,
         annotation_col=GSE89632group,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)

GSE164760 <- read.table("GSE164760.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE164760 <- GSE164760[gene,]

GSE164760group <- read.table("GSE164760-group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE164760 <- GSE164760[,rownames(GSE164760group)]

pheatmap(GSE164760,
         annotation_col=GSE164760group,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)






####配对t检验####
library(reshape2)
library(ggpubr)

exp2 <- read.table("exp2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group3.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


gene2 <- read.table("gene2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene2 <- read.table("Lgene2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene2 <- read.table("gene3.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gene2 <- gene2$gene
exp3 <- exp2[gene2,]
exp3 <- as.data.frame(t(exp3))
exp4 <- exp3[,1:19]
exp4$group <- group$group
exp4 <- melt(exp4)


write.table(exp4, file = "exp4.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp4, file = "exp5.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


exp4 %>%
  ggplot(aes(x=group,y=value,fill=group)) +
  geom_boxplot() + geom_jitter(width = 0.1,alpha = 0.2) +
  facet_wrap(~variable,ncol = 4,nrow = 5) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+ 
  stat_compare_means()

####火山图####

library(ggplot2)
library(ggrepel)
# 设置工作目录
setwd("E:/R/WorkSpace/baimoc/visualization")


# 整理数据集
# 参数'./resource/dataset.txt'，表示载入E:/R/WorkSpace/baimoc/visualization/resource/dataset_heatmap.txt
dataset <- read.table("tumor-normal.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dataset <- read.table("NASH-HC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 设置pvalue和logFC的阈值
#cut_off_adj.P.Val = 0.05
#cut_off_logFC = 0.5
# 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘Stable’，保存到change列
# 这里的change列用来设置火山图点的颜色
dataset$change = ifelse(dataset$adj.P.Val < 0.05 & abs(dataset$logFC) >= 0, 
                        ifelse(dataset$logFC> 0 ,'Up','Down'),
                        'Stable')

table(dataset$change)


# 绘制火山图
ggplot(
  #设置数据
  dataset, 
  aes(x = logFC, 
      y = -log10(adj.P.Val), 
      colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#66C3A5", "#d2dae2","#FD8D62"))+
  
  # 辅助线
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  
  # 坐标轴
  labs(title="GEO.NASH",x="log2(fold change)",
       y="-log10 (Adjust p-value)")+
  theme_bw()+
  
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )

####泛癌分析数据处理####


library(tidyverse)
setwd("TCGA.KIRC")
setwd("GDS")
setwd("cibersort")

col = read.table(file = 'GDC-PANCAN.htseq_fpkm-uq.tsv', sep = '\t', header = TRUE,nrows=2,row.names=1) 

fpkm1 = read.table(file = 'GDC-PANCAN.htseq_fpkm-uq.tsv', sep = '\t', header = TRUE,nrows=10500,row.names=1,skip = 49999) 
colnames(fpkm1) <- colnames(col)

table(substr(colnames(fpkm1),14,16))
fpkm1 <- fpkm1[,substr(colnames(fpkm1),14,16)%in% c("01A","11A")]
table(substr(colnames(fpkm1),14,16))
rownames(fpkm1) <- substr(rownames(fpkm1),1,15)
fpkm <- fpkm1

Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo0 <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] #只要编码RNA
#取行名交集
comgene <- intersect(rownames(fpkm),rownames(Ginfo0))
fpkm <- fpkm[comgene,]
Ginfo0 <- Ginfo0[comgene,]
fpkm$Gene <- as.character(Ginfo0$genename)   #新增Gene Symbol

fpkm <- fpkm[!duplicated(fpkm$Gene),] 

rownames(fpkm) <- fpkm$Gene   #将行名变为Gene Symbol
fpkm <- fpkm[,-ncol(fpkm)]   #去除最后一列
#保存所有患者的fpkm文件
write.table(fpkm, file = "PANCANCER_fpkm_mRNA_all_6.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

gene <- read.table("PDZ list2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

comgene <- intersect(rownames(fpkm),gene$gene)

exp2 <- fpkm[comgene,]

write.table(exp2, file = "PANCANCER_fpkm_mRNA_all_PDZ_6.txt",sep = "\t",row.names = T,col.names = NA,quote = F)







#保存癌症患者的fpkm文件
tumor <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "01A"]
fpkm_01A <- fpkm[,tumor]
write.table(fpkm_01A, file = "LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存正常样本的fpkm文件
normal <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "11A"]
fpkm_11A <- fpkm[,normal]
write.table(fpkm_11A, file = "LIHC_fpkm_mRNA_11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#整理完毕#


#分批数据合并

exp1 <- read.table("PANCANCER_fpkm_mRNA_all_PDZ_1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp2 <- read.table("PANCANCER_fpkm_mRNA_all_PDZ_2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp3 <- read.table("PANCANCER_fpkm_mRNA_all_PDZ_3.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp4 <- read.table("PANCANCER_fpkm_mRNA_all_PDZ_4.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp5 <- read.table("PANCANCER_fpkm_mRNA_all_PDZ_5.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp6 <- read.table("PANCANCER_fpkm_mRNA_all_PDZ_6.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp7 <- rbind(exp1,exp2)

exp8 <- rbind(exp3,exp4)

exp9 <- rbind(exp5,exp6)
expq <- rbind(exp7,exp8)
expw <- rbind(expq,exp9)
exp0 <- expw

write.table(exp0, file = "PANCANCER_fpkm_mRNA_all_PDZ.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

exp0 <- read.table("PANCANCER_fpkm_mRNA_all_PDZ.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- exp0
group_list=factor(ifelse(substr(colnames(exp),14,16) == "01A","Tumor","Normal"),levels = c("Normal","Tumor"))

sample <- colnames(exp)
group <- cbind(sample,group_list)
group <- as.data.frame(group)
rownames(group) <- group$sample
group <- group[,-1,drop=F]
write.table(group, file = "group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


#热图
library(pheatmap)
exp2 <- exp[,rownames(group)]
bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))

pheatmap(exp2,
         annotation_col=group,
         scale = "row",
         show_rownames = T,
         show_colnames =F,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=3,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         
         legend_breaks=seq(-5,5,2),
         
         breaks=bk)

####批量差异分析####
library(limma)

design <- model.matrix(~0+factor(group$group))
row.names(design) <- rownames(group)

colnames(design)=c('Normal','Tumor')
contrast.matrix<-makeContrasts("Normal-Tumor",levels=design)

##step1
fit <- lmFit(exp2,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "Normal-Tumor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

LogFC2 <- read.table("LogFC总.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
P2 <- read.table("P总.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

write.table(LogFC, file = "LogFC.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(P, file = "P.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
P <- P[,1,drop=F]
LogFC <- LogFC[,1,drop=F]


a <-  read.table(file="顺序.txt",header=T,sep="\t",row.names=1)
i=1
a[1,i]

for (i in 1:21) {
      group0 <- group[group$Progect==a[1,i],]
      exp2 <- exp0[,rownames(group0)]
      
      design <- model.matrix(~0+factor(group0$group))
      row.names(design) <- rownames(group0)
      colnames(design)=c('Normal','Tumor')
      contrast.matrix<-makeContrasts("Tumor-Normal",levels=design)
      
      ##step1
      fit <- lmFit(exp2,design)#data的列名要等于design的行名
      ##step2
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2) 
      ##step3
      tempOutput = topTable(fit2, coef=1, n=Inf)
      nrDEG = na.omit(tempOutput) 
      
      write.table(nrDEG,file=paste0("limma/",a[2,i],".txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
      nrDEG <- nrDEG[rownames(LogFC),]
      LogFC$FC <- nrDEG$logFC
      colnames(LogFC)[i+1] <- a[2,i]
      
      P$p <- nrDEG$adj.P.Val
      colnames(P)[i+1] <- a[2,i]
      
  
}

write.table(design, file = "phe.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#循环测试
group0 <- group[group$Progect==a[1,i],]
exp2 <- exp0[,rownames(group0)]

design <- model.matrix(~0+factor(group0$group))
row.names(design) <- rownames(group0)
colnames(design)=c('Normal','Tumor')
contrast.matrix<-makeContrasts("Normal-Tumor",levels=design)

##step1
fit <- lmFit(exp2,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 


write.table(nrDEG,file=paste0("limma/",a[2,i],".txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

nrDEG <- nrDEG[rownames(LogFC),]
LogFC$FC <- nrDEG$logFC
colnames(LogFC)[i+1] <- a[2,i]

P$p <- nrDEG$adj.P.Val
colnames(P)[i+1] <- a[2,i]

####绘制热图展示结果####

library(pheatmap)
LogFC1 <- LogFC[,-1]
LogFC <- LogFC[,-23]
P <- P[,-23]
P1 <- P[,-1]
pmt <- P1
if (!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}


group2 <- read.table("group4.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

pheatmap(LogFC1,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =T,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6,
         display_numbers = pmt)

pheatmap(LogFC1,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =T,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6
         )
####批量火山图####
dataset <- read.table("tumor-normal.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

dataset$change = ifelse(dataset$adj.P.Val < 0.05 & abs(dataset$logFC) >= 0, 
                        ifelse(dataset$logFC> 0 ,'Up','Down'),
                        'Stable')

table(dataset$change)


# 绘制火山图
ggplot(
  #设置数据
  dataset, 
  aes(x = logFC, 
      y = -log10(adj.P.Val), 
      colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#66C3A5", "#d2dae2","#FD8D62"))+
  
  # 辅助线
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  
  # 坐标轴
  labs(title="GEO.NASH",x="log2(fold change)",
       y="-log10 (Adjust p-value)")+
  theme_bw()+
  
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )

####层次聚类####
LogFC <- read.table("LogFC.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
LogFC <- LogFC[,-1]

data(nutrient,package = "flexclust")
bc<-nutrient
bc <- as.data.frame(t(LogFC))
bc <- t(bc)

row.names(bc)<-tolower(row.names(bc))
bc.scaled<-scale(bc)##标准化数据

d<-dist(bc.scaled)###计算欧氏距离
d<-dist(bc)###计算欧氏距离

fit1<-hclust(d,method = "average")
plot(fit1,hang = -1,cex=.8,main = "title")

heatmap(as.matrix(d))

library(NbClust)
devAskNewPage(ask = T)
bc.scaled <- as.data.frame(bc.scaled)

nc<-NbClust(bc.scaled,distance = "euclidean",min.nc = 2,max.nc = 15,
            method = "average")
table(nc$Best.n[1,])
barplot(table(nc$Best.n[1,]),xlab = "No. of cluster")


library(factoextra)
library(igraph)
fviz_dend(fit1)

fviz_dend(fit1,k=2,rect =T,rect_fill = T)##K为聚类个数，rect_border为区域颜色填充

####循环相关性分析####


exp2 <- exp2[,colnames(group)]
#创建空向量
gene_name1<-c()##也可用vector
gene_name2<-c()
cor_r<-c()
pvalue<-c()

LogFC <- exp2[,1,drop=F]
P <- exp2[,1,drop=F]


for (i in 1:142){
  for (r in 1:10){
    c_r=cor(as.numeric(exp2[i,]),as.numeric(group[r,]),method="pearson",use = 'pairwise.complete.obs')
    p=cor.test(as.numeric(exp2[i,]),as.numeric(group[r,]),method ="pearson")[[3]]
    ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
    LogFC[i,r+1] <- c_r
    P[i,r+1] <- p
    
    
  }
}




group <- as.data.frame(t(group))
group <- group[-1,]

i=1
r=1
#循环测试

c_r=cor(as.numeric(exp2[r,]),as.numeric(group[i,]),method="pearson",use = 'pairwise.complete.obs')
p=cor.test(as.numeric(exp2[r,]),as.numeric(group[i,]),method ="pearson")[[3]]
##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
LogFC[i,r+1] <- c_r
P[i,r+1] <- p

colnames(LogFC)[r+1] <- rownames(group[r,])
colnames(P)[r+1] <- rownames(group[r,])




LogFC <- LogFC[,-1]
P <- P[,-1]
colnames(LogFC) <- rownames(group)
colnames(P) <- rownames(group)

write.table(LogFC, file = "cor_r.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(P, file = "P.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


LogFC <- read.table("MAGI3.PRAD.cor.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
P <- read.table("MAGI3.PRAD.P.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)




LogFC$KIRC <- LogFC2$KIRC
P$KIRC <- P2$KIRC

#热图展示结果
pmt <- P
if (!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}


group2 <- read.table("group4.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
bk <- c(seq(-0.3,0,by=0.01),seq(0.001,0.3,by=0.01))
pheatmap(LogFC,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         cluster_cols =F,
         cluster_rows = T,
         border=FALSE,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=5,
         display_numbers = pmt,
         main = "TCGA.KIRC",
         
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks=bk
         )



pheatmap(LogFC,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =T,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6,
         display_numbers = pmt,
         main = "TCGA.KIRC"
      
)
####NAFLD分析####
exp <- read.table("normalize.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gene <- read.table("PDZ list2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gene <- gene$gene
exp2 <- exp[gene,]

write.table(exp2, file = "exp2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
exp2 <- read.table("exp2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp2 <- exp2[,rownames(group)]

library(pheatmap)
pheatmap(exp2,
         annotation_col=group,
         scale = "row",
         show_rownames = T,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6
)


library(limma)

design <- model.matrix(~0+factor(group$group1))
row.names(design) <- rownames(group)

colnames(design)=c('HC','HCC',"HO","NASH","SS")
contrast.matrix<-makeContrasts("HCC-HC",levels=design)

##step1
fit <- lmFit(exp2,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "HCC-HC.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

nrDEG <- nrDEG[rownames(exp2),]
#LogFC <- exp2[,1,drop=F]
#P <- exp2[,1,drop=F]

LogFC$HCC <- nrDEG$logFC

P$HCC <- nrDEG$adj.P.Val

LogFC <- LogFC[,-1]
P <- P[,-1]

pmt <- P
if (!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}
pheatmap(LogFC,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         color = colorRampPalette(c( "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = F,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6,
         display_numbers = pmt,
         border=T)

####GSVA赋分####

library(tidyverse)
library(GSEABase)
library(GSVA)
library(pheatmap)

exp <- read.table("LUAD_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)



geneSet <- read.csv("PDZ protein.gmt",header = F,sep = "\t")
geneSet <- read.csv("genesets.v2023.2.Hs.gmt",header = F,sep = "\t")


geneSet <- read.csv("genesets.v2023.2.Hs.gmt",sep = "\t",check.names = F,stringsAsFactors = F,header = F,col.names = paste("V", 1:1409, sep = ""))

read.table("test.txt", fill = T, col.names = paste("V", 1:1407, sep = ""))

data <- geneSet

data <- t(data)


geneSet <- geneSet[,-2]


geneSet <- geneSet %>%
  column_to_rownames("V1")%>%t()
a <- geneSet
a <- a[1:nrow(a),,drop=F]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}

save(l,file = "脂肪酸相关基因集.Rdata")



exp <- as.matrix(exp)
gs.exp <- gsva(exp, l, kcdf = "Gaussian", min.sz = 10)


group_list=factor(ifelse(substr(colnames(gs.exp),14,16) == "01A","Tumor","Normal"),levels = c("Normal","Tumor"))
gs_exp <- as.data.frame(gs_exp)
gs_exp$group <- group_list
#大型数据矩阵的拆分和合并
exp1 <- read.table("PANCANCER_fpkm_mRNA_all_1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp2 <- read.table("PANCANCER_fpkm_mRNA_all_2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp3 <- read.table("PANCANCER_fpkm_mRNA_all_3.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp4 <- read.table("PANCANCER_fpkm_mRNA_all_4.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp5 <- read.table("PANCANCER_fpkm_mRNA_all_5.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp6 <- read.table("PANCANCER_fpkm_mRNA_all_6.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp1 <- rbind(exp1,exp2)
exp1 <- rbind(exp3,exp4)
exp1 <- rbind(exp5,exp6)



exp1_1 <- exp1[,1:2000]
exp1_2 <- exp1[,2001:4000]
exp1_3 <- exp1[,4001:6000]
exp1_4 <- exp1[,6001:8000]
exp1_5 <- exp1[,8001:10244]


write.table(exp1_1, file = "exp/exp3_1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp1_2, file = "exp/exp3_2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp1_3, file = "exp/exp3_3.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp1_4, file = "exp/exp3_4.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp1_5, file = "exp/exp3_5.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


exp1 <- read.table("exp/exp1_5.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp2 <- read.table("exp/exp2_5.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp3 <- read.table("exp/exp3_5.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp4 <- rbind(exp1,exp2)

exp5 <- rbind(exp4,exp3)

write.table(exp5, file = "exp/exp0_1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp5, file = "exp/exp0_2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp5, file = "exp/exp0_5.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

exp1 <- read.table("exp/exp0_5.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp1 <- read.table("KIRC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp1 <- as.matrix(exp1)
gs.exp <- gsva(exp1, l, kcdf = "Gaussian", min.sz = 10)

group <- read.table("group0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gs_exp
group <- group[colnames(gs.exp),]
gs_exp <- t(gs.exp)
gs_exp <- gs_exp[rownames(group),,drop=F]
gs_exp <- as.data.frame(gs_exp)
gs_exp$Progect <- group$Progect
gs_exp$group <- group$group


gs_exp <- cbind(gs_exp,group)
gs_exp <- gs_exp[,-12]

write.table(gs.exp, file = "PDZ protein gene set GSVA/gs.exp_5.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(gs_exp, file = "GSVA1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


library(reshape2)
library(ggpubr)
gs_exp <- read.table("PDZ protein gene set GSVA.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gs_exp <- read.table("GSVA1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


exp2 <- gs_exp[gs_exp$grade==c("G1","G2","G3","G4"),]
exp4 <- melt(gs_exp)
exp4$Grade <- factor(exp4$Grade,levels = c("G1","G2","G3","G4"))
exp4$pathologic_T <- factor(exp4$pathologic_T,levels = c("T1","T2","T3","T4"))
exp4$tumor_stage <- factor(exp4$tumor_stage,levels = c("stage i","stage ii","stage iii","stage iv"))

#exp0 <- exp4[exp4$Progect=="TCGA-UCEC",]


table(exp4$Progect)

p1 <- ggboxplot(exp4,  x = "gender", y = "value",
                color = "gender", palette = "Dark2",
                add = "jitter",na.rm = TRUE)+ #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme(axis.title.x=element_blank())+ #删除x轴坐标名称
  labs(title = "TCGA.LIHC")+
  theme(plot.title = element_text(hjust = 0.5))

p1


p1 <- ggboxplot(gs_exp,  x = "group", y = "PDZ protein",
                color = "group", palette = "Dark2",
                add = "jitter",na.rm = TRUE)+ #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+  
  theme(axis.title.x=element_blank())+ #删除x轴坐标名称
  labs(title = "TCGA.LIHC")+
  theme(plot.title = element_text(hjust = 0.5))

p1

####PDZK1和MAGI3分析####


library(GSEABase)
library(GSVA)
library(pheatmap)


exp <- read.table("PANCANCER_fpkm_mRNA_all_PDZ.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- read.table("group2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group3 <- read.table("group3.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- exp[,rownames(group)]

exp2 <- exp["PDZK1",]
exp3 <- exp["MAGI3",,drop=F]
exp2 <- as.data.frame(t(exp2))

exp3 <- as.data.frame(t(exp3))


group2 <- cbind(exp2,group)
group3 <- cbind(exp3,group)



exp <- read.table("KIRC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- read.table("KIRP_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- read.table("KIRP_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- read.table("LUSC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp <- as.data.frame(exp)
gs.exp2 <- gs.exp
exp <- as.matrix(exp)
gs.exp <- gsva(exp, l, kcdf = "Gaussian", min.sz = 10)
write.table(gs.exp, file = "GSEA.exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


exp <- read.table("GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group_PDZK1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
sample <- intersect(rownames(group),colnames(exp))

exp <- exp[,sample]
group <- group[sample,]
exp <- as.data.frame(t(exp))


data <- cbind(exp,group)


library(pheatmap)
gs.exp <- gs.exp[,rownames(group2)]
gs.exp <- gs.exp[,rownames(group3)]

gs.exp0 <- gs.exp0[,rownames(group2.0)]
gs.exp0 <- gs.exp0[,rownames(group3.0)]


pheatmap(gs.exp,
         annotation_col=group3,
         scale = "none",
         show_rownames =T ,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=4,
         main = "TCGA.LUSC"
)

pheatmap(gs.exp0,
         annotation_col=group2.0,
         scale = "row",
         show_rownames =T ,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=6,
         main = "TCGA.LIHC KIRC"
)



write.table(group3, file = "group3.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(gs.exp, file = "PDZ protein gene set GSVA/gs.exp_5.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
group2.0 <- read.table("group2.0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group3.2 <- read.table("group3.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group3.0$type <- group2.0$type

conames <- interaction()
gs.exp0 <- cbind(gs.exp,gs.exp2)
group2.0 <- rbind(group2,group2.2)
group3.0 <- rbind(group3,group3.2)
write.table(group2.0, file = "group2.0.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(group3.0, file = "group3.0.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



##相关性分析
c_r=cor(as.numeric(exp2[i,]),as.numeric(group[r,]),method="pearson",use = 'pairwise.complete.obs')
p=cor.test(as.numeric(exp2[i,]),as.numeric(group[r,]),method ="pearson")[[3]]


gs.exp <- read.table("GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- read.table("group2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group3 <- read.table("group3.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group2 <- group2[order(group2$PDZK1),]
group3 <- group3[rownames(group2),]


group <- cbind(group2,group3)
group <- group[,-2]
group <- as.data.frame(t(group))
group.N <- group[group$group_list=="Normal",]
group.T <- group[group$group_list=="Tumor",]
group.N <- as.data.frame(t(group.N))
group.T <- as.data.frame(t(group.T))



gs.exp <- as.data.frame(t(gs.exp))
gs.exp.N <- gs.exp[,colnames(group.N)]
gs.exp.T <- gs.exp[,colnames(group.T)]



c_r=cor(as.numeric(exp2[i,]),as.numeric(group[r,]),method="pearson",use = 'pairwise.complete.obs')

LogFC <- gs.exp.T[,1,drop=F]
P <- gs.exp.T[,1,drop=F]

for (i in 1:234){
  for (r in 1:2){
    c_r=cor(as.numeric(gs.exp[i,]),as.numeric(group[r,]),method="pearson",use = 'pairwise.complete.obs')
    p=cor.test(as.numeric(gs.exp[i,]),as.numeric(group[r,]),method ="pearson")[[3]]
    ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
    LogFC[i,r+4] <- c_r
    P[i,r+4] <- p
    
    
  }
}
LogFC <- LogFC[,-1]
P <- P[,-1]
colnames(LogFC) <- c("PDZK1-T","MAGI3-T","PDZK1-N","MAGI3-N","PDZK1-A","MAGI3-A")
colnames(P) <- c("PDZK1-T","MAGI3-T","PDZK1-N","MAGI3-N","PDZK1-A","MAGI3-A")

write.table(LogFC, file = "cor_r-0.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(P, file = "P-0.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


LogFC3 <- read.table("LIHCcor_r-0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
P3 <- read.table("LIHCP-0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

LogFC2 <- LogFC2[rownames(LogFC),]
P2 <- P2[rownames(P),]

LogFC3 <- LogFC3[rownames(LogFC2),]
P3 <- P3[rownames(P2),]


LogFC <- cbind(LogFC2,LogFC3)

P <- cbind(P2,P3)

LogFC <- LogFC[,c(1,7,2,8,3,9,4,10,5,11,6,12)]

P <- P[,c(1,7,2,8,3,9,4,10,5,11,6,12)]

c_r=cor(as.numeric(gs.exp[,i]),as.numeric(group[,r]),method="pearson",use = 'pairwise.complete.obs')
p=cor.test(as.numeric(gs.exp[,i]),as.numeric(group[,r]),method ="pearson")[[3]]



pmt <- P
if (!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}

library(pheatmap)

pheatmap(LogFC,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=2.5,
         fontsize_col=6,
         display_numbers = pmt,
         main = "BRCA"
         )


LogFC1 <- read.table("PRADcor_r-0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
P1 <- read.table("PRADP-0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

LogFC2 <- read.table("BRCAcor_r-0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
P2 <- read.table("BRCAP-0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

LogFC3 <- read.table("KIRPcor_r-0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
P3 <- read.table("KIRPP-0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


geneset <- read.table("geneset分组.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

set <- geneset$geneset
LogFC0 <- cbind(LogFC1,LogFC2)

P0 <- cbind(P1,P2)

LogFC0 <- LogFC0[,c(1,7,2,8,3,9,4,10,5,11,6,12)]

P0 <- P0[,c(1,7,2,8,3,9,4,10,5,11,6,12)]

LogFC0 <- LogFC0[set,]
P0 <- P0[set,]
group <- geneset
rownames(group) <- group$geneset
group <- group[,-1,drop=F]


pmt <- P0
if (!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}


pheatmap(LogFC0,
         annotation_row = group,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = F,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6,
         display_numbers = pmt,
         border = F,
         main = "PRAD-BRCA"
)
####GSEA制作表达矩阵模板####


exp <- read.table("KIRC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)

exp2 <- exp[,1,drop=F]
colnames(exp2) <- "gene symbol"

exp[,1] <- exp2[,1]
exp2$descreption <- "na"
exp2 <- exp2[,2,drop=F]
exp <- cbind(exp2,exp)
write.table(exp, file = "KIRC.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####相关性点图绘制####
library(ggplot2)
library(ggpubr)
library(ggpmisc)
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))



gs.exp2 <- read.table("GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group_PDZK1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


geneset <- read.table("geneset分组.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


sample <- intersect(colnames(gs.exp2),rownames(group))
gs.exp2 <- gs.exp2[,sample]

group <- group[sample,]
set <- geneset$geneset

gs.exp2 <- gs.exp2[set,]
gs.exp2 <- as.data.frame(t(gs.exp2))
group <- group[rownames(gs.exp2),]

gs.exp2 <- cbind(gs.exp2,group)


data <- gs.exp2
data <- gs.exp

gs.exp <- as.data.frame(t(gs.exp))
sample <- intersect(rownames(gs.exp),rownames(PDZK1))
gs.exp <- gs.exp[sample,]
PDZK1 <- PDZK1[sample,,drop=F]

gs.exp <- cbind(gs.exp,PDZK1)

colnames(data[199])
#KEGG_FATTY_ACID_METABOLISM
#KEGG_CELL_CYCLE
#KEGG_WNT_SIGNALING_PATHWAY



b <- ggplot(data, aes(x = HALLMARK_MYC_TARGETS_V2, y = PDZK1))
b + geom_point()+
  geom_smooth( color="#00AFBB",method = "lm") +
  geom_rug(aes(color="#00AFBB")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  
  ggpubr::stat_cor(aes(color="#00AFBB"))




# Scatter plot with regression line
b + geom_point()+
  geom_smooth(method = "lm", color = "black", fill = "lightgray") 

# Add a loess smoothed fit curve
b + geom_point()+
  geom_smooth(method = "loess", color = "black", fill = "lightgray")


# Change color and shape by groups (cyl)
b + geom_point(aes(color = group_list, shape = group_list))+
  geom_smooth(aes(color = group_list, fill = group_list), method = "lm") +
  geom_rug(aes(color =group_list)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  
  ggpubr::stat_cor(aes(color = group_list))
# Remove confidence region (se = FALSE)
# Extend the regression lines: fullrange = TRUE
b + geom_point(aes(color = cyl, shape = cyl)) +
  geom_rug(aes(color =cyl)) +
  geom_smooth(aes(color = cyl), method = lm, 
              se = FALSE, fullrange = TRUE)+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  ggpubr::stat_cor(aes(color = cyl), label.x = 3)



p1 <- ggboxplot(data,  x = "group_list", y = "KEGG_ERBB_SIGNALING_PATHWAY",
                color = "group_list", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+   #删除x轴坐标名称
  #labs(title = "TGCA-HCC")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))+
  theme(text = element_text(size = 20))
p1


####WGCNA####
#合并数据类型
library(sva)
library(tidyverse)



gs.exp.LUAD <- gs.exp
gs.exp.LUAD <- as.data.frame(gs.exp.LUAD)
gs.exp.LUSC <- as.data.frame(gs.exp.LUSC)


gs.exp.LUSC <- gs.exp
gs.exp.PRAD <- read.table("GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


gs.exp1 <- read.table("GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gs.exp2 <- read.table("GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gs.exp1 <- cbind(gs.exp.PRAD,gs.exp.LUAD)


dataExpr <- read.table("dataExpr.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


data <- cbind(gs.exp1,gs.exp.LUSC)
data <- cbind(gs.exp1,gs.exp2)
boxplot(data)

boxplot(gs.exp1)


merge_eset=inner_join(GSE48452,GSE89632,by="symbol")
boxplot(data,outline=F)


library(WGCNA)

#

library(reshape2)
library(stringr)

# 
options(stringsAsFactors = FALSE)

exprMat <- data

type = "unsigned"

# 相关性计算
# 官方推荐 biweight mid-correlation & bicor
# corType: pearson or bicor
# 为与原文档一致，故未修改
corType = "pearson"

corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- as.data.frame(t(data))

dim(dataExpr)
head(dataExpr)[,1:8]


nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

#层次聚类树

sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

abline(h= 140, col = "red");#手动输入
#Determine cluster under the line
clust= cutreeStatic(sampleTree, cutHeight = 140, minSize = 10)
table(clust)
#clust 1 contains the samples we want to keep.
keepSamples= (clust==1)
dataExpr= dataExpr[keepSamples, ]
nGenes= ncol(dataExpr)
nSamples= nrow(dataExpr)

#计算软阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")


power = sft$powerEstimate
power



net <- blockwiseModules(
  dataExpr,
  power = 6,
  maxBlockSize = ncol(dataExpr),
  corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
  networkType = "unsigned",
  TOMType = "unsigned", 
  minModuleSize = 5,    ##越大模块越少
  mergeCutHeight = 0.2, ##越大模块越少
  numericLabels = TRUE, 
  saveTOMs= TRUE,
  saveTOMFileBase= "femaleMouseTOM",
  verbose = 3)

table(net$colors)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)


exp <- dataExpr
exp <- as.data.frame(t(exp))
exp$functions <- rownames(exp)
exp$group <- moduleColors

colors_group <-exp[,c(1645,1646)] 

write.table(colors_group, file = "colors_group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
MEs = net$MEs

#模块之间相关性图
### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)


#TOM图
TOM=TOMsimilarityFromExpr(dataExpr,power=7)
dissTOM=1-TOM
## draw all genes 

geneTree = net$dendrograms[[1]]
plotTOM = dissTOM^7
diag(plotTOM)=NA

TOMplot(plotTOM,geneTree,moduleColors,
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
        main="Network heapmap plot")



# TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)

## Loading objects:
##   TOM

TOM <- as.matrix(TOM)

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

# 这一部分特别耗时，行列同时做层级聚类
TOMplot(plotTOM, net$dendrograms, moduleColors,
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
        main = "Network heatmap plot")

#关联表型数据
trait <- read.table("phe2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
trait <- trait[rownames(MEs),]
trait <- traitData
# 读入表型数据，不是必须的
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}
traitData <- trait
### 模块与表型数据关联
if (corType=="pearsoon") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)


dev.off()
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 1, 
               ySymbols = colnames(MEs_col), colorLabels = T, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 1, zlim = c(-1,1),
               main = paste("Module-trait relationships"))



#导出用于cytoscape
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)

# Export the network into edge and node list files Cytoscape can read
# threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
# cytoscape中再调整


cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste(exprMat, ".edges.txt", sep=""),
                               nodeFile = paste(exprMat, ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0.5,
                               nodeNames = probes, nodeAttr = moduleColors)



#作图展示


cor <- read.table("PDZK1.KIRC.cor.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
P <- read.table("PDZK1.KIRC.P.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


modTraitCor <- as.matrix(cor)

modTraitP <- as.matrix(P)

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(cor), 
               yLabels = rownames(cor), 
               cex.lab = 1, 
               ySymbols = rownames(cor), colorLabels = T, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 1, zlim = c(-1,1),
               main = paste("Module-trait relationships"))




module = "turquoise";
module = "blue";
module = "brown";

probes = colnames(dataExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule];


modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste( module,collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste( module,collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
)



####整理UCSC.fpkm文件####
#与counts几乎相同，fpkm不需进行log转换
#读取tsv文件
library(tidyverse)

fpkm1 = read.table(file = 'TCGA-LUSC.htseq_fpkm.tsv', sep = '\t', header = TRUE) 
rownames(fpkm1) <- fpkm1[,1]  
fpkm1 = fpkm1[,-1]
table(substr(colnames(fpkm1),14,16))
fpkm1 <- fpkm1[,substr(colnames(fpkm1),14,16)%in% c("01A","11A")]
table(substr(colnames(fpkm1),14,16))
rownames(fpkm1) <- substr(rownames(fpkm1),1,15)
fpkm <- fpkm1

Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] #只要编码RNA
#取行名交集
comgene <- intersect(rownames(fpkm),rownames(Ginfo))
fpkm <- fpkm[comgene,]
Ginfo <- Ginfo[comgene,]
fpkm$Gene <- as.character(Ginfo$genename)   #新增Gene Symbol
fpkm <- fpkm[!duplicated(fpkm$Gene),]   #去重复
rownames(fpkm) <- fpkm$Gene   #将行名变为Gene Symbol
fpkm <- fpkm[,-ncol(fpkm)]   #去除最后一列
#保存所有患者的fpkm文件
write.table(fpkm, file = "LUSC_fpkm_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存癌症患者的fpkm文件
tumor <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "01A"]
fpkm_01A <- fpkm[,tumor]
write.table(fpkm_01A, file = "LUSC_fpkm_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存正常样本的fpkm文件
normal <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "11A"]
fpkm_11A <- fpkm[,normal]
write.table(fpkm_11A, file = "LUSC_fpkm_mRNA_11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#整理完毕#

####GTEx文档####
library(tidyverse)


fpkm <- read.table("gtex_RSEM_gene_fpkm.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

sample <- read.table("kidney and liver.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
sample1 <- rownames(sample)


sample2 <- read.table("prostate and lung.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
sample2 <- rownames(sample2)
exp2 <- fpkm1[sample2,]


fpkm1 <- as.data.frame(t(fpkm1))
exp <- fpkm1[sample1,]
exp <- as.data.frame(t(exp))


Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo0 <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] #只要编码RNA
#取行名交集
fpkm1 <- fpkm

rownames(fpkm1) <- substr(rownames(fpkm1),1,15)
comgene <- intersect(rownames(fpkm1),rownames(Ginfo0))
fpkm1 <- fpkm1[comgene,]
Ginfo0 <- Ginfo0[comgene,]
fpkm1$Gene <- as.character(Ginfo0$genename)   #新增Gene Symbol

fpkm1 <- fpkm1[!duplicated(fpkm1$Gene),] 

rownames(fpkm1) <- fpkm1$Gene   #将行名变为Gene Symbol
fpkm1 <- fpkm1[,-ncol(fpkm1)]   #去除最后一列
write.table(fpkm1, file = "GTEx.exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

fpkm1 <- as.data.frame(t(fpkm1))
exp1 <- fpkm1[sample1,]
exp2 <- fpkm1[sample2,]

exp1 <- as.data.frame(t(exp1))
exp2 <- as.data.frame(t(exp2))

exp1 <- na.omit(exp1)
exp2 <- na.omit(exp2)


write.table(exp1, file = "exp_1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp2, file = "exp_2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

phe <- read.table("GTEX_phenotype.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

phe2 <- phe[colnames(exp1),]
phe3 <- phe[colnames(exp2),]

write.table(phe2, file = "group_1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(phe3, file = "group_2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#GSVA分析
library(tidyverse)
library(GSEABase)
library(GSVA)
library(pheatmap)



exp1 <- as.matrix(exp1)
gs.exp1 <- gsva(exp1, l, kcdf = "Gaussian", min.sz = 10)


exp2 <- as.matrix(exp2)
gs.exp2 <- gsva(exp2, l, kcdf = "Gaussian", min.sz = 10)

write.table(gs.exp1, file = "gs.exp1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(gs.exp2, file = "gs.exp2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



pheatmap(gs.exp2,
         annotation_col=phe3,
         scale = "none",
         show_rownames =T ,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = T,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=4,
         main = "TCGA.LUSC"
)

#limma
library(limma) 

exp1 <- read.table("exp_1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group1 <- read.table("group_1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp2 <- read.table("exp_2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- read.table("group_2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group_list
design <- model.matrix(~0+factor(group2$`_primary_site`))
row.names(design) <- rownames(group2)

colnames(design)=c('Lung','Prostate')
contrast.matrix<-makeContrasts("Lung-Prostate",levels=design)

##step1
fit <- lmFit(gs.exp2,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "Lung-Prostate2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

set <- read.table("set.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
set <- set$set

gs.exp <- as.data.frame(gs.exp2)
gs.exp <- gs.exp[set,]


pheatmap(gs.exp,
         annotation_col=group2,
         scale = "none",
         show_rownames =T ,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = T,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=4,
         fontsize_col=4,
         main = ""
)


exp <- read.table("exp_1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group_1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gene <- read.table("DEG.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- gene$gene

exp2 <- exp[gene,]

exp2 <- na.omit(exp2)
group <- group[,-c(1,4)]
pheatmap(exp2,
         annotation_col=group,
         scale = "none",
         show_rownames =T ,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = T,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=4,
         fontsize_col=4,
         main = ""
)
#主成分降维分析
library (vegan) #加载vegan包
library (ggplot2)#加载ggplot包

color=c( '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000')

otu <- read.delim("normalize.txt", sep = '\t', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
otu <- exp
otu <- na.omit(otu)
otu<-t(otu)
##无效
#otu <- exp0[,(1:149)]

otu<-as.matrix(otu)
otu<- ExpressionSet(assayData = otu)
#处理缺失值和异常值
otu <- filter.NA(otu, thres = 0.25)#排除超过25%的测量缺失的基因
otu <- fill.NA(otu, mode = 'mean')
otu <- filter.std(otu, min.std = 0)

otu <- as.data.frame(otu)
#对比表达谱数据 Hellinger 转化前后的差异
otu_hel <- decostand(otu, method = 'hellinger')
#使用 Hellinger 转化后的数据
pca_sp2 <- rda(otu_hel, scale = FALSE)
pca_sp2 <- rda(otu, scale = FALSE)

#特征值提取




pca_exp2 <- pca_sp2$CA$eig / sum(pca_sp2$CA$eig)
pc1 <- paste('PC1:', round(pca_exp2[1]*100, 2), '%')
pc2 <- paste('PC2:',round(pca_exp2[2]*100, 2), '%')
pca2=pca_sp2[["CA"]][["u"]][,c(1,2)]
map2<-read.table(file="group.txt",header=T,sep="\t",row.names=1)
group <- map2$group1
group <- map2$group2
map2 <- group
map2$group1 <- map2$`_primary_site`

merged2<-merge(pca2,map2,by="row.names",all.x=TRUE)
merged2$`_primary_site` = factor(merged2$`_primary_site`,
                        levels = c("Liver","Kidney"))
merged2$group1 <- merged2$`_primary_site`
merged2$group2 = factor(merged2$group2,
                        levels = c("GSE48452","GSE164760",'GSE164760'))


p <- ggplot(data = merged2, aes(PC1, PC2)) +
  geom_point(aes(color = group1)) +  
  stat_ellipse(aes(fill =group1), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 
  scale_color_manual(values =color[1:length(unique(map2$group1))]) +
  scale_fill_manual(values = color[1:length(unique(map2$group1))]) +
  theme(panel.grid.major = element_line(color = 'gray', linewidth = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5), legend.title = element_blank()) +  
  labs(x = pc1, y = pc2) +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.5)
p

p <- ggplot(data = merged2, aes(PC1, PC2)) +
  geom_point(aes(color = group2)) +  
  stat_ellipse(aes(fill =group2), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 
  scale_color_manual(values =color[1:length(unique(map2$group2))]) +
  scale_fill_manual(values = color[1:length(unique(map2$group2))]) +
  theme(panel.grid.major = element_line(color = 'gray', linewidth = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5), legend.title = element_blank()) +  
  labs(x = pc1, y = pc2) +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.5)
p

PDZK1 <- exp["PDZK1",]
PDZK1 <- t(PDZK1)
PDZK1 <- PDZK1[rownames(group),,drop=F]
PDZK1 <- as.data.frame(PDZK1)
group$PDZK1 <- PDZK1$PDZK1

group$group <- group$`_primary_site`


liver_group <- group[group$group=="Liver",]
kidney_group <- group[group$group=="Kidney",]

liver_exp <- exp[,rownames(liver_group)]
kidney_exp <- exp[,rownames(kidney_group)]


liver_group$type <- factor(ifelse(liver_group$PDZK1>median(liver_group$PDZK1),"High","Low"))
kidney_group$type <- factor(ifelse(kidney_group$PDZK1>median(kidney_group$PDZK1),"High","Low"))

library(limma)
design <- model.matrix(~0+factor(liver_group$type))
row.names(design) <- rownames(liver_group)

colnames(design)=c('High','Low')
contrast.matrix<-makeContrasts("High-Low",levels=design)

##step1
fit <- lmFit(liver_exp,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "Liver_High-Low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#2
design <- model.matrix(~0+factor(kidney_group$type))
row.names(design) <- rownames(kidney_group)

colnames(design)=c('High','Low')
contrast.matrix<-makeContrasts("High-Low",levels=design)

##step1
fit <- lmFit(kidney_exp,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "kidney_High-Low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#PDZK1差异基因热图

exp <- read.table("exp_1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group_1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gene <- read.table("gene_PDZK1_deg.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

genelist <- gene$gene
exp2 <- exp[genelist,]
group <- group[,-c(1,4,5)]
group2 <- rbind(liver_group,kidney_group)

group <- group[rownames(group2),]
group$type <- group2$type

write.table(group, file = "group2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
group <- read.table("group2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group$



library(pheatmap)
exp2 <- exp2[,rownames(group)]
exp2 <- na.omit(exp2)
gene <- gene[rownames(exp2),,drop=F]
rownames(gene) <- gene$gene
gene <- gene[,-1,drop=F]

group <- group[,-4]
group <- group[rownames(group2),]
group$type <- group2$type



pheatmap(exp2,
         annotation_col= group,
         annotation_row = gene,
         scale = "row",
         show_rownames = T,
         show_colnames =T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks=bk,
         cluster_cols =F,
         cluster_rows = F,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6,
         gaps_col = 110,
         gaps_row = c(1075,1293)
         )

bk <- c(seq(-4,0,by=0.01),seq(0.001,4,by=0.01))
pheatmap(LogFC,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         cluster_cols =F,
         cluster_rows = T,
         border=FALSE,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=5,
         display_numbers = pmt,
         main = "TCGA.KIRC",
         
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks=bk
)


#GOBP、KEGG绘图

library(tidyverse)
library(AnnotationDbi)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)



hh <- read.table("common-KEGG.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
hh <- hh[c(1:30),]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Term))
ggplot(hh,aes(y=order,x=Count,fill=PValue))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "#FD8D62",high ="#66C3A5" )+#颜色自己可以换
  labs(
    x = "Gene numbers", 
    y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*PValue))+# 修改点的大小
  scale_color_gradient(low="#66C3A5",high = "#FD8D62")+
  labs(color=expression(PValue,size="Count"), 
       x="Gene Number",y="Pathways")+
  theme_bw()

hh <- read.table("GOBP.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Term))
ggplot(hh,aes(y=order,x=Count,fill=PValue))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "#FD8D62",high ="#66C3A5" )+#颜色自己可以换
  labs(
    x = "Gene numbers", 
    y = "GO term")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()


##肾、肝、乳腺三者的比较
exp <- read.table("GTEx.exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
sample1 <- read.table("kidney and liver.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
sample2 <- read.table("Breast_group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

sample <- read.table("B.K.L_group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

exp0 <- exp[,sample0]

sample <- sample[sample0,]

sample0 <- rownames(sample)

sample0 <- intersect(rownames(sample),colnames(exp))

sample0


write.table(sample, file = "B.K.L_group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp0, file = "B.K.L_exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#PCA
library (vegan) #加载vegan包
library (ggplot2)#加载ggplot包

color=c( '#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000')

otu <- exp0
otu <- na.omit(otu)
otu<-t(otu)
##无效
#otu <- exp0[,(1:149)]

otu<-as.matrix(otu)
otu<- ExpressionSet(assayData = otu)
#处理缺失值和异常值
otu <- filter.NA(otu, thres = 0.25)#排除超过25%的测量缺失的基因
otu <- fill.NA(otu, mode = 'mean')
otu <- filter.std(otu, min.std = 0)

otu <- as.data.frame(otu)
#对比表达谱数据 Hellinger 转化前后的差异
otu_hel <- decostand(otu, method = 'hellinger')
#使用 Hellinger 转化后的数据
pca_sp2 <- rda(otu_hel, scale = FALSE)
pca_sp2 <- rda(otu, scale = FALSE)

#特征值提取




pca_exp2 <- pca_sp2$CA$eig / sum(pca_sp2$CA$eig)
pc1 <- paste('PC1:', round(pca_exp2[1]*100, 2), '%')
pc2 <- paste('PC2:',round(pca_exp2[2]*100, 2), '%')
pca2=pca_sp2[["CA"]][["u"]][,c(1,2)]
map2<-read.table(file="group.txt",header=T,sep="\t",row.names=1)
group <- map2$group1
group <- map2$group2
map2 <- sample
map2$group1 <- map2$`_primary_site`

merged2<-merge(pca2,map2,by="row.names",all.x=TRUE)
merged2$`_primary_site` = factor(merged2$`_primary_site`,
                                 levels = c("Liver","Kidney","Breast"))
merged2$group1 <- merged2$`_primary_site`
merged2$group2 = factor(merged2$group2,
                        levels = c("GSE48452","GSE164760",'GSE164760'))


p <- ggplot(data = merged2, aes(PC1, PC2)) +
  geom_point(aes(color = group1)) +  
  stat_ellipse(aes(fill =group1), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 
  scale_color_manual(values =color[1:length(unique(map2$group1))]) +
  scale_fill_manual(values = color[1:length(unique(map2$group1))]) +
  theme(panel.grid.major = element_line(color = 'gray', linewidth = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5), legend.title = element_blank()) +  
  labs(x = pc1, y = pc2) +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.5)
p

p <- ggplot(data = merged2, aes(PC1, PC2)) +
  geom_point(aes(color = group2)) +  
  stat_ellipse(aes(fill =group2), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) + 
  scale_color_manual(values =color[1:length(unique(map2$group2))]) +
  scale_fill_manual(values = color[1:length(unique(map2$group2))]) +
  theme(panel.grid.major = element_line(color = 'gray', linewidth = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5), legend.title = element_blank()) +  
  labs(x = pc1, y = pc2) +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.5)
p

#breast——PDZK1差异分析
PDZK1 <- exp0["PDZK1",]
PDZK1 <- t(PDZK1)
PDZK1 <- PDZK1[rownames(sample),,drop=F]
PDZK1 <- as.data.frame(PDZK1)
sample$PDZK1 <- PDZK1$PDZK1

sample$group <- sample$`_primary_site`


liver_group <- group[group$group=="Liver",]
kidney_group <- group[group$group=="Kidney",]
breast_group <- sample[sample$group=="Breast",]


liver_exp <- exp[,rownames(liver_group)]
kidney_exp <- exp[,rownames(kidney_group)]
breast_exp <- exp0[,rownames(breast_group)]

liver_group$type <- factor(ifelse(liver_group$PDZK1>median(liver_group$PDZK1),"High","Low"))
kidney_group$type <- factor(ifelse(kidney_group$PDZK1>median(kidney_group$PDZK1),"High","Low"))

breast_group$type <- factor(ifelse(breast_group$PDZK1>median(breast_group$PDZK1),"High","Low"))


library(limma)
design <- model.matrix(~0+factor(breast_group$type))
row.names(design) <- rownames(breast_group)

colnames(design)=c('High','Low')
contrast.matrix<-makeContrasts("High-Low",levels=design)

##step1
fit <- lmFit(breast_exp,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "Breast_High-Low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



#热图
gene <- read.table("gene_PDZK1_DEG.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

genelist <- gene$gene
exp2 <- exp0[genelist,]
group <- group[,-c(1,4,5)]
group2 <- rbind(liver_group,kidney_group)

group <- group[rownames(group2),]
group$type <- group2$type

write.table(group, file = "group2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
group <- read.table("group2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group$
  
  
library(pheatmap)
exp2 <- exp2[,rownames(group)]
exp2 <- na.omit(exp2)
gene <- gene[rownames(exp2),,drop=F]
rownames(gene) <- gene$gene
gene <- gene[,-1,drop=F]

group <- group[,-4]
group <- group[rownames(group2),]
group$type <- group2$type



pheatmap(exp2,
         annotation_col= group,
         annotation_row = gene,
         scale = "row",
         show_rownames = T,
         show_colnames =T,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks=bk,
         cluster_cols =F,
         cluster_rows = F,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6,
         gaps_col = 110,
         gaps_row = c(1075,1293)
)

bk <- c(seq(-4,0,by=0.01),seq(0.001,4,by=0.01))
pheatmap(LogFC,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         cluster_cols =F,
         cluster_rows = T,
         border=FALSE,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=5,
         display_numbers = pmt,
         main = "TCGA.KIRC",
         
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks=bk
)

#GOBP,KEGG绘图
hh <- read.table("B,K,L_KEGG.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
hh <- hh[c(1:30),]
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Term))
ggplot(hh,aes(y=order,x=Count,fill=PValue))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "#FD8D62",high ="#66C3A5" )+#颜色自己可以换
  labs(
    x = "Gene numbers", 
    y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*PValue))+# 修改点的大小
  scale_color_gradient(low="#66C3A5",high = "#FD8D62")+
  labs(color=expression(PValue,size="Count"), 
       x="Gene Number",y="Pathways")+
  theme_bw()

hh <- read.table("GOBP.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Term))
ggplot(hh,aes(y=order,x=Count,fill=PValue))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "#FD8D62",high ="#66C3A5" )+#颜色自己可以换
  labs(
    x = "Gene numbers", 
    y = "GO term")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()




####基因集基因提取####

library(tidyr)
library(reshape2)
geneSet <- read.csv("genesets1.v2023.1.Hs.gmt",header = F,sep = "\t")
geneSet <- geneSet[,-2]
set_1 <- read.table("set_1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

set <- set_1$set
rownames(geneSet) <- geneSet$V2

geneSet <- geneSet[set,]
geneSet <- geneSet[,-c(1,2,3)]
geneSet <- t(geneSet)


geneSet2 <- geneSet %>%
  column_to_rownames("V1")%>%t()
a <- geneSet
a <- a[1:nrow(a),,drop=F]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}

save(l,file = "PDZ_protein.gene_set.Rdata")

geneSet <- as.data.frame(geneSet)
geneSet[[2]]
a <- list()
a <- geneSet[[1]]

a <- unite(geneSet, gene, c(1:87),sep = "\t")
b <- cbind(b[[1]],b[[2]])

geneSet$V1 <- rownames(geneSet)
b <- melt(geneSet,id.vars = "V1",na.rm=T)
write.table(b, file = "gene.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

KIRC.exp <- read.table("KIRC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
LIHC.exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
genelist <- read.table("genelist.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- genelist$gene

KIRC <- KIRC.exp[gene,]
LIHC <- LIHC.exp[gene,]
KIRC.group <- read.table("KIRC.group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
LIHC.group <- read.table("LIHC.group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group2 <- read.table("group2.0.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)



group <- rbind(KIRC.group,LIHC.group)

K_L.exp <- cbind(KIRC,LIHC)
K_L.exp <- na.omit(K_L.exp)
write.table(K_L.exp, file = "K_L.exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
K_L.group

write.table(group, file = "K_L.group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
exp <- read.table("K_L.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp2 <- exp[,rownames(group2)]
exp3 <- exp2[c(1:1500),]

library(pheatmap)

#删除标准差为0的行
exp2 = exp2[apply(exp2, 1, function(x) sd(x)!=0),] 

pheatmap(exp2,
         annotation_col=group2,
         scale = "row",
         show_rownames =T ,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=4,
         fontsize_col=4
)
bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))

pheatmap(exp2,
         annotation_col=group2,
         scale = "row",
         show_rownames = T,
         show_colnames =F,
         cluster_cols = F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=3,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2), colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         
         legend_breaks=seq(-5,5,2),
         
         breaks=bk)


#WGCNA

library(WGCNA)

#

library(reshape2)
library(stringr)

# 
options(stringsAsFactors = FALSE)

exprMat <- data

type = "unsigned"

# 相关性计算
# 官方推荐 biweight mid-correlation & bicor
# corType: pearson or bicor
# 为与原文档一致，故未修改
corType = "pearson"

corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- as.data.frame(t(exp))

dim(dataExpr)
head(dataExpr)[,1:8]


nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

#层次聚类树

sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

abline(h= 140, col = "red");#手动输入
#Determine cluster under the line
clust= cutreeStatic(sampleTree, cutHeight = 140, minSize = 10)
table(clust)
#clust 1 contains the samples we want to keep.
keepSamples= (clust==1)
dataExpr= dataExpr[keepSamples, ]
nGenes= ncol(dataExpr)
nSamples= nrow(dataExpr)

#计算软阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")


power = sft$powerEstimate
power



net <- blockwiseModules(
  dataExpr,
  power = 5,
  maxBlockSize = ncol(dataExpr),
  corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
  networkType = "unsigned",
  TOMType = "unsigned", 
  minModuleSize = 60,    ##越大模块越少
  mergeCutHeight = 0.2, ##越大模块越少
  numericLabels = TRUE, 
  saveTOMs= TRUE,
  saveTOMFileBase= "femaleMouseTOM",
  verbose = 3)

table(net$colors)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)


exp <- dataExpr
exp <- as.data.frame(t(exp))
exp$gene <- rownames(exp)
exp$group <- moduleColors

colors_group <-exp[,c(1018,1019)] 

write.table(colors_group, file = "colors_group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
MEs = net$MEs

#模块之间相关性图
### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)


#TOM图
TOM=TOMsimilarityFromExpr(dataExpr,power=7)
dissTOM=1-TOM
## draw all genes 

geneTree = net$dendrograms[[1]]
plotTOM = dissTOM^7
diag(plotTOM)=NA

TOMplot(plotTOM,geneTree,moduleColors,
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
        main="Network heapmap plot")



# TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)

## Loading objects:
##   TOM

TOM <- as.matrix(TOM)

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

# 这一部分特别耗时，行列同时做层级聚类
TOMplot(plotTOM, net$dendrograms, moduleColors,
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
        main = "Network heatmap plot")

#关联表型数据
trait <- read.table("phe2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
trait <- trait[rownames(MEs),]
trait <- traitData
# 读入表型数据，不是必须的
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}
traitData <- trait
### 模块与表型数据关联
if (corType=="pearsoon") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)


dev.off()
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 1, 
               ySymbols = colnames(MEs_col), colorLabels = T, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 1, zlim = c(-1,1),
               main = paste("Module-trait relationships"))



#导出用于cytoscape
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)

# Export the network into edge and node list files Cytoscape can read
# threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
# cytoscape中再调整


cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste(exprMat, ".edges.txt", sep=""),
                               nodeFile = paste(exprMat, ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0.5,
                               nodeNames = probes, nodeAttr = moduleColors)


#limma分析以及PDZK1相关性分析
library(limma)

KIRC.exp <- read.table("KIRC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
LIHC.exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

genelist <- read.table("genelist.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- genelist$gene

exp1 <- KIRC.exp[gene2,]
exp1 <- na.omit(exp1)
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

normal.group <- group[group$group_list=="Normal",]

tumor.group <- group[group$group_list=="Tumor",]


median(tumor.group$PDZK1)

normal.group$group <- ifelse(normal.group$PDZK1>=median(normal.group$PDZK1),"high","low")
tumor.group$group <- ifelse(tumor.group$PDZK1>=median(tumor.group$PDZK1),"high","low")




exp3 <- KIRC.exp[,rownames(tumor.group)] 
exp3 <- LIHC.exp[,rownames(tumor.group)] 

exp2 <- exp1[,rownames(normal.group)] 


exp <- exp3
group <- tumor.group

design <- model.matrix(~0+factor(group$group))
row.names(design) <- rownames(group)


colnames(design)=c('high','low')
contrast.matrix<-makeContrasts("high-low",levels=design)
exp <- exp[,rownames(group)]
##step1
fit <- lmFit(exp,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "tumor.high-low_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
nrDEG$gene <- rownames(nrDEG)

library(org.Hs.eg.db)

entriz <- mapIds(org.Hs.eg.db, keys = nrDEG$gene , keytype = "SYMBOL", column="ENTREZID")
nrDEG$gene_ID <- entriz
entriz
nrDEG <- na.omit(nrDEG)



gene2 <- read.table("PDZK1_normal.limma.gene.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene2 <- read.table("PDZK1_tumor.limma.gene,FC0.5.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gene2 <- gene2$gene

exp3 <- LIHC.exp[gene2,]
library(reshape2)
library(Hmisc)#加载包
exp2 <- exp3
exp2 <- t(exp2)
res2 <- rcorr(as.matrix(exp2))
res2
r <- res2[["r"]]
p <- res2[["P"]]
r <- melt(r)
p <- melt(p)
r$P <- p$value

write.table(r, file = "tumor.2R.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 可视化
library(PerformanceAnalytics)#加载包
chart.Correlation(exp3, histogram=T, pch=19)
chart
library(corrplot)
#绘制一个下三角的热图，这个包的使用在之前的博客写过，这里一笔带过
exp1 <- t(exp1)
cor_matr = cor(exp3)
cor_matr
corrplot(cor_matr, type="upper", order="hclust", tl.col="black", tl.srt=30,addCoef.col = "grey60")



rownames(cor_matr) = paste(1:14, colnames(cor_matr),sep = ' ')
colnames(cor_matr) = as.character(1:14) # 设置矩阵的列变量名为数字

col3 <- colorRampPalette(c('DodgerBlue3','white', "OrangeRed")) # 返回一个函数col3, 用来设置相关矩阵图的colorbar的分段
corrplot(cor_matr, method = 'circle', diag = F, type = 'full', outline = F,
         col = col3(20), cl.lim = c(-1,1),addgrid.col = NA,
         tl.pos = 'lb',tl.cex = 0.75, tl.col = 'black', tl.srt = 0, tl.offset = 0.5) # 绘制相关矩阵
axis(1,at = 1:14, labels = NA, pos = 14.5, tck = -0.01) # 在图片上边添加坐标轴，设置其刻度为位置，标签，画的位置，刻度的朝向
axis(4,at = 1:14, labels = NA, pos = 0.5, tck = -0.01) # 在图片左边添加坐标轴，设置参数同上

####相关性R分析####

LogFC <- exp[,1,drop=F]
P <- exp[,1,drop=F]
exp <- t(exp)
exp <- as.matrix(exp)
exp <- as.data.frame(exp)
exp <- exp[rownames(group),]
group <- as.matrix(group)


LogFC <- exp[,1,drop=F]
P <- exp[,1,drop=F]

for (i in 1:19620){
  
    c_r=cor(as.numeric(exp[,i]),as.numeric(group[,1]),method="pearson",use = 'pairwise.complete.obs')
    p=cor.test(as.numeric(exp[,i]),as.numeric(group[,1]),method ="pearson")[[3]]
    ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
    LogFC[i,2] <- c_r
    P[i,2] <- p
    
    
  }

LogFC <- LogFC[,-1,drop=F]
LogFC <- as.data.frame(LogFC)
rownames(LogFC) <- rownames(P)
LogFC$P <- P$V2
P <- P[,-1,drop=F]
colnames(LogFC) <- c("PDZK1-T","MAGI3-T","PDZK1-N","MAGI3-N","PDZK1-A","MAGI3-A")
colnames(P) <- c("PDZK1-T","MAGI3-T","PDZK1-N","MAGI3-N","PDZK1-A","MAGI3-A")
LogFC <- data2
write.table(LogFC, file = "LIHC_PDZK1_cor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(P, file = "P-0.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


LogFC$gene <- rownames(LogFC)

library(org.Hs.eg.db)

entriz <- mapIds(org.Hs.eg.db, keys = LogFC$gene , keytype = "SYMBOL", column="ENTREZID")
LogFC$gene_ID <- entriz
entriz
LogFC <- na.omit(LogFC)




####KEGGview绘图####

BiocManager::install("pathview")
library(pathview)
data(gse16873.d)
head(gse16873.d)




data <- read.table("KIRC_PDZK1_cor.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
data2 <- read.table("LIHC_PDZK1_cor.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
data3 <- read.table("KIRCtumor.high-low_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
data4 <- read.table("LIHCtumor.high-low_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


rownames(data) <- data$gene_ID
data2 <- data2[rownames(data),]

data <- cbind(data,data2)
rownames(data) <- data$gene_ID
data <- data[,c(1,5)]
gse16873.d <- read.table("KIRC2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gse16873.d <- as.matrix(gse16873.d)
#必须为矩阵形式
data <- as.matrix(data)
p <- pathview(gene.data = gse16873.d[,2], pathway.id = "04630", species = "hsa",
              out.suffix = "JAK_STAT_LIHC", kegg.native = T)
p


exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp <- read.table("KIRC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- read.table("group2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- group2[rownames(group),]
group <- cbind[group,group2]
group$grade <- group2$neoplasm_histologic_grade
group$stage <- group2$tumor_stage.diagnoses
group0 <- group

group <- group[group$group_list=="Tumor",]

group.G1 <- group[group$grade=="G1",]
group.G2 <- group[group$grade=="G2",]
group.G3 <- group[group$grade=="G3",]
group.G4 <- group[group$grade=="G4",]

group.S1 <- group[group$stage=="stage i",]
group.S2 <- group[group$stage=="stage ii",]
group.S3 <- group[group$stage=="stage iii",]
group.S4 <- group[group$stage=="stage iv",]






exp.G1 <- exp[,rownames(group.G1)]
exp.G2 <- exp[,rownames(group.G2)]
exp.G3 <- exp[,rownames(group.G3)]
exp.G4 <- exp[,rownames(group.G4)]

exp.S1 <- exp[,rownames(group.S1)]
exp.S2 <- exp[,rownames(group.S2)]
exp.S3 <- exp[,rownames(group.S3)]
exp.S4 <- exp[,rownames(group.S4)]




LogFC <- exp[,1,drop=F]
P <- exp[,1,drop=F]



for (i in 1:19620){
  
  c_r=cor(as.numeric(exp.S4[i,]),as.numeric(group.S4[,1]),method="pearson",use = 'pairwise.complete.obs')
  p=cor.test(as.numeric(exp.S4[i,]),as.numeric(group.S4[,1]),method ="pearson")[[3]]
  ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
  LogFC[i,5] <- c_r
  P[i,5] <- p
  
  
}



LogFC <- LogFC[,-1,drop=F]
LogFC <- as.data.frame(LogFC)
rownames(LogFC) <- rownames(P)
LogFC$P <- P$V2
P <- P[,-1,drop=F]
colnames(LogFC) <- c("G1","G2","G3","G4")
colnames(P) <- c("G1","G2","G3","G4")
colnames(LogFC) <- c("S1","S2","S3","S4")
colnames(P) <- c("S1","S2","S3","S4")

LogFC <- data2
write.table(LogFC, file = "KIRC.stage_PDZK1_cor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(P, file = "KIRC.stage_PDZK1_P.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

library(org.Hs.eg.db)

LogFC$gene <- rownames(LogFC)

entriz <- mapIds(org.Hs.eg.db, keys = LogFC$gene , keytype = "SYMBOL", column="ENTREZID")
LogFC$gene_ID <- entriz
entriz
LogFC <- na.omit(LogFC)
rownames(LogFC) <- LogFC$gene_ID

data <- LogFC
p <- pathview(gene.data = data[, 1:4], pathway.id = "04010", species = "hsa",
              out.suffix = "all_2", kegg.native = T)
p

#相关性点图绘制

library(ggplot2)
library(ggpubr)
library(ggpmisc)
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))

b <- ggplot(data, aes(x = KEGG_ERBB_SIGNALING_PATHWAY, y = MAGI3))
# Scatter plot with regression line
b + geom_point()+
  geom_smooth(method = "lm", color = "black", fill = "lightgray") 

# Add a loess smoothed fit curve
b + geom_point()+
  geom_smooth(method = "loess", color = "black", fill = "lightgray")


# Change color and shape by groups (cyl)
b + geom_point(aes(color = group_list, shape = group_list))+
  geom_smooth(aes(color = group_list, fill = group_list), method = "lm") +
  geom_rug(aes(color =group_list)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  
  ggpubr::stat_cor(aes(color = group_list))
# Remove confidence region (se = FALSE)
# Extend the regression lines: fullrange = TRUE
b + geom_point(aes(color = cyl, shape = cyl)) +
  geom_rug(aes(color =cyl)) +
  geom_smooth(aes(color = cyl), method = lm, 
              se = FALSE, fullrange = TRUE)+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  ggpubr::stat_cor(aes(color = cyl), label.x = 3)


b + geom_point(aes(color = cyl, shape = cyl))+
  geom_smooth(aes(color = cyl,fill = cyl),
              method = "lm", fullrange = TRUE)+
  facet_wrap(~cyl)+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  theme_bw()




data2 <- group
data2$grade <- factor(data2$grade,levels = c("G1","G2","G3","G4"))
data2$stage <- factor(data2$stage,levels = c("stage_i","stage_ii","stage_iii","stage_iv"))
data2 <- data2[,-4]
levels(data2$stage) <- c("stage i","stage ii","stage iii","stage iv")

p1 <- ggboxplot(data2,  x = "grade", y = "PDZK1",
                color = "grade", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+   #删除x轴坐标名称
  #labs(title = "TGCA-HCC")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))+
  theme(text = element_text(size = 20))
p1

p1 <- ggboxplot(data2,  x = "stage", y = "PDZK1",
                color = "stage", palette = "Dark2",
                add = "jitter")+             #alpha是调整透明度,也可以放到aes里，区别在于legend中是否显示
  theme_bw()+ theme(legend.position = "none")+#删除标签
  scale_fill_brewer(palette = "Set2")+ 
  stat_compare_means()+   #删除x轴坐标名称
  #labs(title = "TGCA-HCC")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(family = "serif"))+
  theme(text = element_text(size = 20))
p1





write.table(data2, file = "grade_PDZK1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
data2 <- read.table("grade_PDZK1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
data2 <- read.table("stage_PDZK1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)



data1 <- read.table("KIRC_PDZK1_cor.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

data2 <- read.table("LIHC_PDZK1_cor.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

rownames(data1) <- data1$gene_ID
rownames(data2) <- data2$gene_ID

data1$gene_ID <- as.character(data1$gene_ID)
library(pathview)
data1 <- as.matrix(data1)
data <- data1
p <- pathview(gene.data = data[,1], pathway.id = "04630", species = "hsa",
              out.suffix = "JAK_STAT_KIRC", kegg.native = T)




#正常和肿瘤里的差异对照
exp <- read.table("KIRC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
library(limma)

design <- model.matrix(~0+factor(group$group_list))
row.names(design) <- rownames(group)

exp <- exp[,rownames(group)]
colnames(design)=c('Normal','Tumor')
contrast.matrix<-makeContrasts("Tumor-Normal",levels=design)

##step1
fit <- lmFit(exp,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "Tumor-Normal.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

library(org.Hs.eg.db)

nrDEG$gene <- rownames(nrDEG)

entriz <- mapIds(org.Hs.eg.db, keys = nrDEG$gene , keytype = "SYMBOL", column="ENTREZID")
nrDEG$gene_ID <- entriz
entriz
nrDEG <- na.omit(nrDEG)
rownames(nrDEG) <- nrDEG$gene_ID



data1 <- read.table("Tumor-Normal.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

data2 <- data2[rownames(data1),]


data1$logFC2 <- data2$logFC

data <- data1[,c(1,9)]

library(pathview)
data2 <- nrDEG
data <- as.matrix(data)
gids= as.character(df[，1])

class(rownames(data))
data$gene_ID <- as.character(data$gene_ID)
data <- data[,c(8,1)]
rownames(data) <- data$gene_ID
data$logFC <- as.numeric(data$logFC)
data <- as.matrix(data)
p <- pathview(gene.data = data[,1:2], pathway.id = "04630", species = "hsa",
              out.suffix = "Tumor_normal", kegg.native = T,limit = list(gene=2),bins = list(gene=15))
data(gse16873.d)
data <- read.table("KIRC.grade_PDZK1_cor.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
data <- data[,1,drop=F]
data[,8] <- as.character(data[,8])
data <- as.matrix(data)
data[,1] <- scale(data[,1])
p <- pathview(gene.data = gse16873.d[,1], pathway.id = "04630", species = "hsa",
              out.suffix = "Tumor_normal", kegg.native = T,gene.idtype = "entrez")

data <- as.data.frame(data)
data <- data[rownames(gse16873.d),]
data <- as.matrix(data)
data <- na.omit(data)
rownames(gse16873.d)
rownames(data)

coname <- interaction(rownames(data),rownames(gse16873.d))

head(data)







####生存分析预后和GSVA评分的关系####
library(survival)
library(survminer)
surv <- read.table("TCGA-BRCA.survival.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv$OS.time <- surv$OS.time/30

for (i in vars){
  
  splots <- list()
  km_fit <- survfit(Surv(Survival_months,Vital_Status)~mydata[[i]], data=mydata)
  splots[[1]]<-ggsurvplot(km_fit,
                          xlab = "Time,mo",
                          ylab="Proportion Alive",
                          pval = T,
                          conf.int = F,##置信带
                          risk.table = T,
                          legend.title = i,
                          legend.labs = levels(mydata[[i]]),##
                          #surv.median.line = "hv",# 中位生存
                          palette="lancet")
  ## width=6.95,height=6.5
  res<-arrange_ggsurvplots(splots, print = F,
                           ncol = 1, nrow = 1, risk.table.height = 0.25)
  ggsave(paste(i,"All_surv.pdf",sep = "_"), res,width=7,height = 6)
  
  # Arrange multiple ggsurvplots and print the output
}

exp <- read.table("BRCAGSEA.exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 构建共同样本对象 

sample <- intersect(colnames(exp),rownames(surv))
exp <- exp[,sample]
surv <- surv[sample,]
exp <- as.data.frame(t(exp))
surv <- cbind(exp,surv)


# 循环代码

plist<-list()
data <- surv
pvalue0 <- exp[1,]

for(i in colnames(data)[1:234])
{
  group <- ifelse(data[[i]] > median(data[[i]]),"high","low")
  diff <- survdiff(Surv(OS.time,OS) ~ group, data = data)
  pValue = 1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(OS.time,OS) ~ group,data = data)
  p <- ggsurvplot(fit,
                  data = data,
                  conf.int = TRUE,
                  pval = pValue,
                  pval.size = 5,
                  legend.labs = c("High","Low"),
                  legend.title = i,
                  xlab = "Time(months)",
                  ylab = "Overall survival",
                  break.time.by = 20,
                  risk.table.title = "",
                  palette = c("#d7191c","#2b83ba"),
                  risk.table = F,
                  risk.table.height = .25
  )
  #pdf(file=paste0("surv/",i,".pdf"),onefile = FALSE,width = 6.5,height = 5.5)
  #print(p)
  #dev.off()
  plist[[i]]<-p$plot
  pvalue0[i] <- pValue
}


write.table(pvalue0, file = "KIRC_surv_Pvalue.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
save(plist,file = "BRCA_plot.Rdata")

write.table(pvalue0, file = "LIHC_surv_Pvalue.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
save(plist,file = "LIHC_plot.Rdata")

library(cowplot)

a <- plist[[1]]$plot
a
b <- plist[[2]]$plot
# 展示前6个gene的箱线图
plot<-plot_grid(a,b,
                ncol=1,align = 'v')
plot
dev.off()
#出不了图片专用
dev.new()
plist[[6]]

plot<-plot_grid(plist[[1]],plist[[2]],plist[[3]],
                plist[[4]],plist[[5]],plist[[6]],
                plist[[7]],plist[[8]],plist[[9]],
                plist[[10]],plist[[11]],plist[[12]],
                plist[[13]],plist[[14]],
                ncol=4,align = 'v')

plot<-plot_grid(plist[[1]],plist[[2]],plist[[3]],
                plist[[4]],plist[[5]],plist[[6]],
                plist[[7]],plist[[8]],plist[[9]],
                plist[[10]],plist[[11]],plist[[12]],
                plist[[13]],plist[[15]],
                plist[[16]],plist[[17]],plist[[21]],
                
                ncol=4,align = 'v')

plot

save(plist,file = "plist.Rdata")

plot<-plot_grid(plist[[1]],plist[[2]],plist[[3]],
                plist[[4]],plist[[5]],plist[[6]],
                plist[[7]],plist[[8]],plist[[9]],
                plist[[10]],plist[[11]],plist[[12]],
                plist[[13]],plist[[14]],plist[[15]],
                plist[[16]],plist[[17]],plist[[18]],
                plist[[19]],plist[[20]],plist[[21]],
                plist[[22]],
                ncol=4,align = 'v')
a <- plist[[1]]
b <- plist[[2]]

plot_grid(a,b,ncol=1,align = 'v')
library(ggpubr)
library(purrr)
ggarrange(plotlist = plist,nrow = 4,ncol = ceiling(length(plist)/4))

par(mfrow=c(2,7))

plist[[1]]


#循环测试

i <- colnames(data)[125]

group <- ifelse(data[[i]] > median(data[[i]]),"high","low")
diff <- survdiff(Surv(OS.time,OS) ~ group, data = data)
pValue = 1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(OS.time,OS) ~ group,data = data)
p <- ggsurvplot(fit,
                data = data,
                conf.int = TRUE,
                pval = pValue,
                pval.size = 5,
                legend.labs = c("High","Low"),
                legend.title = i,
                xlab = "Time (Months)",
                ylab = "Overall survival",
                break.time.by = 20,
                risk.table.title = "",
                palette = c("#d7191c","#2b83ba"),
                risk.table = F,
                risk.table.height = .25
)

p#res<-arrange_ggsurvplots(p, print = F,
#ncol = 1, nrow = 1, risk.table.height = 0.25)
#ggsave(paste(i,"Overall_survival.pdf",sep = "_"), res,width=7,height = 6)
plist[[i]]<-p
pvalue0[i] <- pValue

plot(plist[[88]])

####GEO数据库数据验证####

#GSVA评分
#GSE73731
exp <- read.table("exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


exp <- log2(exp)
PDZK1 <- exp["PDZK1",]

PDZK1 <- as.data.frame(t(PDZK1))

PDZK1 <- PDZK1[rownames(group),,drop=F]

group$PDZK1 <- PDZK1$PDZK1

exp <- as.data.frame(exp)
gs.exp2 <- gs.exp
exp <- as.matrix(exp)

gs.exp <- gsva(exp, l, kcdf = "Gaussian", min.sz = 10)

write.table(gs.exp, file = "GSEA.exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

BiocManager::install(version = "3.18")


####乳腺癌分型分析####
#PDZK1相关性计算

exp <- read.table("exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("group2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


surv <- read.table("TCGA-BRCA.survival.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


A_group <- group[group$subtype=="LumA",]

B_group <- group[group$subtype=="LumB",]

H_group <- group[group$subtype=="Her2",]

Ba_group <- group[group$subtype=="Basal",]

No_group <- group[group$subtype=="Normal",]
write.table(A_group, file = "A_group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(B_group, file = "B_group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(H_group, file = "H_group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Ba_group, file = "Ba_group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(No_group, file = "No_group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

gs.exp <- exp

gs.exp <- gs.exp[,rownames(No_group)]
gs.exp <- as.data.frame(t(gs.exp))
group <- No_group


LogFC <- gs.exp[,c(1,2,3,4,5),drop=F]
P <- gs.exp[,c(1,2,3,4,5),drop=F]

for (i in 1:234){
  for (r in 1:2){
    c_r=cor(as.numeric(gs.exp[,i]),as.numeric(group[,1]),method="pearson",use = 'pairwise.complete.obs')
    p=cor.test(as.numeric(gs.exp[,i]),as.numeric(group[,1]),method ="pearson")[[3]]
    ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
    LogFC[i,5] <- c_r
    P[i,5] <- p
  }
}


LogFC <- LogFC[,-1]
P <- P[,-1]
colnames(LogFC) <- c("LumA","LumB","Her2","Basal","Normal")
colnames(P) <- c("LumA","LumB","Her2","Basal","Normal")
rownames(LogFC) <- rownames(exp)
rownames(P) <- rownames(exp)

write.table(LogFC, file = "cor_r-All.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(P, file = "P-All.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

logFC <- read.table("cor_r-All.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

P <- read.table("P-All.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


library(pheatmap)
LogFC1 <- LogFC[,-1]
LogFC <- LogFC[,-23]
P <- P[,-23]
P1 <- P[,-1]
pmt <- P
if (!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}

logFC1 <- logFC[c(1:117),]

logFC2 <- logFC[c(118:234),]

pmt1 <- pmt[c(1:117),]
pmt2 <- pmt[c(118:234),]
pmt0 <- pmt

group2 <- read.table("group4.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

p1 <- pheatmap(logFC1,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = F,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6,
         cellwidth = 5, cellheight = 5,
         display_numbers = pmt1)

p2 <- pheatmap(logFC2,
               annotation_col=NA,
               scale = "none",
               show_rownames = T,
               show_colnames =T,
               color = colorRampPalette(c("navy", "white", "red"))(50),
               cluster_cols =F,
               cluster_rows = F,
               fontsize = 10,
               fontsize_row=5,
               fontsize_col=6,
               cellwidth = 5, cellheight = 5,
               display_numbers = pmt2)


cowplot::plot_grid(p1$gtable,p2$gtable,ncol= 2, labels=LETTERS[1:2])





pheatmap(LogFC1,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = F,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6
)


logFC_1 <- read.table("cor_r-1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

P_1 <- read.table("P-1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

pmt <- P_1
if (!is.null(pmt)){
  ssmt <- pmt< 0.01
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}



pheatmap(logFC_1,
         annotation_col=NA,
         scale = "none",
         show_rownames = T,
         show_colnames =T,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = F,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6,
         cellwidth = 15, cellheight = 15,
         display_numbers = pmt)


gs.exp <- read.table("BRCAGSEA.exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

group <- read.table("group3.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

gs.exp <- gs.exp[,rownames(group)]

gs.exp <- as.data.frame(t(gs.exp))

data <- cbind(gs.exp,group)

library(ggplot2)

b <- ggplot(data, aes(x = KEGG_FATTY_ACID_METABOLISM, y = PDZK1))

b + geom_point()+
  geom_smooth(method = "lm", color = "black", fill = "lightgray") 

# Add a loess smoothed fit curve
b + geom_point()+
  geom_smooth(method = "loess", color = "black", fill = "lightgray")


# Change color and shape by groups (cyl)
b <- b + geom_point(aes(color = subtype, shape = subtype))+
  geom_smooth(aes(color = subtype, fill = subtype), method = "lm") +
  geom_rug(aes(color =subtype)) +
  scale_color_manual(values = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231'))+
  scale_fill_manual(values = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231'))+ 
  
  ggpubr::stat_cor(aes(color = subtype))

ggMarginal(
  b,
  type = 'boxplot',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)



# Remove confidence region (se = FALSE)
# Extend the regression lines: fullrange = TRUE
b + geom_point(aes(color = cyl, shape = cyl)) +
  geom_rug(aes(color =cyl)) +
  geom_smooth(aes(color = cyl), method = lm, 
              se = FALSE, fullrange = TRUE)+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  ggpubr::stat_cor(aes(color = cyl), label.x = 3)



library(survival)
library(survminer)
surv <- read.table("TCGA-BRCA.survival.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv$OS.time <- surv$OS.time/30

surv <- surv[rownames(group),]

data <- cbind(gs.exp,surv)

data0 <- data
i <- "HALLMARK_KRAS_SIGNALING_UP"

A_data <- data[data$Var.241=="LumA",]
B_data <- data[data$Var.241=="LumB",]
H_data <- data[data$Var.241=="Her2",]
Ba_data <- data[data$Var.241=="Basal",]
No_data <- data[data$Var.241=="Normal",]

data <- A_data
data <- B_data
data <- H_data
data <- Ba_data
data <- No_data
data <- data0




group <- ifelse(data[[i]] > median(data[[i]]),"high","low")
diff <- survdiff(Surv(OS.time,OS) ~ group, data = data)
pValue = 1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(OS.time,OS) ~ group,data = data)
p6 <- ggsurvplot(fit,
                data = data,
                conf.int = TRUE,
                pval = pValue,
                pval.size = 5,
                legend.labs = c("High","Low"),
                legend.title = i,
                xlab = "Time (Months)",
                ylab = "Overall survival",
                break.time.by = 20,
                risk.table.title = "",
                palette = c("#d7191c","#2b83ba"),
                risk.table = F,
                risk.table.height = .25
)
p1

cowplot::plot_grid(p1$plot,p2$plot,p3$plot,p4$plot,p5$plot,p6$plot,ncol = 3,
                   labels = c('LumA', 'LumB',"Her2","Basal","Normal","All"))


####相关性点图2####

gs.exp <- read.table("KIRCGSEA.exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gs.exp <- read.table("LIHCGSEA.exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gs.exp <- read.table("gs.exp1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


group <- read.table("group2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
PDZK1 <- read.table("PDZK1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


gs.exp <- gs.exp[,rownames(group)]

gs.exp <- gs.exp[rownames(PDZK1),]


gs.exp <- as.data.frame(t(gs.exp))


data <- cbind(gs.exp,PDZK1)

data <- na.omit(data)

library(ggstatsplot)
library(ggExtra)
library(ggsci)

p <- ggscatterstats(data, 
               y =PDZK1, 
               x =KEGG_FATTY_ACID_METABOLISM,
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#009E73",
               yfill = "#D55E00", 
               marginal.type = "boxplot",    #类型可以换成density,boxplot,violin,densigram
               xparams = list(fill = group_list),
               yparams = list(fill = group_list),
               title = ""
               )

p


p <- p+ geom_point(aes(color = group_list, shape = group_list))+
  geom_smooth(aes(color = group_list, fill = group_list), method = "lm") +
  geom_rug(aes(color =group_list)) +
  scale_color_manual(values = c("#009E73", "#D55E00"))+
  scale_fill_manual(values = c("#009E73", "#D55E00"))+ 
  
  ggpubr::stat_cor(aes(color = group_list))

p


b <- ggplot(data, aes(x = KEGG_PPAR_SIGNALING_PATHWAY, y = PDZK1))
b + geom_point()+
  geom_smooth(method = "lm", color = "black", fill = "lightgray") 

# Add a loess smoothed fit curve
b + geom_point()+
  geom_smooth(method = "loess", color = "black", fill = "lightgray")


# Change color and shape by groups (cyl)


b <- b + geom_point(aes(color = VHL, shape = VHL))+
  geom_smooth(aes(color = VHL, fill = VHL), method = "lm") +
  geom_rug(aes(color =VHL)) +
  scale_color_manual(values = c("#2b83ba", "#FC4E07"))+
  scale_fill_manual(values = c("#2b83ba", "#FC4E07"))+ 
  
  ggpubr::stat_cor(aes(color = VHL))

ggMarginal(
  b,
  type = 'boxplot',
  margins = 'both',
  size = 10,
  groupColour = TRUE,
  groupFill = TRUE
)



b <- b + geom_point(aes(color = group_list, shape = group_list))+
  geom_smooth(aes(color = group_list, fill = group_list), method = "lm") +
  geom_rug(aes(color =group_list)) +
  scale_color_manual(values = c("#2b83ba", "#FC4E07"))+
  scale_fill_manual(values = c("#2b83ba", "#FC4E07"))+ 
  
  ggpubr::stat_cor(aes(color = group_list))

b


ggMarginal(
  b,
  type = 'boxplot',
  margins = 'both',
  size = 10,
  groupColour = TRUE,
  groupFill = TRUE
)


####GSEA富集分析####

library(limma) 

exp <- read.table("exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


group1 <- read.table("N_group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

group2 <- read.table("T_group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

group2 <- group[group$group=="tumor",]
data0 <- exp
group <- exp["PDZK1",]
group <- as.data.frame(t(group))

group <- group[rownames(group2),,drop = F]

group2$PDZK1 <- group$PDZK1

write.table(group2, file = "GSE76427group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


data0 <- data0[,rownames(group)]
N_data <- data0[,rownames(group1)]
T_data <- data0[,rownames(group2)]

group <- group2
group$PDZK1 <- t(data["PDZK1",])
colnames(group$PDZK1) <- "PDZK1"
group$group1 <- ifelse(group$PDZK1 > median(group$PDZK1),"high","low")

data<-T_data

design <- model.matrix(~0+factor(group$group1))


rownames(design) <- rownames(group)

colnames(design)=c("high","low")


#design数据表要符合都是数字形式0、1，data看情况要转置
contrast.matrix<-makeContrasts("high-low",levels=design)

##step1
fit <- lmFit(data,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "high-low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)






library(fgsea)
library(biomaRt)
library(enrichplot)
library(tidyverse)
#geneSet制作
geneSet <- read.csv("NASHed2.gmt",header = F,sep = "\t") # 用EXCEL打开删除NA列第二列
geneSet <- geneSet[,-c(1,3)]
geneSet <- geneSet %>%
  column_to_rownames("V2")%>%t()
a <- geneSet
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}
save(l,file = "./NASHed.Rdata")





c <- read.table("顺序指示文件.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
c[1,1]
i=1
for (i in c(1:10)){
  mydata <- read.table(c[2,i],sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  mydata <- mydata[order(mydata$logFC,decreasing=TRUE),]
  FCgenelist <- mydata$logFC
  names(FCgenelist) <- rownames(mydata)
  fgseaRes <- fgsea(pathways = l, 
                    stats = FCgenelist,
                    minSize=5,
                    maxSize=500,
                    nperm=1000)
  result <- fgseaRes[padj < 0.05,]
  result <- result[,-8]
  write.table(result,file=paste0(c[1,i],"_GSEA.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
  
}

#测试
mydata <- read.table("high-low.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mydata <- nrDEG

mydata <- mydata[order(mydata$logFC,decreasing=F),]
FCgenelist <- mydata$logFC
names(FCgenelist) <- rownames(mydata)
fgseaRes <- fgsea(pathways = l, 
                  stats = FCgenelist,
                  minSize=5,
                  maxSize=500,
                  nperm=1000)
result <- fgseaRes
result <- result[,-8]

write.table(result, file = "normal_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


data(examplePathways)
data(exampleRanks)


plotEnrichment(l[["KEGG_P53_SIGNALING_PATHWAY"]],
               ranks) + labs(title="KEGG_P53_SIGNALING_PATHWAY")


library(GseaVis)

data <- fgseaRes
data$pvalue <- data$pval
data$p.adjust <- data$padj

data <- data[,-8]
data$geneList <- data$nMoreExtreme

rownames(data) <- data$pathway
dotplotGsea(data = fgseaRes,topn = 10)

gseaplot2(data, title = "", geneSetID = 1)





gseaNb(object = data,
       geneSetID = 'KEGG_P53_SIGNALING_PATHWAY')




library(org.Hs.eg.db) # human的OrgDB

library(clusterProfiler)
# 可视化
library(enrichplot)
library(ggplot2)


remotes::install_git('https://gitee.com/swcyo/myenrichplot/')
remotes::install_git('https://gitee.com/swcyo/myenrichplot/')

mydata <- read.table("high-low_tumor.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mydata <- read.table("high-low_normal.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mydata <- read.table("high-low.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mydata <- nrDEG

data_sort <- mydata %>%
  arrange(desc(logFC))

gene_list <- data_sort$logFC
names(gene_list) <- rownames(data_sort)
head(gene_list)


gmt <- read.gmt("genesets1.v2023.1.Hs.gmt")

?GSEA
res <- GSEA(
  gene_list,
  pvalueCutoff=1,
  TERM2GENE = gmt
)

gseaplot2(res, title = res$Description[2], geneSetID = 2,pvalue_table = T)


gseaNb(object = res,
       geneSetID = 'KEGG_P53_SIGNALING_PATHWAY')


gseaNb(res,
       geneSetID = 'KEGG_P53_SIGNALING_PATHWAY',
       addPval = T,
       pCol = "steelblue",
       pvalX = 0.65, # 位置
       pvalY = 0.7,
       pHjust = 0, # 对齐方式
       nesDigit = 4, # 小数点位数
       pDigit = 4
)




save(res,file = "res_tumor.Rdata")

write.table(res, file = "result_tumor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



#多条通路一起展示
terms <- c("KEGG_FATTY_ACID_METABOLISM","KEGG_PPAR_SIGNALING_PATHWAY","HALLMARK_FATTY_ACID_METABOLISM")

terms <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING","KEGG_WNT_SIGNALING_PATHWAY","KEGG_MTOR_SIGNALING_PATHWAY",
           "KEGG_P53_SIGNALING_PATHWAY")

terms <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING","KEGG_WNT_SIGNALING_PATHWAY")
           
terms <- c("KEGG_MTOR_SIGNALING_PATHWAY","KEGG_P53_SIGNALING_PATHWAY")

terms <- c("KEGG_JAK_STAT_SIGNALING_PATHWAY","HALLMARK_KRAS_SIGNALING_DN")

terms <- c("HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2")

gseaNb(object = res,
       geneSetID = terms,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addPval = T,
       pvalX = 0.05,pvalY = 0.05,
       rmPrefix = F)

?gseaNb


gseaplot_list <- lapply(terms, function(x){
  gseaNb(object = res,
         geneSetID = x,
         termWidth = 30,
         addPval = T,
         pvalX = 0.75,
         pvalY = 0.6
  )
})
cowplot::plot_grid(plotlist=gseaplot_list, ncol = 2)

####PDZK1泛癌展示####

library(ggplot2)
library(ggpubr)#基于ggplot2的可视化包，主要用于绘制符合出版要求的图形
library(ggsignif)#用于P值计算和显著性标记
library(tidyverse)#数据预处理


exp <- read.table("PANCANCER_fpkm_mRNA_all_PDZ.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("group2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


exp <- exp[,rownames(group)]
  
PDZK1 <- exp["PDZK1",]
PDZK1 <- as.data.frame(t(PDZK1))
group <- cbind(group,PDZK1)

data <- group

p1 <- ggplot(data,aes(x=Progect,y=PDZK1,fill=group))+
  geom_boxplot(width=0.6,alpha=0.8)+
  scale_fill_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8))+
  geom_signif(comparisons = list(c("Tumor","Normal")),#设置需要比较的组
              test = t.test, ##计算方法
              #横线下方的竖线设置
              size=0.8,color="black")

p1


plot_order


#设置x轴分组顺序
plot_order = data[data$group=="Tumor",] %>% 
  group_by(Progect) %>% 
  summarise(m = median(PDZK1)) %>% 
  arrange(desc(m)) %>% 
  pull(Progect)

data$Progect = factor(data$Progect,levels = plot_order)





library(gapminder)
data %>%
  ggplot(aes(x=group,y=PDZK1,fill=group)) +
  geom_boxplot() + geom_jitter(width = 0.1,alpha = 0.2) +
  facet_wrap(~Progect,ncol = 8,nrow = 4) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+ 
  stat_compare_means(label.y = 0.4,label = "p.signif")+
  scale_fill_brewer(palette = "Set2")+ theme(text = element_text(size = 15))

####GSE126694####

library(tidyverse)
library(GEOquery)
gset = getGEO('GSE215286', destdir=".", AnnotGPL = F, getGPL = F)
gset[[1]]
#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
write.table(pdata, file = "phe.txt",sep = "\t",row.names = T,col.names = NA,quote = F)




exp <- read.table("GSE215286_norm_counts_FPKM.tsv",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

an <- read.csv("Human.GRCh38.p13.annot.tsv",sep = "\t",check.names = F,stringsAsFactors = F,header = T)


gene <- an[an$GeneType=="protein-coding",]
rownames(gene) <- gene$GeneID

exp1 <- exp[rownames(gene),]

exp1$gene <- gene$Symbol
rownames(exp1) <- exp1$gene
exp1 <- exp1[,-56]

boxplot(exp1,outline=FALSE)

write.table(exp1, file = "exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
exp1 <- as.data.frame(exp1)
PDZK1 <- exp1["PDZK1",]
PDZK1 <- as.data.frame(t(PDZK1))

library(tidyverse)
library(GSEABase)
library(GSVA)
library(pheatmap)



group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- pdata[,40,drop=F]
group <- group[rownames(PDZK1),,drop=F]
colnames(group) <- c("group")

group$PDZK1 <- PDZK1$PDZK1
exp1 <- as.matrix(exp1)
gs.exp <- gsva(exp1, l, kcdf = "Gaussian", min.sz = 10)
data <- t(gs.exp)
data <- as.data.frame(data)
data <- cbind(data,PDZK1)
data <- cbind(data,group)

write.table(gs.exp, file = "gsva.exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



library(ggstatsplot)
library(ggExtra)
library(ggsci)

b <- ggplot(data, aes(x = KEGG_FATTY_ACID_METABOLISM, y = PDZK1))

#FC4E07
b <- b + geom_point(aes(color = group, shape = group))+
  geom_smooth(aes(color = group, fill = group), method = "lm") +
  geom_rug(aes(color =group)) +
  scale_color_manual(values = c("#FC4E07", "#2b83ba"))+
  scale_fill_manual(values = c("#FC4E07", "#2b83ba"))+ 
  
  ggpubr::stat_cor(aes(color = group))

ggMarginal(
  b,
  type = 'boxplot',
  margins = 'both',
  size = 10,
  groupColour = TRUE,
  groupFill = TRUE
)


#相关性分析

group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
trems <- read.table("trems.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)



exp0 <- gs.exp
group <- group[group$group=="tumor",,drop=F]
exp0 <- as.data.frame(exp0)
exp0 <- exp0[,rownames(group)]
exp0 <- exp0[trems$trems,]
PDZK1 <- PDZK1[rownames(group),,drop=F]


group$PDZK1 <- PDZK1$PDZK1
Rfile <- exp0[,1,drop=F]
Pfile <- exp0[,1,drop=F]

i=1
c_r=cor(as.numeric(exp0[i,]),as.numeric(group[,"PDZK1"]),method="pearson",use = 'pairwise.complete.obs')
p=cor.test(as.numeric(exp0[i,]),as.numeric(group[,"PDZK1"]),method ="pearson")[[3]]
Rfile[i,1] <- c_r

for (i in 1:85){
  
    c_r=cor(as.numeric(exp0[i,]),as.numeric(group[,"PDZK1"]),method="pearson",use = 'pairwise.complete.obs')
    p=cor.test(as.numeric(exp0[i,]),as.numeric(group[,"PDZK1"]),method ="pearson")[[3]]
    ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
    Rfile[i,1] <- c_r
    Pfile[i,1] <- p
    
  
}

colnames(Rfile) <- c("TCGA.KIRC","GSE53757","GSE73731")
colnames(Pfile) <- c("TCGA.KIRC","GSE53757","GSE73731")

write.table(Pfile, file = "Pfile.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Rfile, file = "Rfile.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



#GSEA
library(limma) 

group <- group2
group$PDZK1 <- t(data["PDZK1",])
colnames(group$PDZK1) <- "PDZK1"
group$group1 <- ifelse(group$PDZK1 > median(group$PDZK1),"high","low")





exp1 <- exp1[,rownames(group)]
data<-exp1

design <- model.matrix(~0+factor(group$group1))


rownames(design) <- rownames(group)

colnames(design)=c("high","low")


#design数据表要符合都是数字形式0、1，data看情况要转置
contrast.matrix<-makeContrasts("high-low",levels=design)

##step1
fit <- lmFit(data,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "high-low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)






library(org.Hs.eg.db) # human的OrgDB

library(clusterProfiler)
# 可视化
library(enrichplot)
library(ggplot2)

mydata <- nrDEG

data_sort <- mydata %>%
  arrange(desc(logFC))

gene_list <- data_sort$logFC
names(gene_list) <- rownames(data_sort)
head(gene_list)


gmt <- read.gmt("genesets1.v2023.1.Hs.gmt")

?GSEA
res <- GSEA(
  gene_list,
  pvalueCutoff=1,
  TERM2GENE = gmt
)
save(res,file = "res_tumor.Rdata")

write.table(res, file = "result_tumor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#stage分级分析
exp <- read.table("gsva.exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp0 <- read.table("exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

group <- read.table("group2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


exp0 <- exp0[,rownames(group)]

PDZK1 <- exp0[,"PDZK1",drop=F]
exp0 <- as.data.frame(t(exp0))
data3 <- cbind(PDZK1,group)
data <- cbind(exp,group)




data$Stage <- factor(data$Stage,levels = c("T1a","T1b","T2a","T3a","T3b"))
data$grade <- factor(data$grade,levels = c('Grade1','Grade2','Grade3','Grade4'))

data3$grade <- factor(data3$grade,levels = c('Grade1','Grade2','Grade3','Grade4'))


data2 <- na.omit(data)
data3 <- na.omit(data3)
library(ggpubr) 

p <- ggboxplot(data3, x = "grade", y = "PDZK1",
               color = "grade", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0","#2b83ba","navy", '#eb4b3a'))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
#更换颜色

p + stat_compare_means()


p <- ggboxplot(data, x = "Stage", y = "KEGG_FATTY_ACID_METABOLISM",
               color = "Stage", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0","#2b83ba","navy", '#eb4b3a'))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
#更换颜色

p + stat_compare_means()

p <- ggboxplot(data2, x = "grade", y = "KEGG_FATTY_ACID_METABOLISM",
               color = "grade", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0","#2b83ba","navy", '#eb4b3a'))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
#更换颜色

p + stat_compare_means()







####TCGA stage 分组功能差异####

gs.exp <- read.table("GSEA.exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
trems <- read.table("trems.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("group_stage2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

exp0 <- read.table("GSE73731exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp0 <- read.table("GSE215286exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)



PDZK1 <- exp[c("PDZK1","ULK1"),]
PDZK1 <- as.data.frame(t(PDZK1))



PDZK1 <- PDZK1[rownames(group),,drop=F]
group$PDZK1 <- PDZK1$PDZK1
exp <- gs.exp
exp <- exp[trems$trems,rownames(group)]

exp <- exp[rownames(group),]
exp <- as.data.frame(t(exp))

data <- cbind(exp,group)
data <- cbind(data,PDZK1)
data$stage <- factor(data$stage,levels=c("Stage: 1","Stage: 2","Stage: 3","Stage: 4"))

data$stage <- factor(data$stage,levels=c("stage i/ii","stage iii/iv"))

data$stage <- factor(data$stage,levels=c("Stage: 1/2","Stage: 3/4"))

#KEGG_FATTY_ACID_METABOLISM
#KEGG_PPAR_SIGNALING_PATHWAY
#HALLMARK_FATTY_ACID_METABOLISM
#KEGG_WNT_SIGNALING_PATHWAY
#KEGG_MTOR_SIGNALING_PATHWAY
#KEGG_P53_SIGNALING_PATHWAY
#HALLMARK_IL6_JAK_STAT3_SIGNALING

data$cancer_status <- factor(data$cancer_statu, levels = c("TUMOR FREE","WITH TUMOR"))
data2 <- na.omit(data)

n_fun <- function(x){
  return(data.frame(y = median(x), label = paste0("n = ",length(x))))
}

library(ggpubr) 
p <- ggboxplot(data2, x = "cancer_status", y = "PDZK1",
               color = "cancer_status", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0","#2b83ba", '#eb4b3a',"red"))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
#更换颜色

p + stat_compare_means()+
  stat_summary(fun.data = n_fun, geom = "text")



p2 <- ggboxplot(data, x = "group1", y = "TNFRSF12A",
                color = "group1", palette = "jco",
                add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0", '#eb4b3a'))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
#更换颜色

p2 + stat_compare_means()


plot<-plot_grid(p,p2,
                ncol=,align = 'v')

plot


#mfuzz尝试聚类
sample1<-aggregate(data[,1:85],by=list(data$stage),mean,na.rm= TRUE)#排除最后的一列group

#更改顺序
sample1 <- sample1[c(1,3,5,4,2),]
#设置新行名
row.names(sample1)<-sample1[,1]
sample1<-data.frame(t(sample1[,-1]))

write.table(sample1, file = "stage.mean.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#安装R包
BiocManager::install("Mfuzz")
library("Mfuzz")
library(Mfuzz)
#构建对象
sample1<-as.matrix(sample1)
sample1<- ExpressionSet(assayData = sample1)

#处理缺失值和异常值
sample1 <- filter.NA(sample1, thres = 0.25)#排除超过25%的测量缺失的基因
sample1 <- fill.NA(sample1, mode = 'mean')
sample1 <- filter.std(sample1, min.std = 0)
#标准化
sample1 <- standardise(sample1)

#设置随机种子，设置需要展示的cluster的数量，然后聚类
set.seed(123)
cluster_num <- 6
sample1_cluster <- mfuzz(sample1, c = cluster_num, m = mestimate(sample1))

#作图
mfuzz.plot2(sample1, cl = sample1_cluster, mfrow = c(2, 3),
            time.labels = colnames(sample1),centre=TRUE,x11=F)
#导出基因
dir.create(path="mfuzz",recursive = TRUE)
for(i in 1:6){
  potname<-names(sample1_cluster$cluster[unname(sample1_cluster$cluster)==i])
  write.csv(sample1_cluster[[4]][potname,i],paste0("mfuzz","/mfuzz_",i,".csv"))
}

####脂肪酸代谢相关####

exp0 <- read.table("KIRC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE53757exp <- read.table("GSE53757.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE73731exp <- read.table("GSE73731exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE126694exp <- read.table("GSE126694exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE215286exp <- read.table("GSE215286exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


group <- read.table("group2-2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

library(tidyverse)
library(GSEABase)
library(GSVA)
library(pheatmap)

exp <- read.table("LUAD_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)



geneSet <- read.csv("PDZ protein.gmt",header = F,sep = "\t")
geneSet <- read.csv("genesets.v2023.2.Hs.gmt",header = F,sep = "\t")




exp <- exp0
exp <- GSE53757exp
exp <- GSE73731exp
exp <- GSE126694exp
exp <- GSE215286exp

exp <- as.matrix(exp)
gs.exp <- gsva(exp, l, kcdf = "Gaussian", min.sz = 10)

write.table(gs.exp, file = "TCGA.GSVA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(gs.exp, file = "GSE53757.GSVA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(gs.exp, file = "GSE73731.GSVA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(gs.exp, file = "GSE126694.GSVA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(gs.exp, file = "GSE215286.GSVA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

group <- read.table("GSE53757group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("GSE126694group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

gs.exp <- gs.exp[,rownames(group)]

PDZK1 <- as.data.frame(exp["PDZK1",])

PDZK1 <- as.data.frame(exp["PDZK1",rownames(group)])
colnames(PDZK1) <- "PDZK1"

gs.exp <- as.data.frame(t(gs.exp))

data <- cbind(gs.exp,PDZK1)
#data <- data[,rownames(Rfile)]

Rfile <- as.data.frame(t(gs.exp[1,]))
Pfile <- as.data.frame(t(gs.exp[1,]))



for (i in 1:78){
  
  c_r=cor(as.numeric(data[,i]),as.numeric(data[,"PDZK1"]),method="pearson",use = 'pairwise.complete.obs')
  p=cor.test(as.numeric(data[,i]),as.numeric(data[,"PDZK1"]),method ="pearson")[[3]]
  ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
  Rfile[i,1] <- c_r
  Pfile[i,1] <- p
  
  
}

colnames(Rfile) <- c("TCGA.KIRC","GSE53757","GSE73731")
colnames(Pfile) <- c("TCGA.KIRC","GSE53757","GSE73731")

write.table(Pfile, file = "GSE215286.Pfile.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Rfile, file = "GSE215286.Rfile.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#GSEA
library(limma) 

exp <- read.table("exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


group1 <- read.table("N_group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

group2 <- read.table("T_group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

group2 <- group[group$group=="tumor",]
data0 <- exp
group <- exp["PDZK1",]
group <- as.data.frame(t(group))

group <- group[rownames(group2),,drop = F]

group2$PDZK1 <- group$PDZK1

write.table(group2, file = "GSE76427group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


data0 <- data0[,rownames(group)]
N_data <- data0[,rownames(group1)]
T_data <- data0[,rownames(group2)]

group <- group2
group$PDZK1 <- t(data["PDZK1",])
colnames(group$PDZK1) <- "PDZK1"
group$group1 <- ifelse(group$PDZK1 > median(group$PDZK1),"high","low")

data<-T_data

design <- model.matrix(~0+factor(group$group1))


rownames(design) <- rownames(group)

colnames(design)=c("high","low")


#design数据表要符合都是数字形式0、1，data看情况要转置
contrast.matrix<-makeContrasts("high-low",levels=design)

##step1
fit <- lmFit(data,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "high-low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)






library(fgsea)
library(biomaRt)
library(enrichplot)
library(tidyverse)
library(org.Hs.eg.db) 
library(clusterProfiler)
# 可视化
library(enrichplot)
library(ggplot2)
library(GseaVis)



mydata <- read.table("TCGA.high-low_tumor.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mydata <- read.table("GSE53757high-low_tumor.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mydata <- read.table("GSE73731high-low.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mydata <- read.table("GSE126694.high-low.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mydata <- read.table("GSE215286.high-low.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)




mydata <- read.table("high-low.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mydata <- nrDEG

data_sort <- mydata %>%
  arrange(desc(logFC))

gene_list <- data_sort$logFC
names(gene_list) <- rownames(data_sort)
head(gene_list)


gmt <- read.gmt("genesets.v2023.2.Hs.gmt")

?GSEA
res <- GSEA(
  gene_list,
  pvalueCutoff=1,
  TERM2GENE = gmt
)

gseaplot2(res, title = res$Description[2], geneSetID = 2,pvalue_table = T)


gseaNb(object = res,
       geneSetID = 'KEGG_P53_SIGNALING_PATHWAY')


gseaNb(res,
       geneSetID = 'KEGG_P53_SIGNALING_PATHWAY',
       addPval = T,
       pCol = "steelblue",
       pvalX = 0.65, # 位置
       pvalY = 0.7,
       pHjust = 0, # 对齐方式
       nesDigit = 4, # 小数点位数
       pDigit = 4
)




save(res,file = "TCGA.res_tumor.Rdata")
save(res,file = "GSE53757.res_tumor.Rdata")
save(res,file = "GSE73731.res_tumor.Rdata")
save(res,file = "GSE126694.res_tumor.Rdata")
save(res,file = "GSE215286.res_tumor.Rdata")



write.table(res, file = "TCGA.result_tumor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(res, file = "GSE53757.result_tumor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(res, file = "GSE73731.result_tumor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(res, file = "GSE126694.result_tumor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(res, file = "GSE215286.result_tumor.txt",sep = "\t",row.names = T,col.names = NA,quote = F)




#多条通路一起展示
terms <- c("KEGG_FATTY_ACID_METABOLISM","KEGG_PPAR_SIGNALING_PATHWAY","HALLMARK_FATTY_ACID_METABOLISM")

terms <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING","KEGG_WNT_SIGNALING_PATHWAY","KEGG_MTOR_SIGNALING_PATHWAY",
           "KEGG_P53_SIGNALING_PATHWAY")

terms <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING","KEGG_WNT_SIGNALING_PATHWAY")

terms <- c("KEGG_MTOR_SIGNALING_PATHWAY","KEGG_P53_SIGNALING_PATHWAY")

terms <- c("KEGG_JAK_STAT_SIGNALING_PATHWAY","HALLMARK_KRAS_SIGNALING_DN")

terms <- c("GOBP_FATTY_ACID_CATABOLIC_PROCESS","GOBP_LIPID_MODIFICATION")


terms <- c("GOBP_FATTY_ACID_BETA_OXIDATION","GOBP_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_DEHYDROGENASE")


gseaNb(object = res,
       geneSetID = terms,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.7,0.8),
       addPval = T,
       pvalX = 0.05,pvalY = 0.05,
       rmPrefix = F)

library(corrplot)
res1 <- read.table("TCGA.result_tumor.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
res2 <- read.table("GSE53757.result_tumor.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
res3 <- read.table("GSE73731.result_tumor.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
res4 <- read.table("GSE126694.result_tumor.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
res5 <- read.table("GSE215286.result_tumor.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

res2 <- res2[rownames(res1),]
res3 <- res3[rownames(res1),]
res4 <- res4[rownames(res1),]
res5 <- res5[rownames(res1),]

data <- res1[,"NES",drop=F]
data$GSE53757 <- res2$NES
data$GSE73731 <- res3$NES
data$GSE126694 <- res4$NES
data$GSE215286 <- res5$NES


pmt <- res1[,"pvalue",drop=F]
pmt$GSE53757 <- res2$pvalue
pmt$GSE73731 <- res3$pvalue
pmt$GSE126694 <- res4$pvalue
pmt$GSE215286 <- res5$pvalue

colnames(data) <- c("TCGA.KIRC","GSE53757","GSE73731","GSE126694","GSE215286")
colnames(pmt) <- c("TCGA.KIRC","GSE53757","GSE73731","GSE126694","GSE215286")

data1 <- data[c(1:37),]
data2 <- data[c(38:75),]
pmt1 <- pmt[c(1:37),]
pmt2 <- pmt[c(38:75),]

data1 <- as.matrix(t(data1))
data2 <- as.matrix(t(data2))
pmt1 <- as.matrix(t(pmt1))
pmt2 <- as.matrix(t(pmt2))

P1 <- corrplot(data1,
               method = "circle",           # 图案形状 "square"方框,"circle"圆, "ellipse"椭圆, "number"数字, "shade"阴影花纹, "color"颜色方框, "pie饼图"
               type = "full",               # 绘制范围"full"全部, "lower"下半部分, "upper"半部分
               col=colorRampPalette(c("#839EDB", "white", "#FF8D8D"))(100), # 主体颜色
               bg = "white",                # 背景颜色
               col.lim = c(-3,3),         # 数据颜色的范围，是相关性数据的话，直接is.corr = T就好
               title = "",                  # 标题
               is.corr = F,                 # 输入的矩阵是否是相关性矩阵，如果是的话，数据范围会限制到-1到1
               add = F,                     # 是否在原来的图层上添加图形
               diag = T,         # 下 左 上 右 边距
               addgrid.col = NULL,          # 网格线的颜色，NA为不绘制，NULl为默认的灰色
               addCoef.col = "grey20",          # 当method!="number"时，是否显示相关性数值，显示的颜色
               addCoefasPercent = F,        # 是否把相关性数值改为百分数
               order = "original",          # 排序方式 c("original", "AOE", "FPC", "hclust", "alphabet"), original：原始状态，alphabet：字母顺序 hclust，分层聚类顺序
               hclust.method = c("complete", "ward", "ward.D", "ward.D2", "single", "average","mcquitty", "median", "centroid"), # 当order = "hclust"时，分层聚类的算法
               tl.pos = "lt",               # 坐标轴标签的位置'lt', 'ld', 'td', 'd' or 'n'   # 左边 d中间
               tl.cex = 0.5,                  # 坐标轴标签字体的大小
               tl.col = "black",            # 坐标轴标签字体的颜色
               tl.offset = 0.5,             # 坐标轴标签离图案的距离
               tl.srt = 45,                 # 坐标轴标签旋转角度
               cl.pos = "r",                # 图例位置：r：右边 b：下边 n：不显示
               cl.length = NULL,            # 数字越大，图例的分隔越稠
               cl.cex = 0.8,                # 图例的字体大小
               cl.ratio = 0.15,             # 图例的宽度
               cl.align.text = "l",         # 图例文字的对齐方式 l左对齐 c居中 r右对齐
               cl.offset = 0.5,             # 图例文字距离图例颜色条的距离 居中时无效
               number.cex = 0.7,              # 相关性数字标签的字体大小
               number.font = 2,             # 相关性数字标签的字体
               number.digits = 2,           # 相关性数字标签，保留的小数点位数
               na.label = "",               # 当为NA时，显示的内容
               p.mat = pmt1,           # P值矩阵
               sig.level = 0.05,            # 当p大于sig.level时触发动作
               insig = "pch",               # p值大于sig.level时的方案"pch"图案, "p-value"P值数字, "blank"空白, "n"无操作, "label_sig"星号,
               pch = 4,                     # 当insig = "pch"时的图案形状 4为叉
               pch.col = "red",           # 图案颜色
               pch.cex = 1.5,                 # 图案大小
               plotCI = "n",                # c("n", "square", "circle", "rect"),  # p值置信区间的方案
               lowCI.mat = testRes$lowCI,   # p值置信区间下边界数据
               uppCI.mat = testRes$uppCI,   # p值置信区间上边界数据
)


P2 <- corrplot(data2,
               method = "circle",           # 图案形状 "square"方框,"circle"圆, "ellipse"椭圆, "number"数字, "shade"阴影花纹, "color"颜色方框, "pie饼图"
               type = "full",               # 绘制范围"full"全部, "lower"下半部分, "upper"半部分
               col=colorRampPalette(c("#839EDB", "white", "#FF8D8D"))(100), # 主体颜色
               bg = "white",                # 背景颜色
               col.lim = c(-3,3),         # 数据颜色的范围，是相关性数据的话，直接is.corr = T就好
               title = "",                  # 标题
               is.corr = F,                 # 输入的矩阵是否是相关性矩阵，如果是的话，数据范围会限制到-1到1
               add = F,                     # 是否在原来的图层上添加图形
               diag = T,         # 下 左 上 右 边距
               addgrid.col = NULL,          # 网格线的颜色，NA为不绘制，NULl为默认的灰色
               addCoef.col = "grey20",          # 当method!="number"时，是否显示相关性数值，显示的颜色
               addCoefasPercent = F,        # 是否把相关性数值改为百分数
               order = "original",          # 排序方式 c("original", "AOE", "FPC", "hclust", "alphabet"), original：原始状态，alphabet：字母顺序 hclust，分层聚类顺序
               hclust.method = c("complete", "ward", "ward.D", "ward.D2", "single", "average","mcquitty", "median", "centroid"), # 当order = "hclust"时，分层聚类的算法
               tl.pos = "lt",               # 坐标轴标签的位置'lt', 'ld', 'td', 'd' or 'n'   # 左边 d中间
               tl.cex = 0.5,                  # 坐标轴标签字体的大小
               tl.col = "black",            # 坐标轴标签字体的颜色
               tl.offset = 0.5,             # 坐标轴标签离图案的距离
               tl.srt = 45,                 # 坐标轴标签旋转角度
               cl.pos = "r",                # 图例位置：r：右边 b：下边 n：不显示
               cl.length = NULL,            # 数字越大，图例的分隔越稠
               cl.cex = 0.8,                # 图例的字体大小
               cl.ratio = 0.15,             # 图例的宽度
               cl.align.text = "l",         # 图例文字的对齐方式 l左对齐 c居中 r右对齐
               cl.offset = 0.5,             # 图例文字距离图例颜色条的距离 居中时无效
               number.cex = 0.7,              # 相关性数字标签的字体大小
               number.font = 2,             # 相关性数字标签的字体
               number.digits = 2,           # 相关性数字标签，保留的小数点位数
               na.label = "",               # 当为NA时，显示的内容
               p.mat = pmt2,           # P值矩阵
               sig.level = 0.05,            # 当p大于sig.level时触发动作
               insig = "pch",               # p值大于sig.level时的方案"pch"图案, "p-value"P值数字, "blank"空白, "n"无操作, "label_sig"星号,
               pch = 4,                     # 当insig = "pch"时的图案形状 4为叉
               pch.col = "red",           # 图案颜色
               pch.cex = 1.5,                 # 图案大小
               plotCI = "n",                # c("n", "square", "circle", "rect"),  # p值置信区间的方案
               lowCI.mat = testRes$lowCI,   # p值置信区间下边界数据
               uppCI.mat = testRes$uppCI,   # p值置信区间上边界数据
)
library(ggplot2)

library(cowplot)
plot<-plot_grid(P1,P2,ncol=1,align = 'v')
plot

#不同stage比较
GSVA.exp <- read.table("TCGA.GSVA.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("TCGA_group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


GSVA.exp <- read.table("GSE73731.GSVA.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("GSE73731group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


data <- GSVA.exp
data <- data[,rownames(group)]
data <- as.data.frame(t(data))

data <- cbind(data,group)
data$stage <- factor(data$stage,levels=c("stage i","stage ii","stage iii","stage iv"))
data$stage <- factor(data$stage,levels=c("Stage: 1","Stage: 2","Stage: 3","Stage: 4"))




library(ggpubr) 
p <- ggboxplot(data, x = "stage", y = "GOBP_FATTY_ACID_BETA_OXIDATION",
               color = "stage", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0","#2b83ba", '#eb4b3a',"red","grey"))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
#更换颜色

p + stat_compare_means()


p2 <- ggboxplot(data, x = "group1", y = "TNFRSF12A",
                color = "group1", palette = "jco",
                add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0", '#eb4b3a'))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
#更换颜色



#成组t检验


t_table <- function(data, dvs, iv,
                    var_equal = TRUE,
                    p_adj = "none",
                    alpha = 0.05,
                    paired = FALSE,
                    wilcoxon = FALSE) {
  if (!inherits(data, "data.frame")) {
    stop("data must be a data.frame")
  }  
  if (!all(c(dvs, iv) %in% names(data))) {
    stop("at least one column given in dvs and iv are not in the data")
  }  
  if (!all(sapply(data[, dvs], is.numeric))) {
    stop("all dvs must be numeric")
  }  
  if (length(unique(na.omit(data[[iv]]))) != 2) {
    stop("independent variable must only have two unique values")
  }  
    out <- lapply(dvs, function(x) {
    if (paired == FALSE & wilcoxon == FALSE) {
      tres <- t.test(data[[x]] ~ data[[iv]], var.equal = var_equal)
    }    
      else if (paired == FALSE & wilcoxon == TRUE) {
      tres <- wilcox.test(data[[x]] ~ data[[iv]])
    }
      else if (paired == TRUE & wilcoxon == FALSE) {
      tres <- t.test(data[[x]] ~ data[[iv]],
        var.equal = var_equal,
        paired = TRUE
      )
    }    else {
      tres <- wilcox.test(data[[x]] ~ data[[iv]],
        paired = TRUE
      )
    }
    c(
      p_value = tres$p.value
    )
  })  
  out <- as.data.frame(do.call(rbind, out))
  out <- cbind(variable = dvs, out)
  names(out) <- gsub("[^0-9A-Za-z_]", "", names(out))
  out$p_value <- ifelse(out$p_value < 0.001,
    "<0.001",
    round(p.adjust(out$p_value, p_adj), 3)
  )
  out$conclusion <- ifelse(out$p_value < alpha,
    paste0("Reject H0 at ", alpha * 100, "%"),
    paste0("Do not reject H0 at ", alpha * 100, "%")
  )  
return(out)
}


t_table <- function(data, dvs, iv,
                    var_equal = TRUE,
                    p_adj = "none",
                    alpha = 0.05,
                    paired = FALSE,
                    wilcoxon = FALSE) {
  if (!inherits(data, "data.frame")) {
    stop("data must be a data.frame")
  }  if (!all(c(dvs, iv) %in% names(data))) {
    stop("at least one column given in dvs and iv are not in the data")
  }  if (!all(sapply(data[, dvs], is.numeric))) {
    stop("all dvs must be numeric")
  }  if (length(unique(na.omit(data[[iv]]))) != 2) {
    stop("independent variable must only have two unique values")
  }  
  out <- lapply(dvs, function(x) {
    if (paired == FALSE & wilcoxon == FALSE) {
      tres <- t.test(data[[x]] ~ data[[iv]], var.equal = var_equal)
    }    
    else if (paired == FALSE & wilcoxon == TRUE) {
      tres <- wilcox.test(data[[x]] ~ data[[iv]])
    }
    else if (paired == TRUE & wilcoxon == FALSE) {
      tres <- t.test(data[[x]] ~ data[[iv]],
                     var.equal = var_equal,
                     paired = TRUE
      )
    }    else {
      tres <- wilcox.test(data[[x]] ~ data[[iv]],
                          paired = TRUE
      )
    }
    c(
      p_value = tres$p.value
    )
  })  
  out <- as.data.frame(do.call(rbind, out))
  out <- cbind(variable = dvs, out)
  names(out) <- gsub("[^0-9A-Za-z_]", "", names(out))
  out$p_value <- ifelse(out$p_value < 0.001,
                        "<0.001",
                        round(p.adjust(out$p_value, p_adj), 3)
  )
  out$conclusion <- ifelse(out$p_value < alpha,
                           paste0("Reject H0 at ", alpha * 100, "%"),
                           paste0("Do not reject H0 at ", alpha * 100, "%")
  )  
  return(out)
}


result <- t_table(
  data = data,
  c("GOBP_FATTY_ACID_BETA_OXIDATION","GOBP_FATTY_ACID_CATABOLIC_PROCESS"),
  "stage")
result


length(unique(na.omit(data[["stage"]])))
t.test(data[["GOBP_FATTY_ACID_BETA_OXIDATION"]] ~ data[["stage"]], var.equal = var_equal)


library(tidyverse)
library(rstatix)
library(ggpubr)
set.seed(123)
data("anxiety", package = "datarium")

anx

anxiety %>% sample_n_by(group, size = 1)
anxiety <- anxiety %>%
  gather(key = "time", value = "score", t1, t2, t3) %>%
  convert_as_factor(id, time)
set.seed(123)
anxiety %>% sample_n_by(group, time, size = 1)
stat.test <- anxiety %>%
  group_by(group) %>%
  pairwise_t_test(
    score ~ time, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-df, -statistic, -p) # Remove details
stat.test

bxp <- ggboxplot(
  anxiety, x = "group", y = "score",
  color = "time", palette = "jco"
)

stat.test <- stat.test %>% add_xy_position(x = "stage")
bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08
)

bxp + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = TRUE, tip.length = 0
)


class(data$stage)
data0 <- data

data <- data0[data0$cancer_status!="0",]
data$cancer_status <- factor(data$cancer_status,levels = c("TUMOR FREE","WITH TUMOR"))

stat.test <- data %>%
  pairwise_t_test(
    HALLMARK_FATTY_ACID_METABOLISM ~ stage, paired = FALSE, 
    p.adjust.method = "bonferroni"
  ) 

stat.test
stat.test <- stat.test %>% add_xy_position(x = "stage")

p <- ggboxplot(data, x = "stage", y = "HALLMARK_FATTY_ACID_METABOLISM",
               color = "stage", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0","#2b83ba", '#eb4b3a',"red","grey"))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

p+stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = F, tip.length = 0
)+stat_compare_means(label.y = 1.2)


p+stat_compare_means(label.y = 1.2)


n_fun <- function(x){
  return(data.frame(y = median(x), label = paste0("n = ",length(x))))
}

#HALLMARK_FATTY_ACID_METABOLISM
#KEGG_FATTY_ACID_METABOLISM
#KEGG_PPAR_SIGNALING_PATHWAY
#HALLMARK_IL6_JAK_STAT3_SIGNALING
#KEGG_P53_SIGNALING_PATHWAY
#KEGG_WNT_SIGNALING_PATHWAY
data$laterality
data <- data[data$gender=="male",]
p <- ggboxplot(data, x = "laterality", y = "KEGG_FATTY_ACID_METABOLISM",
               color = "laterality", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0","#2b83ba", '#eb4b3a',"red","grey"))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))


p + stat_compare_means()+
  stat_summary(fun.data = n_fun, geom = "text")


####蛋白质表达谱数据####
#数据读取

group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp0 <- read.csv("初始蛋白表达矩阵.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


gene <- read.csv("gene_mTOR.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- gene$gene

PDZK1 <- exp[c("PDZK1","ULK1"),]

PDZK1 <- as.data.frame(t(PDZK1))

PDZK1 <- PDZK1[-1,]

PDZK1 <- cbind(PDZK1,group)
PDL1 <- exp0[exp0$`Gene Symbol`=="MTOR",]
###
data <- exp[gene,]
data <- na.omit(data)


data <- as.data.frame(t(data))
colnames(data) <- "MTOR"
data <- data[rownames(group),]
data <- cbind(data,group)
data <- t(PDZK1)



data2 <- na.omit(data)
data <- exp[gene,]
data <- data[,rownames(PDZK1)]

data2 <- data2[,rownames(PDZK1)]
library(pheatmap)

bk <- c(seq(-5,-0.1,by=0.1),seq(0,5,by=0.1))

pheatmap(data2,
         annotation_col=PDZK1,
         scale = "row",
         show_rownames = T,
         show_colnames =T,
         color = colorRampPalette(c("navy", "white", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6,
         breaks = bk)





exp <- exp0[!duplicated(exp0$`Gene Symbol`),]
rownames(exp) <- exp$`Gene Symbol`
exp <- exp[,-1]


library("dplyr") 
library(tidyverse)
PDZK1 <- PDZK1 %>% arrange(PDZK1)

PDZK1 <- PDZK1 %>% arrange(group)



PDZK1 <- as.data.frame(t(PDZK1))
class(PDZK1$PDZK1)
PDZK1$PDZK1 <- sort(PDZK1$PDZK1)

PDZK1 <- PDZK1[rownames(group),,drop=F]
colnames(PDZK1) <- "PDZK1"
PDZK1 <- cbind(PDZK1,group)
group2$PDZK1 <- as.numeric(group2$PDZK1)

group2 <- PDZK1[PDZK1$group=="T",]

group2$type <- ifelse(group2$PDZK1 > median(group2$PDZK1), "high", "low")
exp2 <- exp[,rownames(group2)]

median(group$PDZK1)
#配对t检验

PDZK1$PDZK1 <- as.numeric(PDZK1$PDZK1)
class(PDZK1$PDZK1)

t.test (PDZK1~group, data = PDZK1, paired = TRUE)
PDZK1$sample <- c(1:232,1:232)
library(tidyverse)
library(ggplot2)
library(forcats)
library(rstatix)
library(ggpubr)


ggplot(PDZK1, aes(x = group, y = PDZK1)) +
  geom_boxplot(aes(fill = group), show.legend = F, width = 0.6) +  #箱线图
  scale_fill_manual(values = c('#00AFBB', '#E7B800')) +  #设置颜色
  geom_point(size = 2,color='red') + 
  geom_point(size = 3,shape=21) +  #绘制散点
  geom_line(aes(group = sample), color = 'gray70', lwd = 0.2) +  #配对样本间连线
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5), 
        plot.subtitle = element_text(size = 15, hjust = 0.5), 
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black')) +
  labs(x = '', y = 'PDZK1 expression', title = 'T vs TA')+
  stat_compare_means(method = "t.test",paired = TRUE, comparisons=list(c("TA", "T")))#配对t检验


data <- cbind(PDZK1,group)
p <- ggboxplot(data, x = "PDZK1", y = "ULK1",
               color = "ULK1", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#2b83ba", '#eb4b3a',"red","grey"))+ 
  theme_bw() + 
  geom_line(aes(group = ULK1), color = 'gray70', lwd = 0.2) +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))


p + stat_compare_means(method = "t.test",paired = TRUE, comparisons=list(c("TA", "T")))+
  stat_summary(fun.data = n_fun, geom = "text")



#差异分析
library(limma) 


####差异分析
group_list
design <- model.matrix(~0+factor(group2$type))
row.names(design) <- rownames(group2)

colnames(design)=c('high','low')
contrast.matrix<-makeContrasts("high-low",levels=design)

##step1
fit <- lmFit(exp2,design)#data的列名要等于design的行名
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)

write.table(nrDEG, file = "PDZK1.high-low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#GSVA

library(GSEABase)
library(GSVA)

exp <- as.matrix(exp)
gs.exp <- gsva(exp, l, kcdf = "Gaussian", min.sz = 10)

write.table(gs.exp, file = "GSVA.exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


gs.exp <- as.data.frame(t(gs.exp))

gs.exp <- gs.exp[rownames(group),]
data <- cbind(gs.exp,group)


data$PDZK1 <- as.numeric(data$PDZK1)
data$group <- factor(data$group,levels = c("TA","T"))

library(ggstatsplot)
library(ggExtra)
library(ggsci)

#KEGG_P53_SIGNALING_PATHWAY
#KEGG_JAK_STAT_SIGNALING_PATHWAY
#KEGG_MTOR_SIGNALING_PATHWAY
#KEGG_FATTY_ACID_METABOLISM

data <- PDZK1
data$group <- data$group_list

data <- cbind(group,PDZK1)
b <- ggplot(data, aes(x = ULK1, y = PDZK1))
b <- b + geom_point(aes(color = group, shape = group))+
  geom_smooth(aes(color = group, fill = group), method = "lm") +
  geom_rug(aes(color =group)) +
  scale_color_manual(values = c("#2b83ba", "#FC4E07"))+
  scale_fill_manual(values = c("#2b83ba", "#FC4E07"))+ 
  
  ggpubr::stat_cor(aes(color = group))

ggMarginal(
  b,
  type = 'boxplot',
  margins = 'both',
  size = 10,
  groupColour = TRUE,
  groupFill = TRUE
)


data2 <- data[data$group=="T",]
data3 <- data[data$group=="TA",]

Rfile <- as.data.frame(t(data[1,,drop=F]))
Pfile <- as.data.frame(t(data[1,,drop=F]))




for (i in 1:195){
  
  c_r=cor(as.numeric(data3[,i]),as.numeric(data3[,196]),method="pearson",use = 'pairwise.complete.obs')
  p=cor.test(as.numeric(data3[,i]),as.numeric(data3[,196]),method ="pearson")[[3]]
  ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
  Rfile[i,2] <- c_r
  Pfile[i,2] <- p
  
  
}
colnames(Rfile) <- c("T","TA")
colnames(Pfile) <- c("T","TA")
Rfile <- Rfile[-c(196,197),]
Pfile <- Pfile[-c(196,197),]

write.table(Rfile, file = "Rfile.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Pfile, file = "Pfile.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

data <- dataExpr
data <- data[rownames(group),]
data <- cbind(data,group)


b <- ggplot(data, aes(x = PBRM1, y = PDZK1))
b <- b + geom_point(aes(color = group, shape = group))+
  geom_smooth(aes(color = group, fill = group), method = "lm") +
  geom_rug(aes(color =group)) +
  scale_color_manual(values = c("#2b83ba", "#FC4E07"))+
  scale_fill_manual(values = c("#2b83ba", "#FC4E07"))+ 
  
  ggpubr::stat_cor(aes(color = group))

ggMarginal(
  b,
  type = 'boxplot',
  margins = 'both',
  size = 10,
  groupColour = TRUE,
  groupFill = TRUE
)


####mTOR信号通路分析####
#基因集制作
library(tidyverse)
geneSet <- read.csv("genesets.v2023.2.Hs.gmt",sep = "\t",check.names = F,stringsAsFactors = F,header = F,col.names = paste("V", 1:1000,sep = ""))

geneSet <- read.csv("genesets.v2023.2.Hs.gmt",header = F,sep = "\t")


read.table("test.txt", fill = T, col.names = paste("V", 1:1407, sep = ""))

data <- geneSet

data <- t(data)


geneSet <- geneSet[,-2]


geneSet <- geneSet %>%
  column_to_rownames("V1")%>%t()
a <- geneSet
a <- a[1:nrow(a),,drop=F]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}

save(l,file = "mTOR相关基因集.Rdata")


library(GSVA)

exp <- read.table("KIRC_fpkm_mRNA_all.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gs.exp <- read.table("gsva.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

exp <- as.matrix(exp)
gs.exp <- gsva(exp, l, kcdf = "Gaussian", min.sz = 10)
write.table(gs.exp, file = "GSVA.exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)




gs.exp <- as.data.frame(t(gs.exp))

exp <- as.data.frame(exp)
PDZK1 <- exp["PDZK1",]

PDZK1 <- as.data.frame(t(PDZK1))

PDZK1 <- PDZK1[rownames(group),,drop=F]

PDZK1 <- cbind(PDZK1,group)

gs.exp <- gs.exp[rownames(group),]
data <- cbind(gs.exp,group)
data$group <- factor(data$group,levels = c("TA","T"))

library(ggplot2)
library(ggstatsplot)
library(ggExtra)
library(ggsci)

gs.exp <- gs.exp[rownames(group),]

data <- cbind(gs.exp,group)

#CREIGHTON_AKT1_SIGNALING_VIA_MTOR_DN
#PARENT_MTOR_SIGNALING_DN
#PARENT_MTOR_SIGNALING_UP
#REACTOME_MTOR_SIGNALLING
#KEGG_MTOR_SIGNALING_PATHWAY

data <- as.data.frame(t(data))
data <- cbind(data,PDZK1)

b <- ggplot(data, aes(x = RICTOR, y = PDZK1))
b <- b + geom_point(aes(color = group_list, shape = group_list))+
  geom_smooth(aes(color = group_list, fill = group_list), method = "lm") +
  geom_rug(aes(color =group_list)) +
  scale_color_manual(values = c("#2b83ba", "#FC4E07"))+
  scale_fill_manual(values = c("#2b83ba", "#FC4E07"))+ 
  
  ggpubr::stat_cor(aes(color = group_list))


ggMarginal(
  b,
  type = 'boxplot',
  margins = 'both',
  size = 10,
  groupColour = TRUE,
  groupFill = TRUE
)



data2 <- data[data$group_list=="Tumor",]
data3 <- data[data$group_list=="Normal",]

Rfile <- as.data.frame(t(data[1,,drop=F]))
Pfile <- as.data.frame(t(data[1,,drop=F]))

i=1
c_r=cor(as.numeric(data2[,i]),as.numeric(data2[,32]),method="pearson",use = 'pairwise.complete.obs')


for (i in c(1:52)){
  
  c_r=cor(as.numeric(data3[,i]),as.numeric(data3[,53]),method="pearson",use = 'pairwise.complete.obs')
  p=cor.test(as.numeric(data3[,i]),as.numeric(data3[,53]),method ="pearson")[[3]]
  ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
  Rfile[i,2] <- c_r
  Pfile[i,2] <- p
  
  
}
colnames(Rfile) <- c("Tumor","Normal")
colnames(Pfile) <- c("Tumor","Normal")
Rfile <- Rfile[-c(53,54),]
Pfile <- Pfile[-c(53,54),]

write.table(Rfile, file = "mTORgene_Rfile.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Pfile, file = "mTORgene_Pfile.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

gene <- read.csv("gene_mTOR.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- gene$gene
data <- exp[gene,]

data <- data[,rownames(group)]

bk <- c(seq(-10,-0.2,by=0.2),seq(0,10,by=0.2))

pheatmap(data,
         annotation_col=PDZK1,
         scale = "row",
         show_rownames = T,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=6,
         breaks = bk)


data <- as.data.frame(t(data))

data <- cbind(data,PDZK1)
data$group <- factor(data$group,levels = c("normal","tumor"))

gene <- colnames(data)
for (i in 1:22) {
  

              p <- ggboxplot(data, x = "group", y = gene[i],
                             color = "group", palette = "jco",
                             add = "jitter")+
                scale_color_manual(values = c("#2b83ba", '#eb4b3a',"red","grey"))+ 
                theme_bw() + 
                theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
              
              
              p <- p + stat_compare_means()+
                stat_summary(fun.data = n_fun, geom = "text")
              ggsave(filename = paste0("gene/",gene[i],".png"),plot = last_plot(),width = 4, height = 4, units = "in")
              
}

gene <- colnames(data)

p <- ggboxplot(data, x = "group", y = "KEGG_MTOR_SIGNALING_PATHWAY",
               color = "group", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#2b83ba", '#eb4b3a',"red","grey"))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))


p + stat_compare_means()+
  stat_summary(fun.data = n_fun, geom = "text")

ggsave(filename = paste0("passway","KEGG_MTOR_SIGNALING_PATHWAY",".png"),plot = last_plot(),width = 4, height = 4, units = "in")


for (i in 1:52) {
  

        b <- ggplot(data, aes(x = gene[i], y = PDZK1))
        b <- b + geom_point(aes(color = group_list, shape = group_list))+
          geom_smooth(aes(color = group_list, fill = group_list), method = "lm") +
          geom_rug(aes(color =group_list)) +
          scale_color_manual(values = c("#2b83ba", "#FC4E07"))+
          scale_fill_manual(values = c("#2b83ba", "#FC4E07"))+ 
          
          ggpubr::stat_cor(aes(color = group_list))
        
        
        ggMarginal(
          b,
          type = 'boxplot',
          margins = 'both',
          size = 10,
          groupColour = TRUE,
          groupFill = TRUE
        )
}

b <- ggplot(data, aes(x = gene[1], y = PDZK1))
b <- b + geom_point(aes(color = group_list, shape = group_list))+
  geom_smooth(aes(color = group_list, fill = group_list), method = "lm") +
  geom_rug(aes(color =group_list)) +
  scale_color_manual(values = c("#2b83ba", "#FC4E07"))+
  scale_fill_manual(values = c("#2b83ba", "#FC4E07"))+ 
  
  ggpubr::stat_cor(aes(color = group_list))


ggMarginal(
  b,
  type = 'boxplot',
  margins = 'both',
  size = 10,
  groupColour = TRUE,
  groupFill = TRUE
)

####LIHC验证####

exp <- read.table("exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group <- read.table("group2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

gs.exp <- read.table("GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

data <- as.data.frame(t(gs.exp))

data <- data[rownames(group),]
data <- cbind(data,group)
PDZK1 <- exp["PDZK1",]
PDZK1 <- as.data.frame(t(PDZK1))
group <- cbind(group,PDZK1)
class(group$PDZK1)

data$group_list <- factor(data$group_list,levels = c("Normal","Tumor"))


library(ggplot2)
library(ggstatsplot)
library(ggExtra)
library(ggsci)

#KEGG_P53_SIGNALING_PATHWAY
#KEGG_JAK_STAT_SIGNALING_PATHWAY
#KEGG_MTOR_SIGNALING_PATHWAY
#KEGG_FATTY_ACID_METABOLISM



b <- ggplot(data, aes(x = HALLMARK_IL6_JAK_STAT3_SIGNALING, y = PDZK1))
b <- b + geom_point(aes(color = group_list, shape = group_list))+
  geom_smooth(aes(color = group_list, fill = group_list), method = "lm") +
  geom_rug(aes(color =group_list)) +
  scale_color_manual(values = c("#2b83ba", "#FC4E07"))+
  scale_fill_manual(values = c("#2b83ba", "#FC4E07"))+ 
  
  ggpubr::stat_cor(aes(color = group_list))

ggMarginal(
  b,
  type = 'boxplot',
  margins = 'both',
  size = 10,
  groupColour = TRUE,
  groupFill = TRUE
)

ggsave(filename = paste0("passway","KEGG_P53_SIGNALING_PATHWAY",".pdf"),plot = last_plot(),width = 4, height = 4, units = "in")
ggsave(filename = paste0("PDZK1功能相关/","passway","KEGG_P53_SIGNALING_PATHWAY",".png"),plot = b,width = 4, height = 4, units = "in")



