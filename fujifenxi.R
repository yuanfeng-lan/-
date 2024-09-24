####PDZK1差异分析####
setwd("PDZK1_deg")
install.packages("BiocManager") 
if(!require(DESeq2))BiocManager::install('DESeq2')
library(DESeq2)
library(tidyverse)
counts_01A <- read.table("KIRC_counts_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- read.table("KIRC_fpkm_mRNA_01A.txt", sep = "\t",row.names = 1,check.names = F,header = T)
com <- intersect(colnames(counts_01A),colnames(exp))
exp <- exp[,com]
counts_01A <- counts_01A[,com]
identical(colnames(counts_01A),colnames(exp))
gene <- "PDZK1"#每次运行只改这个基因名
med=median(as.numeric(exp[gene,]))

conditions=data.frame(sample=colnames(exp),
                      group=factor(ifelse(exp[gene,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")
write.table(conditions, file = "PDZK1的高低分组.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)

dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds)


save(res,file="res_deseq2_PDZK1.Rda")
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)

write.table(res_deseq2, file = "high-low.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####PDZK1差异分析结果富集分析####
setwd("PDZK1_fuji")
#install.packages("tidyverse")
#install.packages("BiocManager")
BiocManager::install('clusterProfiler')
BiocManager::install('org.Hs.eg.db')
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)

DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

#GO分析
ego <- enrichGO(gene = gene,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
save(ego,ego_res,file = "GO_PDZK1_DEG.Rdata")

#3. 可视化
##3.1 柱状图
barplot(ego, showCategory = 20,color = "pvalue",col = '#66C3A5')
##3.2 气泡图
dotplot(ego, showCategory = 20)
##3.3 分类展示
barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dotplot(ego,showCategory = 10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')


####PDZK1差异分析结果KEGG富集分析####
setwd("PDZK1_KEGG")
#install.packages("tidyverse")
#install.packages("BiocManager")
BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)

DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

#KEGG分析
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk_res <- kk@result
save(kk,kk_res,file = "KEGG_PDZK1_DEG.Rdata")

load("KEGG_PDZK1_DEG.Rdata")

#柱状图
barplot(kk, showCategory = 20,color = "pvalue")
#气泡图
dotplot(kk, showCategory = 20)

dev.off()


####GSEA####
setwd("GSEA_PDZK1")
#install.packages("tidyverse")
#install.packages("BiocManager")
BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05)

DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))


#msigdb_GMTs <- "GSEA.GMTs"
#msigdb <- "GOBP_REGULATION_OF_LIPID_STORAGE.gmt"    #c2.all.v7.0.entrez.gmt 或 c5.all.v7.0.entrez.gmt
#读取上面指定的gmt文件
msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "c5.all.v7.0.entrez.gmt"    #c2.all.v7.0.entrez.gmt 或 c5.all.v7.0.entrez.gmt

kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

geneList = DEG[,3]
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
#降序
geneList = sort(geneList, decreasing = TRUE)

set.seed(1)
KEGG<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
#转换成数据框
KEGG_result_df <- as.data.frame(KEGG)
write.table(KEGG_result_df,file="GSEA_MSigDb_C5_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
save(KEGG,KEGG_result_df,file = "GSEA_deg_PDZK1.rda")
#单个图绘制
library(enrichplot)
gseaplot2(KEGG,1,color="red")
gseaplot2(KEGG,3,color="red",pvalue_table = T)

#汇总结果
gseaplot2(KEGG, geneSetID = c(1,4,21,23,25,43), subplots = 1:3)
gseaplot2(KEGG, geneSetID = c(1,3), subplots = 1:3)
gseaplot2(KEGG, geneSetID = 1:3, subplots = 1)
gseaplot2(KEGG, geneSetID = 1:10, subplots = 1:3)

dev.off()

####GEO数据库差异分析####
exp <- read.table("exp.txt", sep = "\t",row.names = 1,check.names = F,header = T)
gene <- "PDZK1"#每次运行只改这个基因名
med=median(as.numeric(exp[gene,]))

conditions=data.frame(sample=colnames(exp),
                      group=factor(ifelse(exp[gene,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")
write.table(conditions, file = "PDZK1的高低分组.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####富集弦图####


BiocManager::install("GOplot")

#加载GOplot
library(GOplot)

#加载测试数据
data(EC)
#创建circ对象，EC$david为富集分析结果
#EC$genelist为差异表达分析结果
GOBPdavid <- read.table("GOBPdavid.txt", sep = "\t",row.names = 1,check.names = F,header = T)

genelist <- read.table("genelist.txt", sep = "\t",row.names = 1,check.names = F,header = T)
genes <- read.table("gene.txt", sep = "\t",row.names = 1,check.names = F,header = T)

circ <- circle_dat(GOBPdavid, genelist)
gene <- EC$genes
genelist <- EC$genelist
david <- EC$david

process <- EC$process
process <- GOBPdavid[,2]
process <- process[1:10]
#常见chord对象，EC$genes为需要展示的基因包含FC
#EC$process为需要展示的GO条目
chord <- chord_dat(circ, genes, process)

#创建pdf文件，保存弦图
pdf("chord_demo.pdf",height = 14,width = 13)
#绘制弦图
GOChord(chord,   #chord对象
        space = 0.02,  #右侧色块之间的间距
        gene.order = 'logFC',   #基因展示顺序根据logFC来
        gene.space = 0.25,  #基因名字和色块之间的距离
        gene.size = 5  #基因名字大小
)
dev.off()
