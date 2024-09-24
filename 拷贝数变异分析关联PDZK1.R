####KIRC变异景观####

#2024.7.4

library(maftools)
laml.maf = system.file('data_mutations.txt', 'kirc_tcga_pan_can_atlas_2018.tar.gz', package = 'maftools') 
laml.clin = system.file('data_clinical_patient', 'kirc_tcga_pan_can_atlas_2018.tar.gz', package = 'maftools') 

PDZK1 <- read.table("PDZK1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- PDZK1[PDZK1$group=="OV",]

group <- PDZK1
group$group <- ifelse(group$PDZK1>=median(group$PDZK1),"high","low")


comparison <- compareGroups(laml, groups = group)

laml = read.maf(maf = "data_mutations.maf")

data <- laml@data

plotmafSummary(maf=laml, rmOutlier=TRUE, addStat="median", dashboard=TRUE, titvRaw = FALSE)

oncoplot(maf = laml, draw_titv = TRUE,borderCol=NULL)

laml.titv <- titv(maf=laml, plot=FALSE, useSyn=TRUE)
plotTiTv(res=laml.titv)

lollipopPlot(maf=laml, gene="VHL", AACol="HGVSp_Short", showMutationRate=TRUE)

lollipopPlot(maf=laml, gene="PBRM1", AACol="HGVSp_Short", showMutationRate=TRUE)
geneCloud(input=laml, top=20)

plotVaf(maf=laml)

#分组瀑布图
# 从临床数据中提取性别对应的"Tumor_Sample_Barcode"
clin <- read.table("PDZK1.txt", header=T, sep="\t")
clin <- group
group$Tumor_Sample_Barcode <- rownames(group)
clin.high <- subset(clin, group=="high")$Tumor_Sample_Barcode
clin.low <- subset(clin, group=="low")$Tumor_Sample_Barcode
# 使用subsetMaf构建男性和女性的MAF对象
luad.high <- subsetMaf(maf=laml, tsb=clin.high, isTCGA=F)

luad.low <- subsetMaf(maf=laml, tsb=clin.low, isTCGA=F)

# 使用mafCompare比较差异突变基因
fvsm <- mafCompare(m1=luad.high, m2=luad.low, m1Name="high", m2Name="low", minMut=5)
# 结果保存到文件"female_vs_male.tsv"
write.table(fvsm$results, file="low_vs_high.tsv", quote=FALSE, row.names=FALSE, sep="\t")
results <- fvsm$results
genes <- results[1:20,]
genes <- genes$Hugo_Symbol

coOncoplot(m1=luad.high, m2=luad.low, m1Name="high", m2Name="low",genes = genes)

genes <- c("TP53", "TTN", "FLG2", "MUC16", "CSMD3", "HMCN1", "USH2A", "AHNAK")

genes <- c("VHL","PBRM1","TTN","SETD2","BAP1","MTOR","MUC16","KDM5C","DNAH9","CSMD3")












#生存分析
mafSurvival(maf = laml, genes = 'VHL', time = 'OS_MONTHS', Status = 'OS_STATUS', isTCGA = TRUE)

prog_geneset = survGroup(maf = laml, top = 20, geneSetSize = 2, time = "OS_MONTHS", Status = "OS_STATUS", verbose = FALSE)

print(prog_geneset)

write.table(prog_geneset, file = "prog_geneset.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


mafSurvGroup(maf = laml, geneSet = c("VHL", "MTOR"), time = "OS_MONTHS", Status = "OS_STATUS")
mafSurvGroup(maf = laml, geneSet = c("VHL", "BAP1"), time = "OS_MONTHS", Status = "OS_STATUS")
mafSurvGroup(maf = laml, geneSet = c("VHL"), time = "OS_MONTHS", Status = "OS_STATUS")

#致癌通路
OncogenicPathways(maf = laml)

OncogenicPathways(maf = laml)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("maftools")

#计算TMB
tmb_table_wt_log = tmb(maf = laml)
#不取log
tmb_table_wo=tmb(maf = laml,logScale = F)

write.table(tmb_table_wt_log, file = "TMB.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

exp <- read.table("data_mrna_seq_v2_rsem.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)

exp <- na.omit(exp)

exp <- aggregate(exp[,3:512],by=list(exp$Hugo_Symbol),mean,na.rm= TRUE)#排除最后的一列group
exp <- exp[-1,]
rownames(exp) <- exp$Group.1
exp <- exp[,-1]

PDZK1 <- exp["PDZK1",]

write.table(PDZK1, file = "PDZK1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

PDZK1 <- as.data.frame(t(PDZK1))

rownames(tmb_table_wo) <- tmb_table_wo$Tumor_Sample_Barcode
PDZK1 <- PDZK1[rownames(tmb_table_wo),,drop=F]
PDZK1 <- log2(PDZK1)
data <- cbind(tmb_table_wo,PDZK1)
?aggregate

write.table(exp, file = "exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


library(ggstatsplot)
library(ggExtra)
library(ggsci)

p <- ggscatterstats(data, 
                    y =total_perMB, 
                    x =PDZK1,
                    type = "pearson",
                    centrality.para = "mean",                              
                    margins = "both",                                         
                    xfill = "#009E73",
                    yfill = "#D55E00", 
                    marginal.type = "boxplot",    #类型可以换成density,boxplot,violin,densigram
                    
                    title = ""
)

p
data$TMB_group <- ifelse(data$total_perMB > median(data$total_perMB),"high","low")

library(ggpubr) 
p <- ggboxplot(data, x = "TMB_group", y = "PDZK1",
               color = "TMB_group", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0","#2b83ba", '#eb4b3a',"red"))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
#更换颜色

p + stat_compare_means()+
  stat_summary(fun.data = n_fun, geom = "text")


####导入gitic文件####
library(tidyverse)
library(maftools)

laml.gistic = readGistic(gisticAllLesionsFile = './Gistic 2.0 results/all_lesions.conf_99.txt',
                         gisticAmpGenesFile = './Gistic 2.0 results/amp_genes.conf_99.txt', 
                         gisticDelGenesFile = './Gistic 2.0 results/del_genes.conf_99.txt', 
                         gisticScoresFile = './Gistic 2.0 results/scores.gistic', 
                         isTCGA = T)


####PDZK1关联VHL和PBRM1突变####
PDZK1 <- read.table("PDZK1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
sur <- read.table("survival.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
sur <- sur[rownames(PDZK1),]



PDZK1$PDZK1 <- log2(PDZK1$PDZK1)
data <- PDZK1
data$VHL <- factor(data$VHL,levels = c("WT","MUTANT"))
data$PBRM1 <- factor(data$PBRM1,levels = c("WT","MUTANT"))


library(ggpubr) 
p <- ggboxplot(data, x = "VHL", y = "PDZK1",
               color = "VHL", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0","#2b83ba", '#eb4b3a',"red"))+ 
  theme_bw() + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
#更换颜色

p + stat_compare_means()+
  stat_summary(fun.data = n_fun, geom = "text")



n_fun <- function(x){
  return(data.frame(y = median(x), label = paste0("n = ",length(x))))
}

ggplot(mtcars, aes(factor(cyl), mpg, label=rownames(mtcars))) +
  geom_boxplot(fill = "grey80", colour = "#3366FF") + 
  stat_summary(fun.data = n_fun, geom = "text")


####生存分析####

data <- cbind(PDZK1,sur)
data <- na.omit(data)

write.table(data, file = "sur.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

data$group1 <- ifelse(data$PDZK1 > median(data$PDZK1),"high","low")


data$group2 <- ifelse(data$group1=="low"&data$VHL=="MUTANT","VHL mutant+PDZK1 low","other group")

data$group5 <- ifelse(data$group1=="low"&data$PBRM1=="MUTANT","PBRM1 mutant+PDZK1 low","other group")


data$group3 <- paste(data$group1,data$VHL,sep = "+")
data$group4 <- paste(data$group1,data$PBRM1,sep = "+")


library(survival)
fitd <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ group5,
                 data      = data,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ group5, data = data)
summary(fit)

p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))


library(survminer)

ggsurvplot(fit,
           data = data,
           conf.int = TRUE,
           pval = p.lab,
           pval.size = 5,
           #legend.labs = c("High","Low"),
           legend.title = '',
           xlab = "Time (Months)",
           ylab = "Overall survival",
           break.time.by = 20,
           risk.table.title = "",
           palette = c("#d7191c","#2b83ba"),
           risk.table = F,
           risk.table.height = .25
)


##多组生存分析

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

fitd <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ group4,
                 
                 data = data)

p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ group4,
               
               data = data,
               
               type = "kaplan-meier",
               
               error = "greenwood",
               
               conf.type = "plain",
               
               na.action = na.exclude)

# 配对生存分析

ps <- pairwise_survdiff(Surv(OS_MONTHS, OS_STATUS) ~ group4,
                        
                        data = data,
                        
                        p.adjust.method = "none") # 这里不使用矫正，若需要矫正可以将none替换为BH

# 设置颜色

mycol <- brewer.pal(n = 10, "Paired")[c(2,4,6,8,10)]

# 绘制基础图形

## 隐藏类标记

names(fit$strata) <- gsub("group=", "", names(fit$strata))

## 生存曲线图

p <- ggsurvplot(fit = fit,
                
                conf.int = FALSE, # 不绘制置信区间
                
                risk.table = TRUE, # 生存风险表
                
                risk.table.col = "strata",
                
                palette = mycol, # KM曲线颜色
                
                data = data,
                
                xlim = c(0,120), # 时间轴，一般考虑5年（原文）或者10年长度
                legend.labs = c("high PDZK1+PBRM1 MU", "high PDZK1+PBRM1 WT","low PDZK1+PBRM1 MU","low PDZK1+PBRM1 WT"),
                
                size = 1,
                
                break.time.by = 12, # 时间轴的刻度（每年）
                
                legend.title = "",
                
                xlab = "Time (months)",
                
                ylab = "Overall survival",
                
                risk.table.y.text = FALSE,
                
                tables.height = 0.3) # 风险表的高度

## 添加overall pvalue

p.lab <- paste0("log-rank test P",
                
                ifelse(p.val < 0.001, " < 0.001", # 若P值<0.001则标记为“<0.001”
                       
                       paste0(" = ",round(p.val, 3))))

p$plot <- p$plot + annotate("text",
                            
                            x = 0, y = 0.55, # 在y=0.55处打印overall p值
                            
                            hjust = 0,
                            
                            fontface = 4,
                            
                            label = p.lab)

## 添加配对表格

addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                         
                                         round(ps$p.value, 3))))

addTab[is.na(addTab)] <- "-"

df <- tibble(x = 0, y = 0, tb = list(addTab))

p$plot <- p$plot +
  
  geom_table(data = df,
             
             aes(x = x, y = y, label = tb),
             
             table.rownames = TRUE)

##生成图片
p




