####临床数据获取####




####重新性状和共表达模块分析####

library(WGCNA)

#

library(reshape2)
library(stringr)

# 
options(stringsAsFactors = FALSE)
# 打开多线程
enableWGCNAThreads()

## Allowing parallel execution with up to 47 working processes.

# 常规表达矩阵，log2转换后或
# Deseq2的varianceStabilizingTransformation转换的数据
# 如果有批次效应，需要事先移除，可使用removeBatchEffect
# 如果有系统偏移(可用boxplot查看基因表达分布是否一致)，
# 需要quantile normalization

exprMat <- "WGCNA/LiverFemaleClean.txt"
exprMat <- 
  
  
  exprMat <- read.table("normalize.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


exp <- read.table("exp_p.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

exprMat <- exprMat[rownames(exp),]

ids <- read.table("ids.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)




DEGs <- read.table("gene.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- DEGs$gene

group <- read.table("phe.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

sample <- subset(group,group!='SS') 
exprMat <- exprMat[,rownames(group)]
# 官方推荐 "signed" 或 "signed hybrid"
# 为与原文档一致，故未修改 
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

##导入数据##
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, 
                       quote="", comment="", check.names=F)
dataExpr <- exprMat
dim(dataExpr)

## [1] 3600  134
dataExpr <- dataExpr[gene,]
dataExpr <- subset(dataExpr,GSM2385720!='NA')


head(dataExpr)[,1:8]

##                 F2_2    F2_3     F2_14    F2_15    F2_19       F2_20
## MMT00000044 -0.01810  0.0642  6.44e-05 -0.05800  0.04830 -0.15197410
## MMT00000046 -0.07730 -0.0297  1.12e-01 -0.05890  0.04430 -0.09380000
## MMT00000051 -0.02260  0.0617 -1.29e-01  0.08710 -0.11500 -0.06502607
## MMT00000076 -0.00924 -0.1450  2.87e-02 -0.04390  0.00425 -0.23610000
## MMT00000080 -0.04870  0.0582 -4.83e-02 -0.03710  0.02510  0.08504274
## MMT00000102  0.17600 -0.1890 -6.50e-02 -0.00846 -0.00574 -0.01807182
##                F2_23    F2_24
## MMT00000044 -0.00129 -0.23600
## MMT00000046  0.09340  0.02690
## MMT00000051  0.00249 -0.10200
## MMT00000076 -0.06900  0.01440
## MMT00000080  0.04450  0.00167
## MMT00000102 -0.12500 -0.06820

## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0)),]

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))


dataExpr <- as.data.frame(t(dataExpr))

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)

##  Flagging genes and samples with too many missing values...
##   ..step 1

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)

## [1]  134 2697

head(dataExpr)[,1:8]

##       MMT00000051 MMT00000080 MMT00000102 MMT00000149 MMT00000159
## F2_2  -0.02260000 -0.04870000  0.17600000  0.07680000 -0.14800000
## F2_3   0.06170000  0.05820000 -0.18900000  0.18600000  0.17700000
## F2_14 -0.12900000 -0.04830000 -0.06500000  0.21400000 -0.13200000
## F2_15  0.08710000 -0.03710000 -0.00846000  0.12000000  0.10700000
## F2_19 -0.11500000  0.02510000 -0.00574000  0.02100000 -0.11900000
## F2_20 -0.06502607  0.08504274 -0.01807182  0.06222751 -0.05497686
##       MMT00000207 MMT00000212 MMT00000241
## F2_2   0.06870000  0.06090000 -0.01770000
## F2_3   0.10100000  0.05570000 -0.03690000
## F2_14  0.10900000  0.19100000 -0.15700000
## F2_15 -0.00858000 -0.12100000  0.06290000
## F2_19  0.10500000  0.05410000 -0.17300000
## F2_20 -0.02441415  0.06343181  0.06627665

#软阈值筛选
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), metSSd = "average")
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

write.table(dataExpr, file = "dataExpr.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

dataExpr <- read.table("dataExpr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


#接下来我们引入临床数据，可视化展示临床的量化指标，同时综合以上的聚类图。此处用到一个重要的函数numbers2colors。
#此函数主要对数值化的参数进行高低的颜色标记，形成相应的热图。我们直接看下实例：
traitData= read.table("phe.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T);
traitData <- traitData[rownames(dataExpr),]

traitData <- traitData[,c(1,2,3,4,5,8,11,12,20,21,22)]

write.table(traitData, file = "traitData.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


dim(traitData)
names(traitData)
#remove columns that SSld information we do not need.
allTraits= traitData[, -c(31, 16)];
allTraits= allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)

allTraits <- traitData[traitRows,]
# Forma data frame analogous to expression data that will SSld the clinical traits.
femaleSamples= rownames(datExpr);
traitRows= match(femaleSamples, allTraits$Mice);
traitRows <- femaleSamples
datTraits= allTraits[traitRows,c(1:6)];
rownames(datTraits)= allTraits[traitRows, 1]


datTraits$gender <- ifelse(datTraits$gender=='female',1,0)


datTraits <- traitData
datTraits <- as.numeric(datTraits)
#Re-cluster samples
sampleTree2= hclust(dist(dataExpr), metSSd = "average")
#Convert traits to a color representation: white means low, red means high, greymeans missing entry
traitColors= numbers2colors(datTraits, signed = FALSE);
# Plotthe sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2,traitColors,
                    groupLabels= names(datTraits),
                    main ="Sample dendrogram and trait heatmap")
par(oma=c(8,8,9,8),mar=c(12,12,12,12))
plot(hh,xaxt='n')
par(cex=0.6)
axis(1)



#计算软阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThresSSld(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)

## pickSoftThresSSld: will use block size 2697.
##  pickSoftThresSSld: calculating connectivity for given powers...
##    ..working on genes 1 through 2697 of 2697
##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
## 1      1   0.1370  0.825          0.412 587.000  5.95e+02  922.0
## 2      2   0.0416 -0.332          0.630 206.000  2.02e+02  443.0
## 3      3   0.2280 -0.747          0.920  91.500  8.43e+01  247.0
## 4      4   0.3910 -1.120          0.908  47.400  4.02e+01  154.0
## 5      5   0.7320 -1.230          0.958  27.400  2.14e+01  102.0
## 6      6   0.8810 -1.490          0.916  17.200  1.22e+01   83.7
## 7      7   0.8940 -1.640          0.869  11.600  7.29e+00   75.4
## 8      8   0.8620 -1.660          0.827   8.250  4.56e+00   69.2
## 9      9   0.8200 -1.600          0.810   6.160  2.97e+00   64.2
## 10    10   0.8390 -1.560          0.855   4.780  2.01e+00   60.1
## 11    12   0.8020 -1.410          0.866   3.160  9.61e-01   53.2
## 12    14   0.8470 -1.340          0.909   2.280  4.84e-01   47.7
## 13    16   0.8850 -1.250          0.932   1.750  2.64e-01   43.1
## 14    18   0.8830 -1.210          0.922   1.400  1.46e-01   39.1
## 15    20   0.9110 -1.180          0.926   1.150  8.35e-02   35.6
## 16    22   0.9160 -1.140          0.927   0.968  5.02e-02   32.6
## 17    24   0.9520 -1.120          0.961   0.828  2.89e-02   29.9
## 18    26   0.9520 -1.120          0.944   0.716  1.77e-02   27.5
## 19    28   0.9380 -1.120          0.922   0.626  1.08e-02   25.4
## 20    30   0.9620 -1.110          0.951   0.551  6.49e-03   23.5

par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft thresSSld (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft ThresSSld (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft thresSSld与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft ThresSSld (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")



power = sft$powerEstimate
power

## [1] 6

# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}


##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThresSSld = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)

net <- blockwiseModules(
  dataExpr,
  power = power,
  maxBlockSize = ncol(dataExpr),
  corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
  networkType = "unsigned",
  TOMType = "unsigned", 
  minModuleSize = 60,    ##越大模块越少
  mergeCutHeight = 0.01, ##越大模块越少
  numericLabels = TRUE, 
  saveTOMs= TRUE,
  saveTOMFileBase= "femaleMouseTOM",
  verbose = 3)

net <- recutBlockwiseTrees(dataExpr,
                           power = power,
                           maxBlockSize = ncol(dataExpr),
                           corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
                           networkType = "unsigned",
                           TOMType = "unsigned", 
                           minModuleSize = 30,    ##越大模块越少
                           mergeCutHeight = 0.25, ##越大模块越少
                           numericLabels = TRUE, 
                           saveTOMs= TRUE,
                           loadTOMs=TRUE,
                           saveTOMFileBase = paste0(exprMat, ".tom"),
                           verbose = 3)





##  Calculating module eigengenes block-wise from all genes
##    Flagging genes and samples with too many missing values...
##     ..step 1
##  ..Working on block 1 .
##     TOM calculation: adjacency..
##     ..will use 47 parallel threads.
##      Fraction of slow calculations: 0.000000
##     ..connectivity..
##     ..matrix multiplication (system BLAS)..
##     ..normalization..
##     ..done.
##    ..saving TOM for block 1 into file WGCNA/LiverFemaleClean.txt.tom-block.1.RData
##  ....clustering..
##  ....detecting modules..
##  ....calculating module eigengenes..
##  ....checking kME in modules..
##      ..removing 3 genes from module 1 because their KME is too low.
##      ..removing 5 genes from module 12 because their KME is too low.
##      ..removing 1 genes from module 14 because their KME is too low.
##  ..merging modules that are too close..
##      mergeCloseModules: Merging modules wSSse distance is less than 0.25
##        Calculating new MEs...

# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。 
table(net$colors)

## 
##   0   1   2   3   4   5   6   7   8   9  10  11  12  13 
## 135 472 356 333 307 303 177 158 102  94  69  66  63  62

#层级树展示各网络
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)


dataExpr <- read.table("dataExpr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

exp <- dataExpr
exp <- as.data.frame(t(exp))
exp$gene <- rownames(exp)
exp$group <- moduleColors

colors_group <-exp[,c(269,270)] 

write.table(colors_group, file = "colors_group.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(colors_group, file = "colors_group2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)




# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


save(net, moduleColors, file = "step3_genes_modules.Rdata")
save(net, moduleColors, file = "step3_genes_modules2.Rdata")

# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
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
png("step5_TOMplot_Network-heatmap.png",width = 800, height=600) 
TOMplot(plotTOM,geneTree,moduleColors,
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
        main="Network heapmap plot")

save(TOM, file = "TOM.Rdata")

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
        main = "Network heatmap plot, all genes")

#有时基因太多，需要限制基因的数量
nSelect = 1000

nGenes=15026
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), metSSd = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

TOMplot(plotDiss, selectTree, selectColors,
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
        main="Network heapmap plot")


## 如果有表型数据，也可以跟ME数据放一起，一起出图
traitData <- 
  traitData <- read.table("phe2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
traitData <- traitData[rownames(dataExpr),]

traitData$`diagnosis:ch1` <- ifelse(traitData$`diagnosis:ch1`=='NASH',1,0)


MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)



# 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
# 否则需要再计算一遍，比较耗费时间
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
        main = "Network heatmap plot, all genes")

save(TOM, file = "TOM.Rdata")

#导出用于cytoscape
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)

# Export the network into edge and node list files Cytoscape can read
# thresSSld 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
# cytoscape中再调整


cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste(exprMat, ".edges.txt", sep=""),
                               nodeFile = paste(exprMat, ".nodes.txt", sep=""),
                               weighted = TRUE, thresSSld = 0.5,
                               nodeNames = probes, nodeAttr = moduleColors)


#关联表型数据
trait <- read.table("phe3.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
trait <- read.table("GSE89632_phe.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
trait <- read.table("GSE48452_phe2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

trait <- read.table("GSE164760_phe.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

cosample <- intersect(rownames(trait),rownames(MEs_col))
trait <- trait[cosample,]
trait <- trait[rownames(MEs),]
MEs_col2 <- MEs_col[cosample,]
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



if (corType=="pearsoon") {
  modTraitCor = cor(MEs_col2, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col2, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
#win.graph(width=4.875, height=2.5,pointsize=8)

labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = T, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))

## 从上图可以看到MEmagenta与Insulin_ug_l相关
save(modTraitCor, moduleColors, MEs_col, modTraitP, file = "WGCNA.Rdata")
## 模块内基因与表型数据关联

# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达
# 值算出相关系数。
# 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要
# 。

### 计算模块与基因的相关性矩阵

if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}


# 计算性状与基因的相关性矩阵

## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

####LASSO回归重新设计####
#分割数据：

HCCgroup <- read.table("HCCgroup.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

data <- HCCgroup

set.seed(123)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))

HCCtrain <- data[ind==1, ] #the training data set

HCCtest <- data[ind==2, ] #the test data set

train <- rbind(HCtrain,HOtrain,SStrain)
train <- rbind(SStrain,NASHtrain,HCCtrain)

test <- rbind(HCtest,HOtest,SStest)
test <- rbind(SStest,NASHtest,HCCtest)

write.table(train, file = "train.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(train, file = "test.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



####lasso3####

setwd("E:/桌面文件转移/R/TCGAteech/肝脏疾病相关/NASH-HCC时间趋势模式/LASSO/HC-HO-SS")


install.packages('lars')
rm(list=ls())
options(stringsAsFactors = F)

Rdata_dir='Rdata/'
Figure_dir='figures/'
# 加载上一步从RTCGA.miRNASeq包里面提取miRNA表达矩阵和对应的样本临床信息。
load( file = 
        file.path(Rdata_dir,'TCGA-KIRC-miRNA-example.Rdata')
)
dim(expr)
dim(meta)
# 可以看到是 537个病人，但是有593个样本，每个样本有 552个miRNA信息。
# 当然，这个数据集可以下载原始测序数据进行重新比对，可以拿到更多的miRNA信息

# 这里需要解析TCGA数据库的ID规律，来判断样本归类问题。
group_list=ifelse(as.numeric(substr(colnames(expr),14,15)) < 10,'tumor','normal')

table(group_list)
exprSet=na.omit(expr)

table(group_list)
exprSet=na.omit(expr)
dim(exprSet)
load(  file = 
         file.path(Rdata_dir,'TCGA-KIRC-miRNA-survival_input.Rdata')
)
dim(exprSet) ## remove the nomral
head(phe)
exprSet[1:4,1:4]
head(colnames(exprSet))
head(phe$ID)
## 必须保证生存资料和表达矩阵，两者一致
all(substring(colnames(exprSet),1,12)==phe$ID)


library(lars) 
library(glmnet) 
x=t(log2(exprSet+1))#归一化
y=phe$event
#用基因的表达情况预测生死
model_lasso <- glmnet(x, y, family="binomial", nlambda=50, alpha=1)#拉手回归模型
print(model_lasso)
# 列%Dev代表了由模型解释的残差的比例，对于线性模型来说就是模型拟合的R^2(R-squred)。
# 它在0和1之间，越接近1说明模型的表现越好，
# 如果是0，说明模型的预测结果还不如直接把因变量的均值作为预测值来的有效。

head(coef(model_lasso, s=c(model_lasso$lambda[29],0.009)))
#找到最有意义的一个点纳入模型，就是说没有把所有基因都放入模型里面，只是找到了放入基因数使得模型最好的那个点
#这个29是用眼睛看来的，我们可以用下面第三个代码划线来看

plot(model_lasso, xvar = "norm", label = TRUE)

plot(model_lasso,xvar = "lambda", label = T)

cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
#这里选取了1000个，精确些
plot(cv_fit)
#再用得到的最佳的位置去建模，lASSO可以防止过度拟合


#lasso.prob是通过模型预测每个样本的stage值
model_lasso <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(y ,lasso.prob)


write.table(re, file = "lasso模型预测.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#得到预测结果
dat=as.data.frame(re[,1:2])
colnames(dat)=c('event','prob')
dat$event=as.factor(dat$event)#画图时需要factor
library(ggpubr) 
p <- ggboxplot(dat, x = "event", y = "prob",
               color = "event", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()
#得出预测结果

fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)#构建模型
fit$beta
#一倍SE内的更简洁的模型,是22个miRNA
#fit <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
#head(fit$beta)# 这里是40个miRNA
choose_gene=rownames(fit$beta)[as.numeric(fit$beta)!=0]
#选出用于建模的那些基因
length(choose_gene)
myexpr=x[,choose_gene]
mysurv=phe[,c("days","event")]
mysurv=phe2
mysurv$days[mysurv$days< 1] = 1 
# 详细代码参见这个网站https://github.com/jeffwong/glmnet/blob/master/R/coxnet.R#
fit <- glmnet( myexpr, Surv(mysurv$days,mysurv$event), 
               family = "cox") 
fit <- glmnet( myexpr, mysurv$diagnosis, 
               family = "cox")
#cox模型需要有时间和event数据
#用包自带的函数画图
plot(fit, xvar="lambda", label = TRUE)
plot(fit, label = TRUE)
## 如果需要打印基因名，需要修改函数，这里不展开。

library(pheatmap) 
choose_matrix=expr[choose_gene,]
choose_matrix[1:4,1:4]
n=t(scale(t(log2(choose_matrix+1))))  #scale()函数去中心化和标准化
#对每个探针的表达量进行去中心化和标准化
n[n>2]=2 #矩阵n中归一化后，大于2的项，赋值使之等于2（相当于设置了一个上限）
n[n< -2]= -2 #小于-2的项，赋值使之等于-2（相当于设置了一个下限）
n[1:4,1:4]

## http://www.bio-info-trainee.com/1980.html
y <- t(y)
y <- as.data.frame(y)
y$diagnosis<- as.factor(y$diagnosis)
annotation_col = as.data.frame(t(y))
rownames(annotation_col)=colnames(expr)

annotation_col <- ifelse(y$diagnosis=='0','HC','NASH')

y$diagnosis <- annotation_col
y <- order(y$diagnosis)


pheatmap(exp2,show_colnames = F,annotation_col = y,
         filename = 'lasso_genes.heatmap.png')



pheatmap(exp2,
         annotation_col = y,
         scale = "row",
         show_rownames = T,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=10,
         fontsize_col=3)




library(ggfortify)
df=as.data.frame(t(choose_matrix))
df$group=group_list
png('lasso_genes.pca.png',res=120)
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
dev.off()

## 也可以尝试其它主成分分析的R包，视频就不继续没完没了的讲解了。


library("FactoMineR")
library("factoextra")  
## 这里的PCA分析，被该R包包装成一个简单的函数，复杂的原理后面讲解。
dat.pca <- PCA(t(choose_matrix), graph = FALSE) #'-'表示“非”
fviz_pca_ind(dat.pca,repel =T,
             geom.ind = "point", # show points only (nbut not "text")只显示点不显示文本
             col.ind =  group_list, # color by groups 颜色组
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses 集中成椭圆
             legend.title = "Groups"
)

####AUC模型####
install.packages('ROCR')
library(ROCR)
setwd("E:/桌面文件转移/R/TCGAteech/肝脏疾病相关/NASH-HCC时间趋势模式/LASSO/HC-HO-SS/ROC")

re1 <- read.table("HC-HO.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
re2 <- read.table("HC-SS.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
re3 <- read.table("HO-SS.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

re1 <- read.table("SS-HASH.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
re2 <- read.table("SS-HCC.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
re3 <- read.table("NASH-HCC.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

re1 <- as.data.frame(re1)
re2 <- ifelse(re$group1=='2',1,0)

xy <- x
xy <- xy[rownames(re1),]
xy <- as.data.frame(xy)

re=cbind(re1 ,xy)
re$NASH <- re2
re <- re[re$diagnosis!='1',]
xy <- xy[xy$diagnosis!='1',]
re$PDZK1 <- xy$PDZK1
re$FADS2 <- xy$FADS2



pred_min1 <- prediction(re1[,2], re1[,1])
auc_min1 = performance(pred_min1,"auc")@y.values[[1]]

pred_min2 <- prediction(re2[,2], re2[,1])
auc_min2 = performance(pred_min2,"auc")@y.values[[1]]
pred_min3 <- prediction(re3[,2], re3[,1])
auc_min3 = performance(pred_min3,"auc")@y.values[[1]]

#求得AUC值
perf_min1 <- performance(pred_min1,"tpr","fpr")
perf_min2 <- performance(pred_min2,"tpr","fpr")
perf_min3 <- performance(pred_min3,"tpr","fpr")


par(mfrow=c(2,3))
plot(perf_min1,colorize=FALSE, col="red")
plot(perf_min2,colorize=FALSE, col="blue",add=T)
plot(perf_min3,colorize=FALSE, col="orange",add=T)

#绘图
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# y=x
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min1,3)))
text(0.8,0.3, labels = paste0("AUC = ",round(auc_min2,3)))
text(0.8,0.4, labels = paste0("AUC = ",round(auc_min3,3)))



legend("bottomright",
       c(paste0("AUC"),
         paste0("SS vs NASH : ",round(auc_min1,3)), 
         paste0("SS vs HCC  : ",round(auc_min2,3)), 
         paste0("NASH vs HCC: ",round(auc_min3,3))),
       col=c('white',"red", "blue", "orange"),
       lty=1, lwd=2,bty = "n")  

legend("bottomright",
       c(paste0("AUC"),
         paste0("HC vs HO : ",round(auc_min1,3)), 
         paste0("HC vs SS  : ",round(auc_min2,3)), 
         paste0("HO vs SS: ",round(auc_min3,3))),
       col=c('white',"red", "blue", "orange"),
       lty=1, lwd=2,bty = "n")  


plot(perf_min1,colorize=FALSE, col="red")+
  plot(perf_min2,colorize=FALSE, col="blue",add=T)+
  plot(perf_min3,colorize=FALSE, col="orange",add=T)+
  lines(c(0,1),c(0,1),col = "gray", lty = 4 )+
  legend("bottomright",
         c(paste0("AUC"),
           paste0("HC vs HO : ",round(auc_min1,3)), 
           paste0("HC vs SS  : ",round(auc_min2,3)), 
           paste0("HO vs SS: ",round(auc_min3,3))),
         col=c('white',"red", "blue", "orange"),
         lty=1, lwd=2,bty = "n")





# 加AUC值

pred_min1 <- prediction(re[,9], re[,3])
auc_min1 = performance(pred_min1,"auc")@y.values[[1]]

perf_min1 <- performance(pred_min1,"tpr","fpr")

plot(perf_min1,colorize=FALSE, col="blue")

lines(c(0,1),c(0,1),col = "gray", lty = 4 )

text(0.8,0.2, labels = paste0("AUC = ",round(auc_min1,3)))




rm(list = ls())
library(timeROC)
library(survival)

load(file = "../000files/timeROC.RData") # 这是我自己的数据，可以根据数据格式自己造一个



#下面是画图代码

ROC <- timeROC(T=df$futime,   
               delta=df$event,   
               marker=df$riskScore,   
               cause=1,                #阳性结局指标数值
               weighting="marginal",   #计算方法，默认为marginal
               times=c(1, 2, 3),       #时间点，选取1年，3年和5年的生存率
               iid=TRUE)

plot(ROC, 
     time=1, col="red", lwd=2, title = "")   #time是时间点，col是线条颜色
plot(ROC,
     time=2, col="blue", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC,
     time=3, col="orange", add=TRUE, lwd=2)

#添加标签信息
legend("bottomright",
       c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
         paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
         paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2,bty = "n")  

#柱状图美化

re <- read.table("re1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

dat <- dat[rownames(group),]
dat=as.data.frame(re[,1:2])
colnames(dat)=c('group','risk score')
dat$`group`=factor(dat$`group`,levels=c('SS',"NASH","HCC"))#画图时需要factor
dat$`group`=factor(dat$`group`,levels=c('HC',"HO","SS"))#画图时需要factor

dat$`risk score` <- as.numeric(dat$`risk score`)
library(ggpubr) 
p <- ggboxplot(dat, x = "group", y = "risk score",
               color = "group", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()
#得出预测结果



library(ggplot2)
risk <- read.table("lassoprobe.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
phe <- read.table("phe1.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
comname <- intersect(rownames(risk),rownames(phe))
risk <- risk[comname,]

phe <- phe[comname,]

data <- cbind(risk,phe)


filtered_data <- data[data$tumor_stage.diagnoses != "", ]

data$fibrosis_ishak_score <- factor(data$fibrosis_ishak_score,levels = c('0 - No Fibrosis','1,2 - Portal Fibrosis','3,4 - Fibrous Speta','5 - Nodular Formation and Incomplete Cirrhosis','6 - Established Cirrhosis'))
filtered_data$fibrosis_ishak_score <- factor(filtered_data$fibrosis_ishak_score,levels = c('0 - No Fibrosis','1,2 - Portal Fibrosis','3,4 - Fibrous Speta','5 - Nodular Formation and Incomplete Cirrhosis','6 - Established Cirrhosis'))


p = ggplot(data, aes(x=group, y='risk score',color =group, palette = "jco",
                     add = "jitter")) #这个时候再用scale_color就会无事发生


p 


