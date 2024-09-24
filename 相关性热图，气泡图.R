####相关性及p值计算####
exp1 <- read.table("KIRCGSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp2 <- read.table("GSE53757GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp3 <- read.table("GSE73731GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group1 <- read.table("TCGAgroup.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- read.table("GSE53757group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group3 <- read.table("GSE73731group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
trems <- read.table("trems.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


exp1 <- read.table("LIHCGSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp2 <- read.table("GSE62232GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp3 <- read.table("GSE76427GSEA.exp.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

group1 <- read.table("LIHCgroup.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group2 <- read.table("GSE62232group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
group3 <- read.table("GSE76427group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
trems <- read.table("trems.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)



exp <- exp3
group <- group3

exp <- exp[,rownames(group)]
exp <- exp[trems$trems,]


Rfile <- exp[,1,drop=F]
Pfile <- exp[,1,drop=F]

i=1
c_r=cor(as.numeric(exp[i,]),as.numeric(group[,"PDZK1"]),method="pearson",use = 'pairwise.complete.obs')
p=cor.test(as.numeric(exp[i,]),as.numeric(group[,"PDZK1"]),method ="pearson")[[3]]
Rfile[i,1] <- c_r

for (i in 1:85){
  
    c_r=cor(as.numeric(exp[i,]),as.numeric(group[,"PDZK1"]),method="pearson",use = 'pairwise.complete.obs')
    p=cor.test(as.numeric(exp[i,]),as.numeric(group[,"PDZK1"]),method ="pearson")[[3]]
    ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
    Rfile[i,3] <- c_r
    Pfile[i,3] <- p
    
  
}

colnames(Rfile) <- c("TCGA.KIRC","GSE53757","GSE73731")
colnames(Pfile) <- c("TCGA.KIRC","GSE53757","GSE73731")

colnames(Rfile) <- c("TCGA.LIHC","GSE62232","GSE76427")
colnames(Pfile) <- c("TCGA.LIHC","GSE62232","GSE76427")




write.table(Pfile, file = "Pfile.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(Rfile, file = "Rfile.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


Rfile <- read.table("Rfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Pfile <- read.table("Pfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)




##ggcorrplot相关性图可视化
install.packages("hrbrthemes")



library(ggcorrplot)
library(ggtext)
library(RColorBrewer)
library(tidyverse)
library(hrbrthemes)

#将变量的首字母转为大写（主要为了看起来舒服一点，没什么用）
names(mtcars) = str_to_title(names(mtcars))

#计算相关性矩阵
corr_data <- round(cor(mtcars), 1)

#计算对应的p值
p_mat <- cor_pmat(mtcars)

#####相关性热力图#####

data <- Rfile

data <- as.data.frame(t(data))

data1 <- data[,c(1:25)]
data2 <- data[,c(26:50)]
data3 <- data[,c(61:85)]




pmt <- Pfile
pmt <- as.data.frame(t(pmt))


pmt1 <- pmt[,c(1:25)]
pmt2 <- pmt[,c(26:50)]
pmt3 <- pmt[,c(61:85)]
pmt1 <- as.matrix(pmt1)
pmt2 <- as.matrix(pmt2)
pmt3 <- as.matrix(pmt3)




ggcorrplot(data,
           outline.color = "black") +
  scale_fill_gradientn(colors =  brewer.pal(11, "Spectral"),
                       name = NULL)+
  hrbrthemes::theme_ipsum() +
  theme(  
    plot.title = element_markdown(color = "black", size=18),
    plot.subtitle = element_markdown(hjust = 0,vjust = .5, size=14),
    legend.key.height = unit(1, "null"),
    legend.key.width = unit(0.5, "cm"),
    legend.frame = element_rect(color="black", linewidth = 0.25),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = "white"))




ggcorrplot(data,
           outline.color = "black") +
  scale_fill_gradientn(colors =  brewer.pal(11, "Spectral"))

P1 <- ggcorrplot(data1,method = "circle",outline.color = "white",
           ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
           p.mat=pmt1,insig="pch",pch.col = "red", pch.cex = 5, tl.cex = 8,tl.col = 8,show.legend=F)

P2 <- ggcorrplot(data2,method = "circle",outline.color = "white",
                 ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
                 p.mat=pmt2,insig="pch",pch.col = "red", pch.cex = 5, tl.cex = 8,tl.col = 8,show.legend=F)

P3 <- ggcorrplot(data3,method = "circle",outline.color = "white",
                 ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
                 p.mat=pmt3,insig="pch",pch.col = "red", pch.cex = 5, tl.cex = 8,tl.col = 8)



#library(patchwork)

P1+P2+P3



cor <- corr.test(data1,data2,method = "spearman",adjust = "BH",ci = F)
cmt<-cor$r
pmt<-cor$p.adj

pmt1 <- as.matrix(pmt1)
ggcorrplot(data1,method = "circle",outline.color = "white",
           ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
           p.mat=pmt1,insig="pch",pch.col = "red", pch.cex = 3, tl.cex = 12)


#####气泡相关性热图#####

ggcorrplot(data,
           method = "circle",
           hc.order = TRUE,
           outline.color = "grey20") +
  scale_fill_gradientn(colors =  brewer.pal(11, "Spectral"),
                       name = NULL)+
  labs(x=NULL,y=NULL,
       title = "Example of <span style='color:#c1281a'>Correlation heat map</span>",
       subtitle = "draw charts with <span style='color:#03329a'>ggcorrplot()</span>") +
  hrbrthemes::theme_ipsum() +
  theme(  
    plot.title = element_markdown(color = "black", size=18),
    plot.subtitle = element_markdown(hjust = 0,vjust = .5, size=14),
    legend.key.height = unit(1, "null"),
    legend.key.width = unit(0.5, "cm"),
    legend.frame = element_rect(color="black", linewidth = 0.25),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = "white"),
    text = element_text(size = 15))

####右下角气泡相关性热图添加文字####
ggcorrplot(corr_data,
           method = "circle",
           type = "lower",
           lab = TRUE,
           lab_col = "grey20",
           lab_size = 3.5,
           outline.color = "grey20") +
  scale_fill_gradientn(colors =  brewer.pal(11, "Spectral"),
                       name = NULL)+
  labs(x=NULL,y=NULL,
       title = "Example of <span style='color:#c1281a'>Correlation heat map</span>",
       subtitle = "draw charts with <span style='color:#03329a'>ggcorrplot()</span>") +
  hrbrthemes::theme_ipsum() +
  theme(  
    plot.title = element_markdown(color = "black", size=18),
    plot.subtitle = element_markdown(hjust = 0,vjust = .5, size=14),
    legend.key.height = unit(1, "null"),
    legend.key.width = unit(0.5, "cm"),
    legend.frame = element_rect(color="black", linewidth = 0.25),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = "white"))

####给不显著的相关性打叉叉####

ggcorrplot(corr_data,
           type = "upper",
           p.mat = p_mat,
           outline.color = "grey20") +
  scale_fill_gradientn(colors =  brewer.pal(11, "Spectral"),
                       name = NULL)+
  labs(x=NULL,y=NULL,
       title = "Example of <span style='color:#c1281a'>Correlation heat map</span>",
       subtitle = "draw charts with <span style='color:#03329a'>ggcorrplot()</span>") +
  hrbrthemes::theme_ipsum() +
  theme(  
    plot.title = element_markdown(color = "black", size=18),
    plot.subtitle = element_markdown(hjust = 0,vjust = .5, size=14),
    legend.key.height = unit(1, "null"),
    legend.key.width = unit(0.5, "cm"),
    legend.frame = element_rect(color="black", linewidth = 0.25),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = "white"))






NES <- read.table("NES.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
P <- read.table("P.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
trems <- read.table("trems.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

data <- NES
data <- data[trems$trems,]
data <- as.data.frame(t(data))

data1 <- data[,c(1:25)]
data2 <- data[,c(26:50)]
data3 <- data[,c(61:85)]

pmt <- P
pmt <- pmt[trems$trems,]
pmt <- as.data.frame(t(pmt))


pmt1 <- pmt[,c(1:25)]
pmt2 <- pmt[,c(26:50)]
pmt3 <- pmt[,c(61:85)]
pmt1 <- as.matrix(pmt1)
pmt2 <- as.matrix(pmt2)
pmt3 <- as.matrix(pmt3)



P1 <- ggcorrplot(data1,method = "circle",outline.color = "white",
                 ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
                 p.mat=pmt1,insig="pch",pch.col = "red", pch.cex = 3, tl.cex = 8,tl.col = 8,show.legend=F)

P2 <- ggcorrplot(data2,method = "circle",outline.color = "white",
                 ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
                 p.mat=pmt2,insig="pch",pch.col = "red", pch.cex = 3, tl.cex = 8,tl.col = 8,show.legend=F)

P3 <- ggcorrplot(data3,method = "circle",outline.color = "white",
                 ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
                 p.mat=pmt3,insig="pch",pch.col = "red", pch.cex = 3, tl.cex = 8,tl.col = 8)



#library(patchwork)

P1+P2

library(corrplot)

data1 <- as.matrix(data1)
data2 <- as.matrix(data2)
data3 <- as.matrix(data3)
data1 <- data1[,-22]
pmt1 <- pmt1[,-22]
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
         tl.cex = 0.7,                  # 坐标轴标签字体的大小
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
         pch.cex = 3,                 # 图案大小
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
                     tl.cex = 0.7,                  # 坐标轴标签字体的大小
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
                     pch.cex = 3,                 # 图案大小
                     plotCI = "n",                # c("n", "square", "circle", "rect"),  # p值置信区间的方案
                     lowCI.mat = testRes$lowCI,   # p值置信区间下边界数据
                     uppCI.mat = testRes$uppCI,   # p值置信区间上边界数据
)

P3 <- corrplot(data3,
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
               tl.cex = 0.7,                  # 坐标轴标签字体的大小
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
               p.mat = pmt3,           # P值矩阵
               sig.level = 0.05,            # 当p大于sig.level时触发动作
               insig = "pch",               # p值大于sig.level时的方案"pch"图案, "p-value"P值数字, "blank"空白, "n"无操作, "label_sig"星号,
               pch = 4,                     # 当insig = "pch"时的图案形状 4为叉
               pch.col = "red",           # 图案颜色
               pch.cex = 3,                 # 图案大小
               plotCI = "n",                # c("n", "square", "circle", "rect"),  # p值置信区间的方案
               lowCI.mat = testRes$lowCI,   # p值置信区间下边界数据
               uppCI.mat = testRes$uppCI,   # p值置信区间上边界数据
)
library(cowplot)
plot<-plot_grid(P1,P2,P3,ncol=1,align = 'v')
plot


#脂肪酸代谢相关
Rfile1 <- read.table("TCGA.Rfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Rfile2 <- read.table("GSE53757.Rfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Rfile3 <- read.table("GSE73731.Rfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Rfile4 <- read.table("GSE126694.Rfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Rfile5 <- read.table("GSE215286.Rfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

Pfile1 <- read.table("TCGA.Pfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Pfile2 <- read.table("GSE53757.Pfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Pfile3 <- read.table("GSE73731.Pfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Pfile4 <- read.table("GSE126694.Pfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Pfile5 <- read.table("GSE215286.Pfile.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)


Rfile2 <- Rfile2[rownames(Rfile1),,drop=F]
Rfile3 <- Rfile3[rownames(Rfile1),,drop=F]
Rfile4 <- Rfile4[rownames(Rfile1),,drop=F]
Rfile5 <- Rfile5[rownames(Rfile1),,drop=F]

Pfile2 <- Pfile2[rownames(Pfile1),,drop=F]
Pfile3 <- Pfile3[rownames(Pfile1),,drop=F]
Pfile4 <- Pfile4[rownames(Pfile1),,drop=F]
Pfile5 <- Pfile5[rownames(Pfile1),,drop=F]

Rfile <- cbind(Rfile1,Rfile2,Rfile3,Rfile4,Rfile5)

Pfile <- cbind(Pfile1,Pfile2,Pfile3,Pfile4,Pfile5)

colnames(Rfile) <- c("TCGA.KIRC","GSE53757","GSE73731","GSE126694","GSE215286")
colnames(Pfile) <- c("TCGA.KIRC","GSE53757","GSE73731","GSE126694","GSE215286")


library(ggcorrplot)
library(ggtext)
library(RColorBrewer)
library(tidyverse)
library(hrbrthemes)

data <- Rfile
data <- as.data.frame(t(data))
data1 <- data[,c(1:38)]
data2 <- data[,c(39:76)]

pmt <- Pfile

pmt <- as.data.frame(t(pmt))


pmt1 <- pmt[,c(1:38)]
pmt2 <- pmt[,c(39:76)]
pmt1 <- as.matrix(pmt1)
pmt2 <- as.matrix(pmt2)




P1 <- ggcorrplot(data1,method = "circle",outline.color = "white",
                 ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
                 p.mat=pmt1,insig="pch",pch.col = "red", pch.cex = 3, tl.cex = 8,tl.col = 8,show.legend=F)

P2 <- ggcorrplot(data2,method = "circle",outline.color = "white",
                 ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
                 p.mat=pmt2,insig="pch",pch.col = "red", pch.cex = 3, tl.cex = 8,tl.col = 8)




library(patchwork)

P1+P2

