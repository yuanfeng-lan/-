
####多组t检验####

#读取16个m6a甲基化相关基因在CHOL中的表达量
m6a_expr_type=read.table(file="m6a_expr_with_type.txt",header=T,sep="\t",row.names=1)
exp <- read.table("exp.normalize.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- read.table("exp.normalize.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
phe <- read.table("phe.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

GSE66676exp <- read.table("GSE66676exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE66676group <- read.table("GSE66676group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- GSE66676exp
group <- GSE66676group

GSE126848exp <- read.table("GSE126848exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE126848group <- read.table("GSE126848group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- GSE126848exp
group <- GSE126848group

GSE48452exp <- read.table("GSE48452exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
GSE48452group <- read.table("GSE48452group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- GSE48452exp
group <- GSE48452group
group <- group[group$group==c('Nash','Control'),,drop=F]

gene <- c('NOL3','FAM169B','GRIA3','FGF14','ZNF878','FAM168A','FADS2')

exp1 <- exp[gene,]
exp1 <- exp1[-5,]
exp1 <- t(exp1)
exp1 <- exp1[rownames(group),]



exp1$group <- factor(c(rep("Tumor", 366), rep("Normal", 95)), levels = c("Tumor","Normal"))


group <- ifelse(group$diagnosis=='0','HC',
                ifelse(phe$diagnosis=='1','SS','NASH'))

group <- ifelse(group$group=='control','HC','NASH')


exp2 <- exp1
exp2 <- as.data.frame(exp2)

exp2$group <- group

exp2$group <-factor(exp2$group,ordered=TRUE,levels=c("NASH","HC")) #修改因子水平 

exp2 <- subset(exp2,group!='SS')
#获取16个m6a基因的名字，最后一列为样本类型
m6a_sym=names(m6a_expr_type)[1:(ncol(m6a_expr_type)-1)]

#生成一个空向量来存放计算出的p值
pval=c()

#for循环16次计算每个基因的p值
for(gene in exp2){
  #根据type来将样本分成两组
  p=t.test(exp2[,gene]~exp2$group)$p.value
  #存放p值
  pval=c(pval,p)
}
#输出p值看看
pval


#如果没有安装plyr和reshape2这两个R包，先去掉下面两行的#，运行进行安装
#BiocManager::install("plyr")
#BiocManager::install("reshape2")

#加载plyr和reshape2包
library(plyr)
library(reshape2)
#melt对m6a_expr_type数据格式进行转换
ddply(melt(m6a_expr_type),"variable",
      function(x) {
        w <- t.test(value~type,data=x)
        with(w,data.frame(statistic,p.value))
      })

#如果没有安装dplyr，rstatix和reshape2这三个R包，先去掉下面三行的#，运行进行安装
#BiocManager::install("dplyr")
#BiocManager::install("rstatix")
#BiocManager::install("reshape2")

#加载dplyr，rstatix和reshape2这三个R包
library(dplyr)
library(rstatix)
library(reshape2)
result=melt(m6a_expr_type) %>%
  group_by(variable) %>%
  t_test(value ~ type)

#输出result
result

#使用fdr方法对原始p值进行校正
result=melt(m6a_expr_type) %>%
  group_by(variable) %>%
  t_test(value ~ type) %>%
  adjust_pvalue(method = "fdr")

#一次性得到原始p值，FDR校正之后的p值以及转换成对应的***
result=melt(m6a_expr_type) %>%
  group_by(variable) %>%
  t_test(value ~ type) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")
#输出result
result


####多个性状的相关性分析及可视化####

install.packages("Hmisc")
install.packages("PerformanceAnalytics")
dd = as.data.frame(matrix(rnorm(1000),100,10))

head(dd)

# 计算相关系数及显著性
library(Hmisc)#加载包
res2 <- rcorr(as.matrix(exp1))
res2

# 可视化
library(PerformanceAnalytics)#加载包
chart.Correlation(exp1, histogram=T, pch=19)
chart
library(corrplot)
#绘制一个下三角的热图，这个包的使用在之前的博客写过，这里一笔带过
cor_matr = cor(exp1)
cor_matr
corrplot(cor_matr, type="upper", order="hclust", tl.col="black", tl.srt=45,addCoef.col = "grey")


####多组t检验作柱状图####
#载入R包
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
# 添加颜色
color <-c("#5CB85C","#D9534F")
color <-c("#FD8D62","#66C3A5")

#FD8D62
# 添加对比项
my_comparisons <- list(c("NASH", "HC"))
# 提取数据gene symbol
id <- colnames(exp2)
id <- id[-7]#去除最后一行的group
plist<-list()

#绘制拼图
# for循环
for (i in 1:length(id)){
  box<-exp2[,c(id[i],"group")]
  colnames(box)<-c("Expression","type")
  pb1<-ggboxplot(box,
                 x="type",
                 y="Expression",
                 color="type",
                 fill=NULL,
                 ylab = "Relative abundance",
                 add = "jitter",
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size=0.01,
                 font.label = list(size=30), 
                 palette = color)+
    theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
  pb1<-pb1+theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))
  pb1<-pb1+theme(axis.text.y = element_text(size = 12))+ggtitle(id[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1<-pb1+theme(legend.position = "NA")
  pb1<-pb1+stat_compare_means(method="t.test",hide.ns = F,comparisons =my_comparisons,label="p.signif")
  plist[[i]]<-pb1
} 
library(cowplot)
# 展示前6个gene的箱线图
plot<-plot_grid(plist[[1]],plist[[2]],plist[[3]],
                plist[[4]],plist[[5]],plist[[6]],
                ncol=3)
plot


#多组t检验，方差分析
exp <- read.table("normalize.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- read.table("choose_gene.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
gene <- gene$gene
group <- read.table("group.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
exp <- exp[,rownames(group)]
exp <- exp[gene,]
exp <- t(exp)
data <- cbind(exp,group)

group$group1 <- factor(group$group1,levels = c("SS","NASH","HCC"))




library(ggpubr) 
p <- ggboxplot(data, x = "group1", y = "CYP7A1",
               color = "group1", palette = "jco",
               add = "jitter")+
  scale_color_manual(values = c("#1a9781", "#48bad0", '#eb4b3a'))+ 
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

 p2 + stat_compare_means()


plot<-plot_grid(p,p2,
                ncol=,align = 'v')

plot
####统计描述####

##
## 批量T检验
##

## Author : Cdudu
## Data   : 2020 7/27


#清空目前环境中的变量
rm(list = ls()) 


library(readxl)
library(dplyr)
#读入数据
dat<-read_excel('sample.xlsx')
dat <- read.table("phe.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dat <- read.table("phe2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

dat$age<-as.numeric(dat$age)
#将分组转化为因子型
dat$diagnosis<-as.factor(dat$diagnosis)
summary(dat)

#创建新表用于存放Mean，SD和显著性标记
dat2<-data.frame(t1=as.character(1:5)) 

#输入dat中应变量（待检验指标）的起始列数
Star<-2

#输入dat中应变量（待检验指标）的终止列数
Over<-23
Over<-20
group_by(dat, diagnosis) %>%
  summarise(
    count = n(),
    mean = mean(dat[[2]], na.rm = TRUE),
    sd = sd(dat[[2]], na.rm = TRUE)
  )

count <- count(dat, diagnosis=='0')

#Mean,SD&P值计算
for ( i in c(Star:Over)){   
  
  means<-tapply(dat[[i]],dat$diagnosis,mean, na.rm = TRUE)#na.rm计算时去除NA值
  means<-sprintf('%.1f',round(means,1))  #小数点位数，注意'%.4f'和round()中的数字都要做对应修改
  SD<-tapply(dat[[i]],dat$diagnosis,sd,na.rm = TRUE)
  SD<-sprintf('%.1f',round(SD,1))       #小数点位数，注意'%.4f'和round()中的数字都要做对应修改
  M.t<-t.test(dat[[i]]~diagnosis,data=dat,var.equal=T) 
  pvalue<-M.t[[3]]
  if(pvalue>0.05){
    a<-paste(means,'±',SD)
    a[3]<-'NS'  
    a[4] <- M.t[["p.value"]]
  }           
  else if(pvalue>0.01){
    a<-paste(means,'±',SD)
    a[3]<-'*'
    a[4] <- M.t[["p.value"]]
  }
  else {
    a<-paste(means,'±',SD)
    a[3]<-'**'
    a[4] <- M.t[["p.value"]]
  }
  dat2[i-(Star-1)]<-a
  
}


names(dat2)[i-(Star-1)]<-names(dat[,i])

dat4 <- dat[,2:20]

colnames(dat2) <- colnames(dat4)
#行列转置
dat3<-t(dat2)  

#导出表格
write.csv(dat3,'数据统计Result2.csv')
M.t<-t.test(dat[[2]]~diagnosis,data=dat,var.equal=T) 



c <- with(dat, shapiro.test(dat[[3]][diagnosis == "0"]))

d <- with(dat, shapiro.test(dat[[3]][diagnosis == "2"]))

M.t<-wilcox.test(dat[[2]]~diagnosis,data=dat,var.equal=T,exact=FALSE)


dat2<-data.frame(t1=as.character(1:5)) 
for ( i in c(Star:Over)){   
  
  means<-tapply(dat[[i]],dat$diagnosis,mean, na.rm = TRUE)#na.rm计算时去除NA值
  means<-sprintf('%.1f',round(means,1))  #小数点位数，注意'%.4f'和round()中的数字都要做对应修改
  SD<-tapply(dat[[i]],dat$diagnosis,sd,na.rm = TRUE)
  SD<-sprintf('%.1f',round(SD,1))       #小数点位数，注意'%.4f'和round()中的数字都要做对应修改
  c <- with(dat, shapiro.test(dat[[i]][diagnosis == "0"]),na.rm = TRUE)
  d <- with(dat, shapiro.test(dat[[i]][diagnosis == "2"]))
  if(c$p.value>0.05&d$p.value>0.05){
    M.t<-t.test(dat[[i]]~diagnosis,data=dat,var.equal=T) 
    f <- 't.test'} else{
      M.t<-wilcox.test(dat[[i]]~diagnosis,data=dat,var.equal=T,exact=FALSE)
      f <- 'wilcox.test'
    }
  pvalue<-M.t[[3]]
  if(pvalue>0.05){
    a<-paste(means,'±',SD)
    a[3] <- f
    a[4]<-'NS'  
    a[5] <- M.t[["p.value"]]
  }           else if(pvalue>0.01){
    a<-paste(means,'±',SD)
    a[3] <- f
    a[4]<-'*'
    a[5] <- M.t[["p.value"]]
  }else {
    a<-paste(means,'±',SD)
    a[3] <- f
    a[4]<-'**'
    a[5] <- M.t[["p.value"]]
  }
  dat2[i-(Star-1)]<-a
  
}

#循环语句调整
i=20

M.t<-t.test(dat[[i]]~diagnosis,data=dat,var.equal=T) 
M.t<-wilcox.test(dat[[i]]~diagnosis,data=dat,var.equal=T,exact=FALSE)
c <- with(dat, shapiro.test(dat[[i]][diagnosis == "0"]),na.rm = TRUE)
d <- with(dat, shapiro.test(dat[[i]][diagnosis == "2"]))
if(c$p.value>0.05&d$p.value>0.05){
  M.t<-t.test(dat[[i]]~diagnosis,data=dat,var.equal=T) 
  f <- 't.test'} else{
  M.t<-wilcox.test(dat[[i]]~diagnosis,data=dat,var.equal=T,exact=FALSE)
    f <- 'wilcox.test'
  }
pvalue<-M.t[[3]]
if(pvalue>0.05){
  a<-paste(means,'±',SD)
  a[3] <- f
  a[4]<-'NS'  
  a[5] <- M.t[["p.value"]]
}           else if(pvalue>0.01){
  a<-paste(means,'±',SD)
  a[3] <- f
  a[4]<-'*'
  a[5] <- M.t[["p.value"]]
}else {
  a<-paste(means,'±',SD)
  a[3] <- f
  a[4]<-'**'
  a[5] <- M.t[["p.value"]]
}
dat2[i-(Star-1)]<-a



B = list()
for (i in c(Star:Over))
{
  d <- tryCatch(
    {with(dat, shapiro.test(dat[[i]][diagnosis == "2"])) },
    warning = function(w) { message('Waring @ ',i) ; return(NA) },
    error = function(e) { message('Error @ ',i) ; return(NA) },
    finally = { message('next...') }
  )}





for ( i in c(Star:Over))tryCatch({   
  
  means<-tapply(dat[[i]],dat$diagnosis,mean, na.rm = TRUE)#na.rm计算时去除NA值
  means<-sprintf('%.1f',round(means,1))  #小数点位数，注意'%.4f'和round()中的数字都要做对应修改
  SD<-tapply(dat[[i]],dat$diagnosis,sd,na.rm = TRUE)
  SD<-sprintf('%.1f',round(SD,1))       #小数点位数，注意'%.4f'和round()中的数字都要做对应修改
  
  c <- with(dat, shapiro.test(dat[[i]][diagnosis == "0"])) 
  d <- with(dat, shapiro.test(dat[[i]][diagnosis == "2"])) 
  if(c$p.value>0.05&d$p.value>0.05){
    M.t<-t.test(dat[[i]]~diagnosis,data=dat,var.equal=T) 
    f <- 't.test'}else{M.t<-wilcox.test(dat[[i]]~diagnosis,data=dat,var.equal=T,exact=FALSE) 
  f <- 'wilcox.test'}
  
  pvalue<-M.t[[3]]
  
  if(pvalue>0.05){
    a<-paste(means,'±',SD)
    a[3] <- f
    a[4]<-'NS'  
    a[5] <- M.t[["p.value"]]
  }           
  else if(pvalue>0.01){
    a<-paste(means,'±',SD)
    a[3] <- f
    a[4]<-'*'
    a[5] <- M.t[["p.value"]]
  }
  else {
    a<-paste(means,'±',SD)
    a[3] <- f
    a[4]<-'**'
    a[5] <- M.t[["p.value"]]
  }
  dat2[i-(Star-1)]<-a
  
},
error = function(e) {dat2[i-(Star-1)+1]<-b})

b <- c('0','0','0','0','0')
