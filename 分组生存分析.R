setwd("survival")
surv <- exp_sur
surv <- read.table("surv3.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T,comment.char = "")
surv <- read.csv("survgroup2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

surv$OS.time <- surv$OS.time/30
surv$OS <- as.numeric(unlist(surv$OS))#OS和OS.time数据numeric化，可数化。
surv$PDZK1 <- as.numeric(unlist(surv$PDZK1))#OS和OS.time数据numeric化，可数化。

#median或者平均数mean
median(surv$altered)
surv$group <- ifelse(surv$altered > median(surv$altered),"altered High","altered Low")
surv$group <- factor(surv$group, levels = c("altered Low","altered High")) 
class(surv$group)
table(surv$group)


median(surv$TNFRSF12A)
surv$group <- ifelse(surv$TNFRSF12A > median(surv$TNFRSF12A),"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)

median(surv$CYP7A1)
surv$group <- ifelse(surv$CYP7A1 > median(surv$CYP7A1),"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 



install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)

#3. 绘制生存曲线
#方法1
###3. 设置颜色，坐标
plot(fit, conf.int = T,
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time(Months)",
     ylab = "Survival probablity(%)"
)
###添加标签
legend("topright",
       title = "Group",
       c("Low risk+low fifbrosis", "other group"),
       lwd = 2, lty = 1,
       col = c("blue", "red"))
###添加P值
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
text(500, 0.1, p.lab)
#前两个参数表示粘贴的标签横纵坐标位置
dev.off()




km_fit <- survfit(Surv(OS.time,OS) ~ group, data=surv)

p <-ggsurvplot(km_fit,
                        xlab = "Time (Months)",
                        ylab="Proportion Alive",
                        pval = T,
                        conf.int = T,##置信带
                        risk.table = F,
                        legend.title = "",
                        legend.labs = levels(surv),##
                        #surv.median.line = "hv",# 中位生存
                        palette="lancet")


p


####方法2做生存分析图####
#install.packages("survminer")
library(survminer)


ggsurv <- ggsurvplot(
  km_fit,                     
  data = surv,             
  risk.table = F,       
  pval = TRUE,             
  conf.int = TRUE,         
  palette = c("#2b83ba","#d7191c"),
  xlim = c(0,500),         
  xlab = "Time in days",   
  break.time.by = 100,     
  ggtheme = theme_light(), 
  risk.table.y.text.col = T,
  risk.table.height = 0.25, 
  risk.table.y.text = FALSE,
  ncensor.plot = TRUE,      
  ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  
  legend.labs = c("Male", "Female")    
)
ggsurv



font.title = c(16, "bold", "darkblue"),
font.subtitle = c(15, "bold.italic", "purple"),
font.caption = c(14, "plain", "orange"),
font.x = c(14, "bold.italic", "red"),
font.y = c(14, "bold.italic", "darkred"),
font.tickslab = c(12, "plain", "darkgreen")




ggsurvplot(km_fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = F, # 显示风险表
           risk.table.col = "strata",
           palette = c("#2b83ba","#d7191c"), # 配色采用jco
           legend.labs = c("low", "high"), # 图例
           size = 0.5,
           xlim = c(0,120), # x轴长度，一般为0-10年
           break.time.by = 20, # x轴步长为20个月
           legend.title = "TNFRSF12A",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Months)", # 修改x轴标签
           #ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)

dev.off()


####分类群组####



####循环生存分析作图####
### for循环批量绘图


library(survival)
library(survminer)
surv <- read.table("TME_data.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
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

exp <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
survival <- read.table("surv3.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- read.table("choose_gene.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

group <- group[rownames()]




gene <- gene$gene

exp <- exp[gene,]
exp <- t(exp)
exp <- exp[rownames(survival),]

surv <- cbind(exp,survival)
write.table(data, file = "surv.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#方法2
plist<-list()
data <- surv

for(i in colnames(data)[1:22])
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
}
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

i <- colnames(data)[1]

group <- ifelse(data[[i]] > median(data[[i]]),"high","low")
diff <- survdiff(Surv(OS.time,OS) ~ risk, data = data)
pValue = 1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(OS.time,OS) ~ risk,data = data)
p <- ggsurvplot(fit,
                data = data,
                conf.int = TRUE,
                pval = pValue,
                pval.size = 5,
                legend.labs = c("High","Low"),
                legend.title = 'risk score',
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


p <- ggsurvplot(fit,
                data = surv,
                conf.int = TRUE,
                pval = pValue,
                pval.size = 5,
                legend.labs = c("High","Low"),
                
                xlab = "Time(months)",
                ylab = "Overall survival",
                break.time.by = 20,
                risk.table.title = "",
                palette = c("#d7191c","#2b83ba"),
                risk.table = F,
                risk.table.height = .25
)

ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           #risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = c("#d7191c","#2b83ba"), # 配色采用jco
           legend.labs = c("Low risk", "High risk"), # 图例
           size = 1,
           xlim = c(0,120), # x轴长度，一般为0-10年
           break.time.by = 20, # x轴步长为20个月
           legend.title = "",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Months)", # 修改x轴标签
           #ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)



####多组生存分析####
install.packages("RcolorBrewer")
install.packages("tibble")
install.packages("ggpp")
install.packages("survival")
install.packages("survminer")

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

surv <- read.table("surv2.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
fibro <- read.table("TCGA.fibro数据.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
surv <- read.table("survivalgroup.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)

surv <- surv[rownames(fibro),]
surv <- cbind(surv,fibro)
surv <- na.omit(surv)
surv$`fibrosis score` <- factor(surv$`fibrosis score`,levels = c('F0','F1,F2','F3,F4',"F5","F6"))


surv$OS.time <- surv$OS.time/30




dat$group <- sapply(dat$PAM50,function(x) {
  
  switch(x,
         
         "Basal" = "A", # 将Basal标记为A
         
         "Her2" = "B", # 将Her2标记为B
         
         "LumA" = "C", # 将LumA标记为C
         
         "LumB" = "D", # 将LumB标记为D
         
         "Normal" = "E")}) # 将Normal标记为E

##将生存时间转换为月份

dat$OS.time <- dat$OS.time * 12

# 生存分析
data <- surv
colnames(data) <- c("risk score" ,    "fibrosis_score", "OS"  , "OS.time" )

fitd <- survdiff(Surv(OS.time,OS) ~ fibrosis_score,
                 
                 data = data)

p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(Surv(OS.time, OS)~ fibrosis_score,
               
               data = data,
               
               type = "kaplan-meier",
               
               error = "greenwood",
               
               conf.type = "plain",
               
               na.action = na.exclude)

# 配对生存分析

ps <- pairwise_survdiff(Surv(OS.time, OS)~ fibrosis_score,
                        
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
                legend.labs = c("F0", "F1,F2","F3,F4","F5","F6"),
                
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
pdf("km_curve_with_pairwise_logrank.pdf", width = 4.5, height = 6)

print(p)

dev.off()