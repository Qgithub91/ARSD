rm(list = ls())

library(tidyverse)
library(readxl)
library(ggimage)
library(plyr)
library(Hmisc)
library(minpack.lm)
library(stats4)

meta <- readxl::read_xlsx("I:/rawdata/meta.xlsx", sheet = "meta") %>% 
  filter(regroup %in% c("H", "L")) %>% 
  select(sample, group) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))
otu.taxa <- read.table("I:/abundance_division/otu.taxa.tsv", header = T, row.names = 1, sep = "\t") %>% 
  as.data.frame()
taxonomy <- readxl::read_xlsx("I:/rawdata/otu_table.xlsx", sheet = "taxonomy")  
otu <- readxl::read_xlsx("I:/abundance_division/otu_table.xlsx", sheet = "otutable")  

# 数据处理及分类
target.abundance <- "RT" # AT或RT
target.group <- "High"

target.otu <- otu %>% 
  filter(otus %in% otu.taxa[otu.taxa$staxa == target.abundance, "otus"]) %>%
  select(colnames(otu.taxa[,1:(ncol(otu.taxa)-2)])) %>%  # 减2是因为最后两列是taxa和staxa
  select(otus, meta[meta$group == target.group,]$sample) %>% 
  column_to_rownames("otus") %>% 
  # 将样本发现率低于10个样品的otu过滤掉
  mutate(counts = rowSums(.>0)) %>% 
  filter(counts >2) %>% 
  select(-counts)

spp <- target.otu %>% 
  t() %>% 
  as.data.frame()

##将 Sloan 等（2006）的中性模型拟合到一个物种或分类群的丰度表，并返回几个拟合统计数据。或者，将根据它们在元群落中的丰度返回每个分类群的预测出现频率
#用非线性最小二乘法（Non-linear least squares，NLS）拟合模型参数
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
m.fit  #获取 m 值
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  #获取模型的 R2

#输出 3 个统计结果数据表，包括各物种或分类群的平均相对丰度（p.csv）、出现频率（freq.csv）和预测的出现频率（freq.pred.csv）
# write.table(p, paste0(target.group, ".", format(Sys.time(), "%Y%m%d"), ".p.tsv"), sep = "\t", row.names = F, quote = F)
# write.table(freq, paste0(target.group, ".", format(Sys.time(), "%Y%m%d"), ".freq.tsv"), sep = "\t", row.names = F, quote = F)
# write.table(freq.pred, paste0(target.group, ".", format(Sys.time(), "%Y%m%d"), ".freq.pred.tsv"), sep = "\t", row.names = F, quote = F)

#p 是平均相对丰度（mean relative abundance）
#freq 是出现频率（occurrence frequency）的观测值
#freq.pred 是出现频率（occurrence frequency）的预测值，即中性模型的拟合值

#绘制统计图，基础绘图法
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'#出现频率低于中性群落模型预测的部分
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'#出现频率高于中性群落模型预测的部分

#绘制统计图，ggplot2绘图法
#绘图数据整理
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
head(bacnlsALL)
#根据freq、Lower、Upper来分配颜色，并转为group组
bacnlsALL$group<-with(bacnlsALL,ifelse(freq<Lower,"#A52A2A",
                                       ifelse(freq>Upper,"#29A6A6","#f6d378")))
#开始绘图了
p1<-ggplot(data = bacnlsALL)+
  #添加freq.pred线和置信区间的线
  geom_line(aes(x=log(p),y=freq.pred),size=0.6,linetype =1)+
  geom_line(aes(x=log(p),y=Lower),size=0.6,linetype=2)+
  geom_line(aes(x=log(p),y=Upper),size =0.6,linetype = 2)+
  #下面代码是更改散点的大小和颜色
  geom_point(aes(x=log(p),y=freq,color =group),size =1)+
  xlab("log10(mean relative abundance)")+
  ylab("occurrence frequency")+
  scale_colour_manual(values =c("#A52A2A","#29A6A6","#f6d378"))+
  #下面两行代码可以更改标签的位置与大小
  annotate("text",x=-5,y=0,label=paste("Nm=",round(N,3)),size = 4)+
  annotate("text",x=-5,y=0.10,label = paste("R2=",round(Rsqr,3)),size=4)+
  #下面代码更改横纵坐标轴字体大小
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x=element_line(size =0.5,colour = "black"),
    axis.line.y =element_line(size =0.5,colour = "black"),
    axis.ticks = element_line(color = "black" ),
    axis.text =element_text(color ="black",size =12),
    legend.position = "none",
    legend.background= element_blank(),
    legend.key = element_blank(),
    legend.text=element_text(size =12),
    text =element_text(family="sans",size = 12),
    panel.border =element_rect(color ="black",fill = NA, size = 1)#添加图表边框
  )
p1
#整合一下数据框，添加饼图
data<- as.data.frame(table(bacnlsALL$group))
colnames(data)<- c("group","nums")
#重新添加一个type列
data$type<-(c("mid","low","high"))
#计算不同类型的的占比
data$percentage<- round(data$nums/sum(data$nums)*100, 1)
#将小数转为百分制
data$label<-paste(data$type,
                  paste(data$percentage,"%",sep=''))
#绘制饼图
p2<-ggplot(data,aes(x="",y= nums,fill = group))+
  geom_bar(stat ="identity",width=1)+
  coord_polar(theta ="y")+
  scale_fill_manual(
    values =c("#f6d378","#A52A2A","#29A6A6"),
    labels = data$label
  )+
  theme_void()+
  theme(
    panel.background =element_blank(),
    panel.grid =element_blank(),
    legend.background =element_blank(),
    legend.key = element_blank(),
    legend.text=element_text(size = 8)
  )
p2
#开始拼图
p<-p1+geom_subview(subview =p2,
                   #位置大小
                   x=-6,
                   y=0.3,
                   #下面两个是饼图大小
                   w = 2.5,
                   h=2.5)
p

ggsave(paste0(target.abundance, ".", target.group, ".", format(Sys.time(), "%Y%m%d"), ".assembly.pdf"), plot = p, height = 4, width = 5, dpi = 300)


