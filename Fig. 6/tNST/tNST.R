rm(list = ls())

pacman::p_load(tidyverse, NST, ggsignif)

meta <- readxl::read_xlsx("J:/rawdata/meta.xlsx", sheet = "meta") %>% 
  filter(regroup %in% c("H", "L")) %>% 
  select(sample, group) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))
otu.taxa <- read.table("J:/abundance_division/otu.taxa.tsv", header = T, row.names = 1, sep = "\t") %>% 
  as.data.frame()
taxonomy <- readxl::read_xlsx("J:/rawdata/otu_table.xlsx", sheet = "taxonomy")  
otu <- readxl::read_xlsx("J:/rawdata/otu_table.xlsx", sheet = "otutable")  

# 数据处理及分类
target.abundance <- "AT" # AT或RT

target.otu <- otu %>% 
  filter(otus %in% otu.taxa[otu.taxa$staxa == target.abundance, "otus"]) %>% 
  select(colnames(otu.taxa[,1:(ncol(otu.taxa)-2)])) %>%  # 减2是因为最后两列是taxa和staxa
  column_to_rownames("otus") %>% 
  select(meta$sample) %>% 
  # 将样本发现率低于10个样品的otu过滤掉
  mutate(counts = rowSums(.>0)) %>% 
  filter(counts >4) %>% 
  select(-counts) %>% 
  t() %>% 
  as.data.frame()

target.meta <- meta %>% 
  column_to_rownames("sample")

#计算 NST，以下仅展示了一些重要的参数项，更多详情 ?tNST
#指定 OTU 丰度表 target.otu 和样本分组 target.meta
#dist.method 用于指定相异度量，这里以 jaccard 相异指数为例
#null.model 用于选择零模型算法，这里使用默认值 PF
#rand = 1000，随机化次数，这里使用默认值 1000
#nworker 可用于指定多线程，以加快运算效率
set.seed(123)
tnst <- tNST(comm = target.otu, group = target.meta, dist.method = 'jaccard', null.model = 'PF', 
             rand = 1000, nworker = 12, abundance.weighted = T, output.rand = T, SES = T, 
             RC = T, between.group = F)

#输出结果 tnst 以列表类型存储，包含了计算的两两样本对之间的 NST
#names(tnst)
#tnst$ndex.pair  #所有样本对之间的 NST
#tnst$index.pair.grp  #仅组内样本对之间的 NST

#例如查看两组群落中，组内样本对的 NST
nst_group <- tnst$index.pair.grp
nst_group

# 查看随机性占比NST.i.bray
tnst[[2]]

#输出主要的统计结果
write.table(nst_group, paste0(target.abundance, ".", format(Sys.time(), "%Y%m%d"), ".nst.group.tsv"), sep = "\t", row.names = F, quote = F)
write.table(tnst[[2]], paste0(target.abundance, ".", format(Sys.time(), "%Y%m%d"), ".nst.proportion.tsv"), sep = "\t", row.names = F, quote = F)

# Bootstrapping test, 置换检验
tnstbt <- nst.boot(nst.result=tnst, group=target.meta, rand=999, trace=TRUE,
                   two.tail=FALSE, out.detail=TRUE, between.group=FALSE, nworker=1)

tnstbt[[2]]#NST
write.table(tnstbt[[2]], paste0(target.abundance, ".", format(Sys.time(), "%Y%m%d"), ".Bootstrapping.test.tsv"), sep = "\t", row.names = F, quote = F)

# Anova检验
tnstpaov <- nst.panova(nst.result=tnst, group = target.meta, rand = 999, trace = TRUE)

tnstpaov
write.table(tnstpaov, paste0(target.abundance, ".", format(Sys.time(), "%Y%m%d"), ".tnstpaov.tsv"), sep = "\t", row.names = F, quote = F)

# 绘制箱线图
# 颜色设置
fill_colors <- c("High" = "#7A1B6D", "Low" = "#91D540")
point_colors <- c("High" = "#4C1252", "Low" = "#5C8A2D")  # 更深的颜色用于散点

# 计算图形参数
y_max <- max(as.numeric(nst_group$NST.ij.ruzicka)) * 1.15  # 扩展y轴范围

p <- ggplot(nst_group, aes(x = group, y = NST.ij.ruzicka, fill = group))+
  geom_boxplot(alpha = 0.7, colour = 'black', width = 0.6)+
  geom_jitter(aes(fill = group), width = 0.2, shape = 21, size = 2.5, alpha = 0.7, colour = 'black') +
  geom_signif(comparisons = list( c("High", "Low")), test = "wilcox.test",
              map_signif_level = TRUE, textsize = 5, vjust = 0.1,
              y_position = y_max * 0.95, tip_length = 0.01, na.rm = TRUE)+
  #facet_wrap(~area, ncol = 1) + #使用facet_wrap函数，其中~region表示根据region变量横向排列，ncol用来控制列数
  #scale_fill_brewer(palette = "Pastel2") + 
  labs(x= "", y = "tNST") + #修改图例标题
  scale_fill_manual(values = fill_colors, name = "Area") +  # 设置名称为Area并保留图例
  scale_color_manual(values = point_colors, guide = "none") +  # 隐藏color的图例
  # geom_hline(yintercept = c(2, -2), color = "red", size = 1)+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ 
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16,colour = "black"),
        legend.text = element_text(size = 14,colour = "black"), 
        legend.title = element_text(size = 16,colour = "black"))

p

ggsave(paste0(target.abundance, ".", format(Sys.time(), "%Y%m%d"), ".nst.boxplot.pdf"), p, width = 5, height = 5, dpi = 300)

