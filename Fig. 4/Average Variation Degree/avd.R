rm(list = ls())

library(tidyverse)
library(ggsignif)

meta <- readxl::read_xlsx("J:/rawdata/meta.xlsx", sheet = "meta") %>% 
  filter(regroup %in% c("H", "L")) %>% 
  select(sample, group) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))
otu.taxa <- read.table("J:/abundance_division/otu.taxa.tsv", header = T, row.names = 1, sep = "\t") %>% 
  select(otus, meta$sample, taxa, staxa) %>% 
  as.data.frame()

# 数据处理及分类，去除低频率物种
target.abundance <- "AT" # AT或RT

target.otu <- otu.taxa %>% 
  filter(staxa == target.abundance) %>% 
  select(-taxa, -staxa) %>% 
  filter(rowSums(select_if(., is.numeric)) > 0) %>% 
  mutate(across(-otus, ~ round(.x / sum(.x)*100, 10))) %>% 
  column_to_rownames("otus") %>% 
  # 将样本发现率低于10个样品的otu过滤掉
  mutate(counts = rowSums(.>0)) %>% 
  filter(counts >4) %>% 
  select(-counts)

otu <- target.otu
group <- meta

#计算各 OTU 的变异度
ai <- abs(otu-apply(otu, 1, mean))/apply(otu, 1, sd)

#由于此时计算的是单个样本的 AVD，即 k=1
avd <- apply(ai,2,sum)/(1*nrow(otu))

#将AVD合并到分组信息中
group$avd <- avd

#异常值剔除
# 根据 target.abundance 的值进行过滤
if (target.abundance == "AT") {
  group <- group %>% filter(avd < 0.6)
} else {
  group <- group %>% filter(avd < 1)
}

# 绘制箱线图
# 颜色设置
fill_colors <- c("High" = "#7A1B6D", "Low" = "#91D540")
point_colors <- c("High" = "#4C1252", "Low" = "#5C8A2D")  # 更深的颜色用于散点

# 计算图形参数
y_max <- max(as.numeric(group$avd)) * 1.15  # 扩展y轴范围

p <- ggplot(group, aes(x = group, y = avd, fill = group))+
  geom_boxplot(alpha = 0.7, colour = 'black', width = 0.6)+
  geom_jitter(aes(fill = group), width = 0.2, shape = 21, size = 2.5, alpha = 0.7, colour = 'black') +
  geom_signif(comparisons = list( c("High", "Low")), test = "t.test",
              map_signif_level = TRUE, textsize = 5, vjust = 0.1,
              y_position = y_max * 0.95, tip_length = 0.01, na.rm = TRUE)+
  #facet_wrap(~area, ncol = 1) + #使用facet_wrap函数，其中~region表示根据region变量横向排列，ncol用来控制列数
  #scale_fill_brewer(palette = "Pastel2") + 
  labs(x= "", y = "AVD index") + #修改图例标题
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

ggsave(paste0(target.abundance, ".", format(Sys.time(), "%Y%m%d"), ".avd.boxplot.pdf"), p, width = 5, height = 5, dpi = 300)

