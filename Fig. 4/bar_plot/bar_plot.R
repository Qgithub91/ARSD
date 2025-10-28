rm(list = ls())

library(tidyverse)
library(paletteer)

meta <- readxl::read_xlsx("J:/rawdata/meta.xlsx", sheet = "meta") %>% 
  filter(regroup %in% c("H", "L")) %>% 
  select(sample, group) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))
otu.taxa <- read.table("J:/abundance_division/otu.taxa.tsv", header = T, row.names = 1, sep = "\t") %>% 
  as.data.frame()
taxonomy <- readxl::read_xlsx("J:/rawdata/otu_table.xlsx", sheet = "taxonomy")  

# 根据丰度划分结果分组
# abundance.otu <- otu.taxa[otu.taxa$staxa == "AT",] 
# rare.otu <- otu.taxa[otu.taxa$staxa == "RT",]

# 选择目标群落进行处理，
target.abundance <- "AT" # AT或RT
target.category <- "genus" # species或genus
target.group <- "Low" # High或Low

species <- otu.taxa %>% 
  filter(staxa == target.abundance) %>% 
  mutate(across(where(is.numeric), ~ round(.x / sum(.x)*100, 10))) %>% 
  left_join(taxonomy, by = "otus") %>% 
  select(!!sym(target.category), where(is.numeric)) %>%  # !!sym() 主要用于在动态生成变量名时使用。
  group_by(!!sym(target.category)) %>%
  summarise_all(sum) %>% 
  ungroup()
  
###### 筛选前10个丰度最高的物种 ######
top.species <- species %>% 
  # 计算平均值，根据平均值排序，取前9个物种。
  mutate(ave = rowMeans(.[,-1])) %>% 
  arrange(desc(ave)) %>% 
  head(14) %>% 
  select(-ave)

others <- species %>%
  # anti_join()函数返回一个数据框中与另一个数据框中没有匹配值的所有行。
  anti_join(top.species, by = target.category) %>% 
  # 将category列中的所有值替换为“Others”。
  # 使用 := 来进行动态赋值
  mutate(!!sym(target.category) := "Others") %>% 
  # 将category列中的所有值根据Others加和
  group_by(!!sym(target.category)) %>%
  summarise_all(sum)

top.species <- rbind(top.species, others)

###### 绘制堆积柱状图 ######
iq.species <- top.species %>% 
  pivot_longer(-!!sym(target.category), names_to = "sample", values_to = "value") %>% 
  filter(sample %in% meta$sample[meta$group == target.group]) %>%
  mutate(!!sym(target.category) := factor(!!sym(target.category), unique(!!sym(target.category))))

iq.bar.plot <- ggplot(iq.species, aes(x = sample, y = value, fill = !!sym(target.category))) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() + # 经典主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Relative abundance (%)",
       fill = "Category") +
  scale_fill_manual(values = as.character(paletteer_d("awtools::bpalette")[1:15]))+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(color = "black"),  # 添加刻度线
    panel.background = element_blank(),
    axis.line = element_line(),  # 确保X轴和Y轴有线条
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 18),
    legend.text = element_text(color = "black", size = 14),
    legend.title = element_text(color = "black", size = 16)
  ) +
  # 确保Y轴从0开始并与X轴相交，expand = expansion(mult = c(0, 0))消除了默认的扩展空间，使得图表的边界紧贴数据的最小值和最大值。
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, NA))

iq.bar.plot
ggsave(paste0(target.abundance, ".", target.group, ".", target.category, ".bar.plot.pdf"),
       iq.bar.plot, width = 10, height = 6, dpi = 300)
