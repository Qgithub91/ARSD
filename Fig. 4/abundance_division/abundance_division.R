rm(list = ls())

library(readxl)
library(tidyverse)

meta <- read_xlsx("I:/rawdata/meta.xlsx", sheet = "meta") %>% 
  filter(group %in% c("H", "L")) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))
otu <- read_xlsx("I:/rawdata/otu_table.xlsx", sheet = "otutable") %>% 
  select(otus, meta$sample) %>% 
  filter(rowSums(select_if(., is.numeric)) > 0) %>% 
  mutate(across(-otus, ~ round(.x / sum(.x)*100, 10))) %>% 
  column_to_rownames("otus")

#（i）稀有类群（rare taxa，ART），在所有样本中丰度均 ≤0.1% 的 OTU
otu_ART <- otu[apply(otu, 1, function(x) max(x) <= 0.0001), ]

#（ii）丰富类群（abundant taxa，AT），在所有样本中丰度均 ≥1% 的 OTU
otu_AAT <- otu[apply(otu, 1, function(x) min(x) >= 0.01), ]

#（iii）中间类群（moderate taxa，MT），在所有样本中丰度均 >0.1% 且 <1% 的 OTU
otu_MT <- otu[apply(otu, 1, function(x) min(x) > 0.0001 & max(x) < 0.01), ]

#（iv）条件稀有类群（conditionally rare taxa，CRT），在所有样本中丰度均 <1%，且仅在部分样本中丰度 <0.1% 的 OTU
otu_CRT <- otu[apply(otu, 1, function(x) min(x) < 0.0001 & max(x) < 0.01), ]
otu_CRT <- otu_CRT[which(! rownames(otu_CRT) %in% rownames(otu_ART)), ]  #CRT 和 ART 是没有重叠的，在所有样本中丰度均 ≤0.1% 的 OTU 是不可取的

#（v）条件丰富类群（conditionally abundant taxa，CAT），在所有样本中丰度均 >0.1%，且仅在部分样本中丰度 >1% 的 OTU
otu_CAT <- otu[apply(otu, 1, function(x) min(x) > 0.0001 & max(x) > 0.01), ]
otu_CAT <- otu_CAT[which(! rownames(otu_CAT) %in% rownames(otu_AAT)), ]  #CAT 和 AAT 是没有重叠的，在所有样本中丰度均 ≥1% 的 OTU 是不可取的

#（vi）条件稀有或丰富类群（conditionally rare or abundant taxa，CRAT），丰度跨越从稀有（最低丰度 ≤0.1%）到丰富（最高丰度 ≥1%）的 OTU
otu_CRAT <- otu[apply(otu, 1, function(x) min(x) <= 0.0001 & max(x) >= 0.01), ]

#备注：这 6 个类群没有重叠，总数即等于 OTU 表的总数，相对丰度总和 100%
otu[which(rownames(otu) %in% rownames(otu_ART)),'taxa'] <- 'ART'
otu[which(rownames(otu) %in% rownames(otu_AAT)),'taxa'] <- 'AAT'
otu[which(rownames(otu) %in% rownames(otu_MT)),'taxa'] <- 'MT'
otu[which(rownames(otu) %in% rownames(otu_CRT)),'taxa'] <- 'CRT'
otu[which(rownames(otu) %in% rownames(otu_CAT)),'taxa'] <- 'CAT'
otu[which(rownames(otu) %in% rownames(otu_CRAT)),'taxa'] <- 'CRAT'

# 将6个类群继续细分为AT、RT和MT。
taxa.otu <- otu %>%
  mutate(staxa = case_when(
    taxa %in% c('CAT', 'AAT', 'CRAT') ~ 'AT',
    taxa == 'MT' ~ 'SMT', # 为避免与taxa中的MT重复，将MT改为SMT， 后期需要在AI中修改
    taxa %in% c('CRT', 'ART') ~ 'RT',
    TRUE ~ 'UNKNOWN' # 默认情况，当没有匹配的情况
  )) %>% 
  rownames_to_column("otus")

write.table(taxa.otu, 'otu.taxa.tsv', col.names = NA, sep = '\t', quote = FALSE)

# 一级分类统计、保存
otu_stat1 <- taxa.otu %>% 
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, .))) %>%  #将大于0的替换为1，表示存在，否则为0.
  select(-staxa) %>% 
  group_by(taxa) %>% 
  summarise(across(where(is.numeric), ~ sum(.x)))

write.table(otu_stat1, 'otu_stat1.tsv', col.names = NA, sep = '\t', quote = FALSE)

# 二级分类统计、保存
otu_stat2 <- taxa.otu %>% 
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, .))) %>%  #将大于0的替换为1，表示存在，否则为0.
  select(-taxa) %>% 
  group_by(staxa) %>% 
  summarise(across(where(is.numeric), ~ sum(.x)))

write.table(otu_stat2, 'otu_stat2.tsv', col.names = NA, sep = '\t', quote = FALSE)

# 一级分类绘图
otu_stat1 <- reshape2::melt(otu_stat1, id = 'taxa')

# svg("otu.taxa.svg", family = "MSYH", width = 18)
ggplot(otu_stat1, aes(variable, value, fill = taxa)) +
  geom_col(position = 'fill', width = 0.6) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'gray', fill = 'transparent')) +
  theme(axis.text.x = element_text(angle = 90, color = "black"))+
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = '样本', y = '不同丰富或稀有类群的占比')
# dev.off()

# 二级分类丰度绘图
abun.stat1 <- taxa.otu %>% 
  group_by(taxa) %>% 
  summarise(across(where(is.numeric), ~ sum(.x))) %>% 
  mutate(taxa = factor(taxa, levels = c("AAT", "CAT", "CRAT", "ART", "CRT", "MT"))) %>% 
  pivot_longer(cols = -taxa, names_to = "sample", values_to = "abundance")

fill_colors <- c("AAT" = "#420D3B", "CAT" = "#7A1B6D", "CRAT" = "#B32A9D",
                 "ART" = "#91D540", "CRT" = "#B8E67A",
                 "MT" = "#FDDB9F")

abun.plot1 <- ggplot(abun.stat1, aes(x = sample, y = abundance, fill = taxa))+
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() + # 经典主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Relative abundance (%)",
       fill = "Category") +
  scale_fill_manual(values = fill_colors, name = "Taxa2")+
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

abun.plot1

ggsave(paste0("taxa2.", ".barplot.pdf"), abun.plot1, width = 13, height = 6, dpi = 300)
