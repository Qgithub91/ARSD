rm(list = ls())

pacman::p_load(
  tidyverse,
  spaa,
  EcolUtils,
  ggsignif,
  paletteer
)

meta <- readxl::read_xlsx("I:/rawdata/meta.xlsx", sheet = "meta") %>% 
  filter(regroup %in% c("H", "L")) %>% 
  select(sample, group) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))
otu.taxa <- read.table("I:/abundance_division/otu.taxa.tsv", header = T, row.names = 1, sep = "\t") %>% 
  select(otus, meta$sample, taxa, staxa) %>% 
  as.data.frame()
otu <- readxl::read_xlsx("I:/rawdata/otu_table.xlsx", sheet = "otutable")  
taxonomy <- readxl::read_xlsx("I:/rawdata/otu_table.xlsx", sheet = "taxonomy") 

# 预处理
treat.otu <- otu %>% 
  select(colnames(otu.taxa[,1:(ncol(otu.taxa)-2)])) %>%  # 减2是因为最后两列是taxa和staxa
  column_to_rownames("otus") %>% 
  # 将样本发现率低于10个样品的otu过滤掉
  mutate(counts = rowSums(.>0)) %>% 
  filter(counts >4) %>% 
  select(-counts) %>% 
  t() %>% 
  as.data.frame()

abun.otu <- otu %>% 
  select(colnames(otu.taxa[,1:(ncol(otu.taxa)-2)])) %>%  # 减2是因为最后两列是taxa和staxa
  column_to_rownames("otus") %>% 
  # 将样本发现率低于10个样品的otu过滤掉
  mutate(counts = rowSums(.>0)) %>% 
  filter(counts >4) %>% 
  select(-counts) %>% 
  mutate(across(everything(), ~ round(.x / sum(.x)*100, 10))) %>% 
  t() %>% 
  as.data.frame()

#niche.width.method = 'levins'，基于 Levins（1968）的公式计算生态位宽度指数；若要计算 Shannon 生态位宽度指数可修改此参数
#n = 1000，随机化重排 1000 次
#probs = c(0.025, 0.975)，计算双侧 95% 置信区间为准划分
set.seed(123)
spec_gen <- spec.gen(treat.otu, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
niche.width <- niche.width(treat.otu, method = "levins")
abun.niche.width <- niche.width(abun.otu, method = "levins")

# niche.width处理及绘图
niche.width2plot <- abun.niche.width %>% 
  pivot_longer(cols = everything(), names_to = "otus", values_to = "niche.width") %>% 
  left_join(otu.taxa, by = "otus") %>% 
  filter(staxa %in% c("AT", "RT")) %>%
  left_join(taxonomy, by = "otus") %>% 
  select(otus, staxa, niche.width, kingdom:species)

write.table(niche.width2plot, paste0("all.", format(Sys.time(), "%Y%m%d"), ".niche.width.tsv"), sep = "\t", row.names = F, quote = F)

# 绘制箱线图
# 颜色设置
fill_colors <- c("AT" = "#7A1B6D", "RT" = "#91D540")
point_colors <- c("AT" = "#4C1252", "RT" = "#5C8A2D")  # 更深的颜色用于散点

# 计算图形参数
y_max <- max(niche.width2plot$niche.width) * 1.15  # 扩展y轴范围

niche.width.boxplot <- ggplot(niche.width2plot, aes(x = staxa, y = niche.width, fill = staxa)) +
  geom_boxplot(alpha = 0.7, colour = 'black', width = 0.6) +
  geom_jitter(aes(fill = staxa), width = 0.2, shape = 21, size = 2.5, alpha = 0.01) +
  geom_signif(comparisons = list( c("AT", "RT")), test = "wilcox.test",
              map_signif_level = TRUE, textsize = 5, vjust = 0.1, 
              y_position = y_max * 0.95, tip_length = 0.01, na.rm = TRUE)+
  theme_minimal() +
  coord_cartesian(ylim = c(min(niche.width2plot$niche.width) * 0.9, y_max))+
  theme(legend.position = "none") +
  # geom_text(aes(label = sample, vjust = 0, color = group), size = 3.5, show.legend = F) + # 添加数据点的标签
  scale_fill_manual(values = fill_colors, name = "Group") +  # 设置名称为Area并保留图例
  scale_color_manual(values = point_colors, guide = "none") +  # 隐藏color的图例
  labs(x = "", y = "niche breadth") +
  theme_bw()+
  theme(panel.border = element_rect(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ 
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16,colour = "black"),
        legend.text = element_text(size = 14,colour = "black"), 
        legend.title = element_text(size = 16,colour = "black"))

niche.width.boxplot

ggsave(paste0("niche.width.", format("all.", Sys.time(), "%Y%m%d"), ".boxplot.pdf"), niche.width.boxplot, width = 5, height = 4, dpi = 300)

# spec.gen处理及绘图
all.spec.gen <- spec_gen %>% 
  rownames_to_column("otus") %>% 
  left_join(otu.taxa, by = "otus") %>% 
  filter(staxa %in% c("AT", "RT")) %>%
  left_join(taxonomy, by = "otus") %>% 
  select(otus, staxa, observed:sign, kingdom:species)
write.table(all.spec.gen, paste0("all.", format(Sys.time(), "%Y%m%d"), ".spec.gen.tsv"), sep = "\t", row.names = F, quote = F)

count.spec.gen <- all.spec.gen %>% 
  group_by(staxa, sign) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(
    sign = recode(
      sign,
      "NON SIGNIFICANT" = "Neutral",
      "SPECIALIST" = "Specialist",
      "GENERALIST" = "Generalist",
      .default = "Unknown"  # 处理未预期的值
    ) 
    )%>% 
  group_by(staxa) %>%
  mutate(
    total_count = sum(count),  # 计算每个 staxa 的总数量
    relative_abundance = count / total_count * 100,  # 计算百分比相对丰度
  ) %>%
  ungroup()

spec.gen.bar.plot <- ggplot(count.spec.gen, aes(x = staxa, y = relative_abundance, fill = sign)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() + # 经典主题
  labs(x = "", y = "Relative counts (%)",
       fill = "Sign") +
  scale_fill_manual(values = as.character(paletteer_d("awtools::bpalette")[c(5,3,7)]))+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
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

spec.gen.bar.plot
ggsave(paste0("all.", format(Sys.time(), "%Y%m%d"), ".spec.gen.bar.plot.pdf"), plot = spec.gen.bar.plot, height = 5, width = 5, dpi = 300)
