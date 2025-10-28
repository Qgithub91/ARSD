rm(list = ls())

library(tidyverse)
library(ggsankey)

data <- read.table("J:/abundance_division/otu_stat1.tsv", header = T, row.names = 1, sep = "\t") %>% 
  pivot_longer(cols = -1, names_to = "sample", values_to = "abundance") %>% 
  pivot_wider(names_from = "taxa", values_from = "abundance")
meta <- readxl::read_xlsx("J:/rawdata/meta.xlsx", sheet = "meta") %>% 
  filter(regroup %in% c("H", "L")) %>% 
  select(sample, group) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))
otu.taxa <- read.table("J:/abundance_division/otu.taxa.tsv", header = T, row.names = 1, sep = "\t") %>% 
  select(otus, meta$sample, taxa, staxa) %>% 
  as.data.frame()
taxonomy <- readxl::read_xlsx("J:/rawdata/otu_table.xlsx", sheet = "taxonomy")  

group.handle <- otu.taxa %>% 
  select(-taxa, -staxa) %>% 
  pivot_longer(cols = -1, names_to = "sample", values_to = "abundance") %>% 
  pivot_wider(names_from = "otus", values_from = "abundance") %>% 
  left_join(meta, by = "sample") %>% 
  group_by(group) %>% 
  summarise(across(where(is.numeric), mean, na.rm = T)) %>% 
  ungroup() %>% 
  # 使用mutate/across对除了第一列之外的所有列进行操作,将大于0的值替换为group
  mutate(across(where(is.numeric), ~ if_else(. > 0, as.character(group), as.character(.x)))) %>% 
  pivot_longer(cols = -1, names_to = "otus", values_to = "groups") %>% 
  filter(groups != 0)

bind.data <- otu.taxa %>% 
  select(otus, taxa, staxa) %>% 
  left_join(taxonomy, by = "otus") %>% 
  right_join(group.handle, by = "otus") %>% 
  select(taxa, staxa, kingdom, group)

sangkey.data <- bind.data %>% 
  make_long(kingdom, staxa, taxa, group) %>% 
  mutate(node = factor(node, levels = c(rev(c("Bacteria", "Eukaryota", "Archaea")),
                                        rev(c("AT", "SMT", "RT")),
                                        rev(c("AAT", "CAT", "CRAT", "MT", "CRT", "ART")),
                                        rev(c("High", "Low")))))

ggplot(sangkey.data, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = node,
               label = node)) +
  geom_sankey(flow.alpha = 0.5,  #条带不透明度
              node.color = 1,
              smooth = 7, #条带平滑度
              width = 0.2) + #节点宽度
  geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_viridis_d(option = "A", alpha = 0.95) +
  theme_sankey(base_size = 16) +
  theme(legend.position = 'none') #隐藏图例

ggsave("sangkey.pdf",width = 6, height = 4, dpi = 300)

# 统计bind.data各个节点的数量
count.data <- otu.taxa %>% 
  left_join(taxonomy, by = "otus") %>% 
  select(taxa, staxa, kingdom) %>% 
  mutate(kingdom2staxa = paste(kingdom, staxa, sep = "-")) %>% 
  mutate(staxa2taxa = paste(staxa, taxa, sep = "-"))

node.count <- lapply(count.data[,c("taxa", "staxa", "kingdom")], table)

node.count

# 统计两点对接的数量
# 使用paste函数连接两列，并创建新列node2nextnode
nodes.count <- lapply(count.data[,c("kingdom2staxa", "staxa2taxa")], table)

nodes.count
