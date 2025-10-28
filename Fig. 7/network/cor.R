rm(list = ls())

pacman::p_load(tidyverse, Hmisc, igraph,
               microeco)

meta <- readxl::read_xlsx("J:rawdata/meta.xlsx", sheet = "meta") %>% 
  filter(regroup %in% c("H", "L")) %>% 
  select(sample, group) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))
otu.taxa <- read.table("J:/abundance_division/otu.taxa.tsv", header = T, row.names = 1, sep = "\t") %>% 
  as.data.frame()
taxonomy <- readxl::read_xlsx("J:/rawdata/otu_table.xlsx", sheet = "taxonomy")  

# 数据处理及分类
target.group <- "Low" # High或Low

target.otu <- otu.taxa %>% 
  select(otus, meta[meta$group == target.group,]$sample) %>% 
  column_to_rownames("otus") %>% 
  filter(rowSums(across(where(is.numeric))) != 0) %>%
  # filter(rowSums(select_if(., is.numeric)) > 0.01) %>%
  # 将样本发现率低于10个样品的otu过滤掉
  mutate(counts = rowSums(.>0)) %>% 
  filter(counts >2) %>% 
  select(-counts)

#计算相关性
cor <- rcorr(as.matrix(t(target.otu)), type = "spearman")

#相关系数 r 值和显著性 p 值矩阵
r <- cor$r
diag(r) <- 0
p <- cor$P

#阈值筛选
#将 spearman 相关系数低于 r.constant 的关系剔除，即 r>=r.constant
#该模式下，一定要注意负值的选择是否是合适的，因为很多情况下可能负相关无意义
r.constant <- 0.98
r[abs(r) < r.constant] <- 0
r[abs(r) == 1] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.01] <- -1
p[p<0.01 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p

#再转换为对称矩阵，igraph 只能识别这种样式的邻接矩阵类型
z1 <- cor$r
z1[z1 != 0] <- 0
z1[rownames(z),colnames(z)] <- z
z1[colnames(z),rownames(z)] <- t(z)

#输出为excel文件
z2 <- z1
z2[is.na(z2)] =0
z2 <- z2 %>% 
  as.data.frame() %>% 
  select_if(function(x) sum(x)>0) %>% 
  filter(rowSums(.)>0) %>% 
  rownames_to_column()
writexl::write_xlsx(as.data.frame(z2), paste(target.group, format(Sys.time(), "%Y%m%d"), "cor", r.constant, "xlsx", sep = "."), col_names = T)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物丰度和功能基因丰度间的 spearman 相关系数
g <- graph.adjacency(z1, weighted = TRUE, mode = 'undirected')
#g

# 删除自相关
g <- simplify(g)

#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

# #查看网络图
# # plot(g)
# #计算群体结构（short random walks）；
# c<- cluster_walktrap(g)
# #使用默认颜色列表；
# V(g)$color <-c$membership+1
# #绘制网络图；
# layout<- layout.fruchterman.reingold
# #也可以尝试其他的Layout;
# layout<- layout_as_tree
# layout<- layout_nicely
# layout<- layout_on_sphere
# #……
# #绘制网络图；
# plot.igraph(g,vertex.size=4,
#             vertex.label=NA,
#             edge.curved=F,
#             edge.size=1.5,
#             layout= layout.fruchterman.reingold)

##网络文件输出，略，邻接矩阵、边列表等获取方式参考上文即可
#例如 gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
#write.graph(g, 'network.0.9.gml', format = 'gml')
#write.graph(g, 'network.0.9.graphml', format = 'graphml')

#边列表
edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation
)
head(edge_list)

writexl::write_xlsx(edge_list, paste(target.group, format(Sys.time(), "%Y%m%d"), "network.edge", r.constant, "xlsx", sep = "."), col_names = T)
# write.table(edge_list, 'network.edge_list(0.6).txt', sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(
  label = names(V(g))
) %>% 
  left_join(taxonomy, by = c("label" = "otus")) %>% 
  left_join(otu.taxa, by = c("label" = "otus")) %>% 
  select(label, kingdom, staxa)

head(node_list)

writexl::write_xlsx(node_list, paste(target.group, format(Sys.time(), "%Y%m%d"), "network.node", r.constant, "xlsx", sep = "."), col_names = T)
# write.table(node_list, 'network.node_list(0.6).txt', sep = '\t', row.names = FALSE, quote = FALSE)


