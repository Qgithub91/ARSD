rm(list = ls())

pacman::p_load(tidyverse, Hmisc, igraph,
               microeco)

meta <- readxl::read_xlsx("I:/工作/专利、文章撰写/高温大曲稀有和丰富物种分析/数据分析/原始数据/meta.xlsx", sheet = "meta") %>% 
  filter(regroup %in% c("H", "L")) %>% 
  select(sample, group) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))
otu.taxa <- read.table("I:/工作/专利、文章撰写/高温大曲稀有和丰富物种分析/数据分析/丰度划分/otu.taxa.tsv", header = T, row.names = 1, sep = "\t") %>% 
  as.data.frame()
taxonomy <- readxl::read_xlsx("I:/工作/专利、文章撰写/高温大曲稀有和丰富物种分析/数据分析/原始数据/otu_table.xlsx", sheet = "taxonomy")  

# 数据处理及分类
target.group <- "Low" # High或Low

target.otu <- otu.taxa %>% 
  select(otus, meta[meta$group == target.group,]$sample) %>% 
  column_to_rownames("otus") %>% 
  mutate(across(everything(), ~ round(.x / sum(.x), 10))) %>%
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

#-------------------------------
# RMT 阈值计算（参考 WWTP_community_network 方法）
#-------------------------------
calculate_rmt_threshold_WWTP <- function(cor_mat,
                                         thr_seq = seq(0.3, 1, by = 0.05),
                                         min_nodes = 20,
                                         verbose = TRUE) {
  results <- data.frame(thr = thr_seq, max_eig = NA_real_, n_nodes = NA_integer_)
  
  for (i in seq_along(thr_seq)) {
    thr <- thr_seq[i]
    A <- abs(cor_mat)
    A[A < thr] <- 0
    diag(A) <- 0
    keep <- which(rowSums(A) > 0)
    n_k <- length(keep)
    if (n_k < min_nodes) {
      if (verbose) message("thr=", thr, " 保留节点太少 (", n_k, "), 跳过")
      next
    }
    subA <- A[keep, keep, drop = FALSE]
    ev <- eigen(subA, symmetric = TRUE, only.values = TRUE)$values
    results$max_eig[i] <- max(ev)
    results$n_nodes[i] <- n_k
    if (verbose) message("thr=", thr, " => 节点数=", n_k, " 最大特征值=", round(max(ev),3))
  }
  
  # 简单选择第一个最大特征值小于中位数的阈值作为 RMT 阈值
  med_me <- median(results$max_eig, na.rm = TRUE)
  sel_idx <- which(!is.na(results$max_eig) & results$max_eig < med_me)[1]
  if (is.na(sel_idx)) {
    thr_selected <- NA
    warning("未找到明显的转折点阈值")
  } else {
    thr_selected <- results$thr[sel_idx]
  }
  
  if (verbose) {
    plot(results$thr, results$max_eig, type="b", pch=19,
         xlab="correlation threshold", ylab="max eigenvalue",
         main=sprintf("RMT threshold scan (selected = %s)", thr_selected))
    abline(v = thr_selected, col="blue", lty=2)
  }
  
  list(threshold = thr_selected, results = results)
}

# 计算 RMT 阈值
# rmt_res <- calculate_rmt_threshold_WWTP(r)
# r.constant <- rmt_res$threshold
# cat("WWTP 方法计算得到的相关性阈值 =", r.constant, "\n")

#阈值筛选
#将 spearman 相关系数低于 r.constant 的关系剔除，即 r>=r.constant
#该模式下，一定要注意负值的选择是否是合适的，因为很多情况下可能负相关无意义
r.constant <- 0.7
r[abs(r) < r.constant] <- 0
r[abs(r) == 1] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.001] <- -1
p[p<0.001 & p>=0] <- 1
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
g <- graph_from_adjacency_matrix(as.matrix(z1), weighted = TRUE, mode = 'undirected')
#g

# 删除自相关
g <- simplify(g)

#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

#-------------------------------
# 删除小独立簇和高度密集可疑簇
#-------------------------------
comps <- components(g)

# 1 删除小组件（节点数 < min_comp_size）
min_comp_size <- 5
small_comp_ids <- which(comps$csize < min_comp_size)
small_nodes <- names(comps$membership[comps$membership %in% small_comp_ids])
if(length(small_nodes) > 0){
  message(sprintf("删除 %d 个小组件节点（组件小于 %d）", length(small_nodes), min_comp_size))
  g <- delete.vertices(g, small_nodes)
}

# 2 删除高度密集可疑簇（dens >= 0.8 且节点数 >= 10）
comps <- components(g)
for(comp_id in unique(comps$membership)){
  nodes_in_comp <- names(comps$membership[comps$membership==comp_id])
  subg <- induced_subgraph(g, nodes_in_comp)
  n_nodes <- vcount(subg)
  dens <- edge_density(subg)
  if(n_nodes >= 10 & dens >= 0.8){
    message(sprintf("删除可疑高度密集子图：节点数=%d, 密度=%.3f", n_nodes, dens))
    g <- delete.vertices(g, nodes_in_comp)
  }
}


# 3 删除度低孤立节点（可选）
g <- delete.vertices(g, names(degree(g)[degree(g) < 2]))

#边列表
edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation
)
head(edge_list)

writexl::write_xlsx(edge_list, paste(target.group, target.abundance, format(Sys.time(), "%Y%m%d"), "network.edge", r.constant, "xlsx", sep = "."), col_names = T)
# write.table(edge_list, 'network.edge_list(0.6).txt', sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(
  ID = names(V(g)),
  label = names(V(g))
) %>% 
  left_join(taxonomy, by = c("label" = "otus")) %>% 
  left_join(otu.taxa, by = c("label" = "otus")) %>% 
  select(ID, label, kingdom, staxa)

head(node_list)

writexl::write_xlsx(node_list, paste(target.group, target.abundance, format(Sys.time(), "%Y%m%d"), "network.node", r.constant, "xlsx", sep = "."), col_names = T)
# write.table(node_list, 'network.node_list(0.6).txt', sep = '\t', row.names = FALSE, quote = FALSE)


