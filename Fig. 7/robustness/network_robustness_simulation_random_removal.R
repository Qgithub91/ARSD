rm(list = ls())

pacman::p_load(tidyverse, Hmisc, igraph, microeco, readxl, writexl)

### ----------------------------
### 1. 读入数据与基本设置
### ----------------------------
meta <- readxl::read_xlsx("I:/工作/专利、文章撰写/高温大曲稀有和丰富物种分析/数据分析/原始数据/meta.xlsx", sheet = "meta") %>% 
  filter(regroup %in% c("H", "L")) %>% 
  select(sample, group) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))

otu.taxa <- read.table("I:/工作/专利、文章撰写/高温大曲稀有和丰富物种分析/数据分析/丰度划分/otu.taxa.tsv",
                       header = T, row.names = 1, sep = "\t") %>% as.data.frame()

taxonomy <- readxl::read_xlsx("I:/工作/专利、文章撰写/高温大曲稀有和丰富物种分析/数据分析/原始数据/otu_table.xlsx",
                              sheet = "taxonomy")

target.group <- "Low"  # High或Low
# target.abundance <- "RT" # AT或RT

### ----------------------------
### 2. 数据过滤与相关性计算
### ----------------------------
target.otu <- otu.taxa %>% 
  # filter(staxa == target.abundance) %>% 
  select(otus, meta[meta$group == target.group, ]$sample) %>% 
  column_to_rownames("otus") %>% 
  mutate(across(everything(), ~ round(.x / sum(.x), 10))) %>%
  filter(rowSums(across(where(is.numeric))) != 0) %>% 
  mutate(counts = rowSums(. > 0)) %>% 
  filter(counts > 2) %>% 
  select(-counts)

cor <- rcorr(as.matrix(t(target.otu)), type = "spearman")
r <- cor$r
diag(r) <- 0
p <- cor$P

### ----------------------------
### 3. 计算 RMT 阈值
### ----------------------------
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
    if (n_k < min_nodes) next
    subA <- A[keep, keep, drop = FALSE]
    ev <- eigen(subA, symmetric = TRUE, only.values = TRUE)$values
    results$max_eig[i] <- max(ev)
    results$n_nodes[i] <- n_k
    if (verbose) message("thr=", thr, " => 节点数=", n_k, " 最大特征值=", round(max(ev),3))
  }
  med_me <- median(results$max_eig, na.rm = TRUE)
  sel_idx <- which(!is.na(results$max_eig) & results$max_eig < med_me)[1]
  thr_selected <- ifelse(is.na(sel_idx), NA, results$thr[sel_idx])
  list(threshold = thr_selected, results = results)
}

rmt_res <- calculate_rmt_threshold_WWTP(r)
r.constant <- rmt_res$threshold
cat("RMT 计算得到的阈值 =", r.constant, "\n")

### ----------------------------
### 4. 根据阈值和p值筛选边
### ----------------------------
r[abs(r) < r.constant] <- 0
r[abs(r) == 1] <- 0
p <- p.adjust(p, method = 'BH')
p[p >= 0.001] <- -1
p[p < 0.001 & p >= 0] <- 1
p[p == -1] <- 0

z <- r * p
z[is.na(z)] <- 0

### ----------------------------
### 5. 构建 igraph 网络并清理
### ----------------------------
g <- graph_from_adjacency_matrix(z, weighted = TRUE, mode = 'undirected')
g <- simplify(g)
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

# 删除小组件
comps <- components(g)
min_comp_size <- 5
small_comp_ids <- which(comps$csize < min_comp_size)
small_nodes <- names(comps$membership[comps$membership %in% small_comp_ids])
if(length(small_nodes) > 0){
  g <- delete.vertices(g, small_nodes)
}

# 删除高度密集可疑簇
comps <- components(g)
for(comp_id in unique(comps$membership)){
  nodes_in_comp <- names(comps$membership[comps$membership==comp_id])
  subg <- induced_subgraph(g, nodes_in_comp)
  if(vcount(subg) >= 10 & edge_density(subg) >= 0.8){
    g <- delete.vertices(g, nodes_in_comp)
  }
}
g <- delete.vertices(g, names(degree(g)[degree(g) < 2]))

### ----------------------------
### 6. 生成 network.raw 矩阵 (用于后续模拟)
### ----------------------------
network.raw <- as_adjacency_matrix(g, attr = "weight", sparse = FALSE)
diag(network.raw) <- 0

# 物种相对丰度（按网络节点匹配）
sp.ra <- colMeans(t(target.otu))
sp.ra2 <- sp.ra[names(V(g))]
sp.ra2 <- sp.ra2 / sum(sp.ra2, na.rm = TRUE)

### ----------------------------
### 7. 随机移除鲁棒性模拟部分 (原样保留)
### ----------------------------
rand.remov.once <- function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm <- sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw <- netRaw
  net.Raw[id.rm,] <- 0; net.Raw[,id.rm] <- 0
  net.stength <- if (abundance.weighted) net.Raw * sp.ra else net.Raw
  sp.meanInteration <- colMeans(net.stength)
  id.rm2 <- which(sp.meanInteration <= 0)
  remain.percent <- (nrow(netRaw) - length(id.rm2)) / nrow(netRaw)
  remain.percent
}

rmsimu <- function(netRaw, rm.p.list, sp.ra, abundance.weighted=T, nperm=100){
  t(sapply(rm.p.list, function(x){
    remains <- sapply(1:nperm, function(i){
      rand.remov.once(netRaw, rm.percent=x, sp.ra, abundance.weighted)
    })
    remain.mean <- mean(remains)
    remain.sd <- sd(remains)
    remain.se <- sd(remains)/(nperm^0.5)
    c(remain.mean=remain.mean, remain.sd=remain.sd, remain.se=remain.se)
  }))
}

Weighted.simu <- rmsimu(network.raw, seq(0.05,1,by=0.05), sp.ra2, TRUE, 100)
Unweighted.simu <- rmsimu(network.raw, seq(0.05,1,by=0.05), sp.ra2, FALSE, 100)

dat1 <- data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),
                   rbind(Weighted.simu,Unweighted.simu),
                   weighted=rep(c("weighted","unweighted"),each=20),
                   treat=rep(target.group,40))

write.table(dat1, paste0(target.group, ".", format(Sys.time(), "%Y%m%d"), ".random_removal.tsv"),
            sep = "\t", row.names = F, quote = F)

### ----------------------------
### 8. 可视化
### ----------------------------
ggplot(dat1[dat1$weighted=="weighted",], aes(x=Proportion.removed, y=remain.mean, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  theme_light()+
  labs(x="Proportion of species removed", y="Proportion of species remained")

ggplot(dat1[dat1$weighted=="unweighted",], aes(x=Proportion.removed, y=remain.mean, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  theme_light()+
  labs(x="Proportion of species removed", y="Proportion of species remained")
