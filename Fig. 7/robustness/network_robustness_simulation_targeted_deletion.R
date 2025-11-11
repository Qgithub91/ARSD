rm(list = ls())

pacman::p_load(tidyverse)

meta <- readxl::read_xlsx("I:/工作/专利、文章撰写/高温大曲稀有和丰富物种分析/数据分析/原始数据/meta.xlsx", sheet = "meta") %>% 
  filter(regroup %in% c("H", "L")) %>% 
  select(sample, group) %>% 
  mutate(group = factor(group, levels = c("H", "L"), labels = c("High", "Low")))
otu.taxa <- read.table("I:/工作/专利、文章撰写/高温大曲稀有和丰富物种分析/数据分析/丰度划分/otu.taxa.tsv", header = T, row.names = 1, sep = "\t") %>% 
  select(otus, meta$sample, staxa) %>% 
  as.data.frame()
taxonomy <- readxl::read_xlsx("I:/工作/专利、文章撰写/高温大曲稀有和丰富物种分析/数据分析/原始数据/otu_table.xlsx", sheet = "taxonomy")  

# 数据处理及分类，去除低频率物种
# target.abundance <- "RT" # AT或RT
target.group <- "High"

## ------------------ 筛选目标 OTU ------------------
target.otu <- otu.taxa %>% 
  # filter(staxa == target.abundance) %>% 
  select(otus, meta$sample[meta$group==target.group]) %>% 
  column_to_rownames("otus") %>% 
  mutate(across(everything(), ~ round(.x / sum(.x), 10))) %>%
  filter(rowSums(across(where(is.numeric))) != 0) %>% 
  mutate(counts = rowSums(.>0)) %>% 
  filter(counts > 2) %>% 
  select(-counts)

comm <- t(target.otu)
sp.ra <- colMeans(comm)   # species relative abundance

## ------------------ 计算 Spearman 相关性矩阵 ------------------
cor_res <- rcorr(as.matrix(comm), type = "spearman")
r <- cor_res$r
p <- cor_res$P
diag(r) <- 0

## ------------------ RMT阈值计算 ------------------
calculate_rmt_threshold_WWTP <- function(cor_mat,
                                         thr_seq = seq(0.3, 1, 0.05),
                                         min_nodes = 20,
                                         verbose = TRUE) {
  results <- data.frame(thr=thr_seq, max_eig=NA_real_, n_nodes=NA_integer_)
  for(i in seq_along(thr_seq)){
    thr <- thr_seq[i]
    A <- abs(cor_mat)
    A[A<thr] <- 0
    diag(A) <- 0
    keep <- which(rowSums(A)>0)
    if(length(keep) < min_nodes) next
    subA <- A[keep, keep, drop=FALSE]
    results$max_eig[i] <- max(eigen(subA, symmetric=TRUE, only.values=TRUE)$values)
    results$n_nodes[i] <- length(keep)
  }
  med_me <- median(results$max_eig, na.rm=TRUE)
  sel_idx <- which(!is.na(results$max_eig) & results$max_eig < med_me)[1]
  thr_selected <- ifelse(is.na(sel_idx), NA, results$thr[sel_idx])
  list(threshold=thr_selected, results=results)
}

# rmt_res <- calculate_rmt_threshold_WWTP(r)
# r.constant <- rmt_res$threshold
# cat("RMT阈值 =", r.constant, "\n")

## ------------------ 阈值筛选 ------------------
r[abs(r) < 0.7] <- 0
r[abs(r) == 1] <- 0

p <- p.adjust(p, method = 'BH')
p[p >= 0.001] <- 0
p[p < 0.001] <- 1

z <- r * p
z1 <- matrix(0, nrow=nrow(z), ncol=ncol(z), dimnames=dimnames(z))
z1[] <- z

## ------------------ 构建 igraph 网络 ------------------
g <- graph_from_adjacency_matrix(as.matrix(z1), weighted = TRUE, mode='undirected')
g <- simplify(g)
g <- delete.vertices(g, names(degree(g)[degree(g)==0]))

# 删除小组件和高度密集可疑簇
comps <- components(g)
min_comp_size <- 5
small_nodes <- names(comps$membership[comps$membership %in% which(comps$csize < min_comp_size)])
if(length(small_nodes)>0) g <- delete.vertices(g, small_nodes)

comps <- components(g)
for(comp_id in unique(comps$membership)){
  nodes_in_comp <- names(comps$membership[comps$membership==comp_id])
  subg <- induced_subgraph(g, nodes_in_comp)
  if(vcount(subg)>=10 & edge_density(subg)>=0.8) g <- delete.vertices(g, nodes_in_comp)
}

g <- delete.vertices(g, names(degree(g)[degree(g)<2]))
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

network.raw <- as.matrix(as_adjacency_matrix(g, attr='weight'))
sp.ra2 <- sp.ra[rownames(network.raw)]

## ------------------ 读取关键种 ------------------
node.attri <- read.csv(paste0(target.group, ".0.7.20251106.zp.result.tsv"), header=TRUE, sep="\t")
module.hub <- as.character(node.attri$name[node.attri$z > 2.5 & node.attri$p <= 0.62])

module.hub.saving <- data.frame(otus=module.hub) %>% 
  left_join(taxonomy, by=c("otus"="otus")) %>% 
  left_join(otu.taxa[,c("otus","staxa")], by=c("otus"="otus"))

# write.table(module.hub.saving, paste0(target.group, ".", format(Sys.time(), "%Y%m%d"), ".module.hub.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

## ------------------ 关键种移除模拟 ------------------
rand.remov2.once <- function(netRaw, rm.num, keystonelist, sp.ra, abundance.weighted=TRUE){
  rm.num2 <- min(rm.num, length(keystonelist))
  id.rm <- sample(keystonelist, rm.num2)
  net.new <- netRaw[!rownames(netRaw) %in% id.rm, !colnames(netRaw) %in% id.rm]
  if(nrow(net.new)<2) return(0)
  sp.ra.new <- sp.ra[!names(sp.ra) %in% id.rm]
  net.strength <- if(abundance.weighted) net.new*sp.ra.new else net.new
  sp.meanInter <- colMeans(net.strength)
  while(length(sp.meanInter)>1 & min(sp.meanInter)<=0){
    id.remain <- which(sp.meanInter>0)
    net.new <- net.new[id.remain, id.remain]
    sp.ra.new <- sp.ra.new[id.remain]
    net.strength <- if(abundance.weighted) net.new*sp.ra.new else net.new
    sp.meanInter <- if(length(net.strength)>1) colMeans(net.strength) else 0
  }
  length(sp.ra.new)/length(sp.ra)
}

rmsimu <- function(netRaw, rm.p.list, keystonelist, sp.ra, abundance.weighted=TRUE, nperm=100){
  t(sapply(rm.p.list, function(x){
    remains <- replicate(nperm, rand.remov2.once(netRaw, x, keystonelist, sp.ra, abundance.weighted))
    c(remain.mean=mean(remains), remain.sd=sd(remains), remain.se=sd(remains)/sqrt(nperm))
  }))
}

Weighted.simu <- rmsimu(network.raw, 1:length(module.hub), module.hub, sp.ra2, TRUE, 100)
Unweighted.simu <- rmsimu(network.raw, 1:length(module.hub), module.hub, sp.ra2, FALSE, 100)

currentdat <- data.frame(
  Number.hub.removed = rep(1:length(module.hub),2),
  rbind(Weighted.simu, Unweighted.simu),
  weighted = rep(c("weighted","unweighted"), each=length(module.hub)),
  treat = rep(target.group, 2*length(module.hub))
)

write.table(currentdat, paste0(target.group, ".", format(Sys.time(), "%Y%m%d"), ".target_removal.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

## ------------------ 绘图 ------------------
ggplot(currentdat %>% filter(weighted=="weighted"), aes(Number.hub.removed, remain.mean, group=treat, color=treat)) +
  geom_line() + geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd), size=0.2) +
  scale_color_manual(values=c("red","blue")) +
  labs(x="Number of module hubs removed", y="Proportion of species remained") + theme_light()

ggplot(currentdat %>% filter(weighted=="unweighted"), aes(Number.hub.removed, remain.mean, group=treat, color=treat)) +
  geom_line() + geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd), size=0.2) +
  scale_color_manual(values=c("blue","red")) +
  labs(x="Number of module hubs removed", y="Proportion of species remained") + theme_light()
