rm(list = ls())

pacman::p_load(tidyverse, ggprism)

robustness.group <- "target_removal"

high.robustness <- read.table(paste0("High.20250411.", robustness.group, ".tsv"), sep = "\t", header = T) %>% 
  filter(weighted == "weighted", Number.hub.removed == 5)
low.robustness <- read.table(paste0("Low.20250411.", robustness.group, ".tsv"), sep = "\t", header = T) %>% 
  filter(weighted == "weighted", Number.hub.removed == 5)

robustness <- rbind(high.robustness, low.robustness)

robustness.plot <- ggplot(robustness, aes(x = treat, y = remain.mean, fill = treat))+
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7)+
  geom_errorbar(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd), 
                position = position_dodge(0.7), width = 0.15)+ #调整误差线长度
  # 确保Y轴从0开始并与X轴相交，expand = expansion(mult = c(0, 0))消除了默认的扩展空间，使得图表的边界紧贴数据的最小值和最大值。
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 0.6))+
  labs(y = "Proportion of species remained", x = "Area")+
  scale_fill_manual(values = c("#7A1B6D", "#91D540"), guide = "none") +
  theme_bw()+#主题设置
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black", size = 14),
        axis.title = element_text(color="black", size = 16),
        legend.text = element_text(color="black", size = 14),
        legend.title = element_text(color="black", size = 16))    #设置分面字体)

robustness.plot

ggsave(paste0(robustness.group, ".", format(Sys.time(), "%Y%m%d"), ".robustness.plot.pdf"), robustness.plot, width = 4, height = 4, dpi = 300)
