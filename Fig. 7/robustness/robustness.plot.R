rm(list = ls())

pacman::p_load(tidyverse, ggprism)

robustness.group <- "target_removal" #target_removal, random_removal

high.robustness <- read.table(paste0("High.20251106.", robustness.group, ".tsv"), sep = "\t", header = T) %>% 
  # filter(weighted == "weighted", Proportion.removed == 0.5)
  filter(weighted == "weighted", Number.hub.removed == 5)
low.robustness <- read.table(paste0("Low.20251106.", robustness.group, ".tsv"), sep = "\t", header = T) %>% 
  # filter(weighted == "weighted", Proportion.removed == 0.5)
  filter(weighted == "weighted", Number.hub.removed == 5)

robustness <- rbind(high.robustness, low.robustness)

robustness.plot <- ggplot(robustness, aes(x = treat, y = remain.mean, fill = treat))+
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7)+
  geom_errorbar(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd), 
                position = position_dodge(0.5), width = 0.15)+ #调整误差线长度
  # 确保Y轴从0开始并与X轴相交，expand = expansion(mult = c(0, 0))消除了默认的扩展空间，使得图表的边界紧贴数据的最小值和最大值。
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 1.1))+
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


# 计算Welch t-test的步骤

# 步骤1：从数据框提取统计量 ----------------------------------------
high_data <- high.robustness  # 提取High组数据行
low_data <- low.robustness    # 提取Low组数据行

# 步骤2：计算样本量n (通过标准误公式反推：SE = SD/sqrt(n) => n = (SD/SE)^2)
calc_n <- function(sd, se) (sd / se)^2

n_high <- calc_n(high_data$remain.sd, high_data$remain.se)  # High组样本量
n_low <- calc_n(low_data$remain.sd, low_data$remain.se)     # Low组样本量

# 步骤3：手动计算 Welch t-test --------------------------------------
# 定义计算函数
welch_t_test <- function(m1, m2, s1, s2, n1, n2) {
  se_diff <- sqrt( (s1^2/n1) + (s2^2/n2) )
  t_value <- (m1 - m2) / se_diff
  df <- ( (s1^2/n1 + s2^2/n2)^2 ) / ( (s1^4)/(n1^2*(n1-1)) + (s2^4)/(n2^2*(n2-1)) )
  p_value <- 2 * pt(abs(t_value), df, lower.tail = FALSE)
  return(list(t = t_value, df = df, p = p_value))
}

# 执行计算
result <- welch_t_test(
  m1 = high_data$remain.mean,
  m2 = low_data$remain.mean,
  s1 = high_data$remain.sd,
  s2 = low_data$remain.sd,
  n1 = n_high,
  n2 = n_low
)

# 步骤4：格式化输出结果 -------------------------------------------
cat(
  "\nWelch Two Sample t-test (Direct from Data Frame)\n",
  "-----------------------------------------------\n",
  "t-value\t\t:", round(result$t, 4), "\n",
  "Degrees of freedom\t:", round(result$df, 1), "\n",
  "p-value\t\t:", format.pval(result$p, eps = 0.0001), "\n"
)

# 步骤5：验证：通过t.test函数直接验证 -------------------------------
# 生成符合统计量的模拟数据
set.seed(123)
high_sim <- rnorm(n_high, high_data$remain.mean, high_data$remain.sd)
low_sim <- rnorm(n_low, low_data$remain.mean, low_data$remain.sd)

# 执行R内置的Welch检验
validation <- t.test(high_sim, low_sim, var.equal = FALSE)
print(validation)

