# 准备：安装和加载包 ---------------------------------------------------------------
  
# 需要的包列表
pkgs <- c("survival", "ggsurvfit", "survminer", "rms", 
            "tidyverse", "broom", "glue", "gtsummary", 
            "flextable", "officer", "showtext", "RColorBrewer",
            "tidycmprsk", "MASS",
            "randomForestSRC", "TH.data", "pec")
  
# 检查并安装未安装的包
for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
  
# 加载所有包
lapply(pkgs, library, character.only = TRUE)
  
# 加载showtext包，可以让R图形支持更多中文字体
library(showtext)
  
# 自动开启showtext功能，后续所有图形都支持中文
showtext_auto()
  
# 设置全局字体为'华文楷体'（STKaiti），适用于基础plot、par等绘图函数
par(family = 'STKaiti')
  

# 5.3.1 Weibull基线生存曲线-------------------------------------------------------------------------

library(TH.data)
data("GBSG2")
  
  
# 拟合Weibull分布生存模型(无协变量)
library(survival)
wb <- survreg(Surv(time, cens) ~ 1, data = GBSG2)

# 设置一系列生存概率
surv <- seq(0.99, 0.01, by = -0.01)

# 计算每个生存概率对应的生存时间
t <- predict(wb, 
             type = "quantile", 
             p = 1-surv, 
             newdata = data.frame(1))

# 整理结果为数据框
surv_wb <- data.frame(time = t, surv = surv)

# 查看前几行结果
head(surv_wb)

# 为可视化准备(添加上下界和标准误占位列，便于ggsurvplot_df绘图)
surv_wb <- data.frame(time = t, surv = surv, 
                      upper = NA, lower = NA, std.err = NA) 

library(survminer)
ggsurvplot_df(fit = surv_wb, surv.geom = geom_line)
  
# 5.3.2 估计Weibull 模型-------------------------------------------------------------------------

wb_fit <- survreg(Surv(time, cens) ~ horTh + tsize, GBSG2)

summary(wb_fit)  
 
# 5.3.3 从 Weibull 模型中计算生存时间指标-------------------------------------------------------------------------
  
# 计算肿瘤大小变量 tsize 的中位数（排除缺失值）
med_tsize <- median(GBSG2$tsize, na.rm = TRUE)

# 创建一个新的数据框
new_df <- data.frame(
  horTh = levels(GBSG2$horTh),
  tsize = med_tsize
)

new_df

# 预测中位数生存时间（50th percentile）
pred_median <- predict(wb_fit, 
                       newdata = new_df, 
                       type = "quantile",  #type = "quantile"表示预测的是生存时间的分位数
                       p = 0.5, 
                       se.fit = TRUE)


new_df_pred <- data.frame(
  horTh = new_df$horTh,
  tsize = new_df$tsize,
  median_time = pred_median$fit
)
new_df_pred

# 执行预测：0.25, 0.5, 0.75 分位数的生存时间
pred_mat <- predict(
  wb_fit,
  type = "quantile",
  p = c(0.25, 0.5, 0.75),
  newdata = new_df
)

# 整理为数据框，直接加上变量名
pred_df <- as.data.frame(pred_mat)
colnames(pred_df) <- c("q25", "q50", "q75")
pred_df$horTh <- new_df$horTh  # 添加分组标签

# 查看结果
pred_df
  
# 5.3.4 绘制不同组别的生存曲线-------------------------------------------------------------------------

# 1.拟合模型
wb_fit <- survreg(Surv(time, cens) ~ horTh + tsize, GBSG2)

# 2.虚拟患者数据
patients <- expand.grid(
  horTh = c("no", "yes"),
  tsize = as.numeric(quantile(GBSG2$tsize, c(0.25, 0.5, 0.75)))
)

# 3.生存概率
probs <- seq(0.99, 0.01, -0.01)

# 4.计算生存时间
times <- predict(wb_fit, 
                 type = "quantile", 
                 p = 1 - probs, 
                 newdata = patients)
  
# 5. 转成长数据
plotdata <- data.frame(
  horTh = rep(patients$horTh, each = length(probs)),
  tsize = rep(patients$tsize, each = length(probs)),
  time = as.vector(t(times)),
  surv = rep(probs, times = nrow(patients))
)

str(plotdata)

# 6. 画图
ggplot(plotdata, aes(x = time, y = surv, 
                     color = factor(tsize), 
                     linetype = horTh)) +
  geom_line(size = 1) +
  labs(x = "Time", y = "Survival probability",
       color = "Tumor size", linetype = "Hormone therapy",
       title = "Survival curves for 6 patients") +
  theme_minimal()

  
# 5.4 其他AFT模型-------------------------------------------------------------------------

# 指数模型
fit_exp <- survreg(Surv(time, cens) ~ horTh + tsize, data = GBSG2, dist = "exponential")
summary(fit_exp)

# 对数正态模型
fit_lnorm <- survreg(Surv(time, cens) ~ horTh + tsize, 
                     data = GBSG2, dist = "lognormal")
summary(fit_lnorm)


# 检查 $\log(T)$的分布

# 提取完整观测（未删失）样本
GBSG2 %>%
  filter(cens == 1) %>%                # 仅保留未删失数据
  mutate(logtime = log(time)) %>%     # 计算 log(生存时间)
  ggplot(aes(x = logtime)) +          # ggplot 不能用 %>% 链接
  geom_histogram(binwidth = 0.2, fill = "#69b3a2", color = "white") +
  labs(
    title = "log(生存时间) 的直方图",
    x = "log(生存时间)",
    y = "频数"
  ) +
  theme_minimal()


# 拟合三种分布，比较 AIC 或 BIC，选择值最低的模型。
fit_weibull <- survreg(Surv(time, cens) ~ horTh + tsize, data = GBSG2, dist = "weibull")

fit_exp <- survreg(Surv(time, cens) ~ horTh + tsize, data = GBSG2, dist = "exponential")

fit_lnorm <- survreg(Surv(time, cens) ~ horTh + tsize, data = GBSG2, dist = "lognormal")

AIC(fit_weibull, fit_exp, fit_lnorm)