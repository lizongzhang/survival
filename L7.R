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

# 7.2.2 加载包 --------------------------------------------------------------------
library(randomForestSRC)  # 随机生存森林
library(TH.data)         # GBSG2 数据集
library(tidyverse)       # 数据处理和可视化
library(survival)        # 生存分析工具
library(pec)             # 生存概率预测

# 加载数据
data(GBSG2, package = "TH.data")

table(GBSG2$cens)

# 7.2.3 数据预处理 --------------------------------------------------------------------
gbsg2 <- GBSG2 %>% #对象名用小写，方便代码录入
  mutate(
    horTh = as.factor(horTh),  #在GBSG2中horTh是Factor,此处是示范如何将数值型转换为因子
    menostat = as.factor(menostat),
    tgrade = as.factor(tgrade)
  )

# 7.2.4 拟合RSF模型-------------------------------------------------------------------

## 设置模型参数-------------------------------------------------------------------

set.seed(123)   #让每次运行都能得到相同的结果
rsf_model <- rfsrc(Surv(time, cens) ~ horTh + age + menostat + 
                     tsize + tgrade + pnodes + progrec + estrec,
                   data = gbsg2,
                   ntree = 500,            # 树的数量
                   mtry = 3,               # 变量个数的平方根 
                   nodesize = 15,          # 终端节点最小样本数
                   importance = TRUE)      # 计算变量重要性

# 查看模型摘要
print(rsf_model)


## 计算 C-index-------------------------------------------------------------------

# C-index (判别能力)
cindex <- 1 - rsf_model$err.rate[length(rsf_model$err.rate)]
cindex


## 循环语句调参-------------------------------------------------------------------

# ntree    # 森林中树的数量，通常500-2000, 主要影响方差收敛，对C-index提升作用不大
# mtry     # 每次分裂考虑的变量数，可尝试sqrt、log2、全部变量
# nodesize # 终端节点最小样本数，调小可提升复杂度，防止欠拟合

# 网格搜索
library(randomForestSRC)
cindex_list <- c()
for (m in c(2, 3, 4, 5)) {
  for (n in c(5, 10, 15, 20)) {
    set.seed(123) 
    model <- rfsrc(Surv(time, cens) ~ horTh + age  + menostat + 
                     tsize + tgrade + pnodes + progrec + estrec,
                   data = gbsg2, ntree = 500, mtry = m, nodesize = n)
    cindex <- 1 - model$err.rate[length(model$err.rate)]
    print(paste("mtry:", m, "nodesize:", n, "C-index:", round(cindex, 3)))
    cindex_list <- c(cindex_list, cindex)
  }
}

# [1] "mtry: 2 nodesize: 15 C-index: 0.698"

set.seed(123)
rsf_model <- rfsrc(Surv(time, cens) ~ horTh + age + menostat + tsize + tgrade + pnodes + progrec + estrec,
                   data = gbsg2,
                   ntree = 500,            # 树的数量
                   mtry = 2,               # 变量个数的平方根 
                   nodesize = 15,          # 终端节点最小样本数
                   importance = TRUE
                   ) 
cindex <- 1 - rsf_model$err.rate[length(rsf_model$err.rate)]
cindex

# 7.2.5 RSF模型估计结果解读 ---------------------------------------------------

# 7.2.6 评估变量重要性 -----------------------------------------------------------

## 计算变量重要性 -----------------------------------------------------------

var_importance <- data.frame(
  Variable = names(rsf_model$importance),
  Importance = rsf_model$importance
) %>%
  arrange(desc(Importance))

var_importance



## 绘制变量重要性的条形图 ----------------------------------------------------

# 
ggplot(var_importance, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = "变量重要性（RSF）", x = "变量", y = "重要性") +
  theme_minimal()

# 7.2.7 风险分数 ------------------------------------------------------------

## 计算风险分数 ------------------------------------------------------------

gbsg2$risk_score  <- predict(rsf_model, newdata = gbsg2)$predicted

summary(gbsg2$risk_score)


## 风险分数直方图 ------------------------------------------------------------

gbsg2 %>% 
  ggplot(aes(x = risk_score)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
  labs(title = "RSF 风险分数的分布", x = "风险分数", y = "频数") +
  theme_minimal()

## 高中低风险组分层 --------------------------------------------------------------

gbsg2$risk_group <- cut(gbsg2$risk_score, 
                        breaks = quantile(gbsg2$risk_score, 
                                          probs = c(0, 0.33, 0.66, 1)), 
                        labels = c("Low", "Medium", "High"))

table(gbsg2$risk_group)

## 高中低风险组风险分数箱线图 --------------------------------------------------------------

gbsg2 %>% 
  filter(!is.na(risk_group)) %>%
  ggplot(aes(risk_group, risk_score, color = risk_group)) +
  geom_boxplot() + 
  labs(x = "Risk Group", y = "Risk Score")

# 7.2.8 高中低风险组分组KM曲线 -----------------------------------------------------

# Kaplan-Meier估计
surv_fit <- survfit(Surv(time, cens) ~ risk_group, data = gbsg2)

# 低中高风险组分组KM曲线
library(survminer)
ggsurvplot(surv_fit, 
           data = gbsg2,
           palette = "set1",  # 低、中、高风险颜色
           ggtheme = theme_minimal(),
           legend.title = "Risk Group",
           legend.labs = c("Low", "Medium", "High"),
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           conf.int = TRUE)  # 添加置信区间

#log-rank 检验
survdiff(Surv(time, cens) ~ risk_group, data = gbsg2)

# 7.2.9 预测患者的生存概率 --------------------------------------------------------

library(reshape2)
library(ggplot2)

# 假设 rsf_model 已经训练好，将gbsg2[1:5, ]的前五个样本视作 5个新样本
library(pec)
times <- seq(0, 2000, 100)
surv_pred <- predictSurvProb(rsf_model, 
                             newdata = gbsg2[1:5, ], 
                             times = times)
# 转成长数据框
df_long <- data.frame(
  time = rep(times, 5),
  surv = as.vector(t(surv_pred)),
  patient = rep(paste0("患者", 1:5), 
                each = length(times))
)

head(df_long)
tail(df_long)

# 7.2.10 绘制患者的生存曲线 ---------------------------------------------------------------

# 绘制5名患者的生存曲线
library(ggplot2)
ggplot(df_long, aes(x = time, y = surv, color = patient)) +
  geom_line(size = 1.2) +
  labs(title = "前5个患者的预测生存曲线", 
       x = "时间（天）", y = "生存概率", 
       color = "患者") +
  theme_minimal(base_size = 15) +
  scale_color_brewer(palette = "Set1")

# 7.4.1 GBM-Survival数据背景 ------------------------------------------------------------

data(GBSG2, package = "TH.data")

# 7.4.2 加载包-------------------------------------------------------------------------

# 加载必要的包
library(survival)        # 生存分析工具
library(TH.data)         # GBSG2 数据集
library(tidyverse)       # 数据处理和可视化
library(gbm)             # Generalized Boosted Regression Modeling

# 7.4.3 数据预处理-------------------------------------------------------------------------

# 数据预处理
gbsg2 <- GBSG2 %>%             # 保存到一个新的数据框，对象名称小写
  mutate(
    horTh = as.factor(horTh),  # R举例如何转成因子
    menostat = as.factor(menostat),
    tgrade = as.factor(tgrade)
  )

# 7.4.4 拟合GBM-Survival模型-------------------------------------------------------------------------

gbm_model <- gbm(Surv(time, cens) ~ horTh + age + menostat + tsize + 
                   tgrade + pnodes + progrec + estrec,
                 data = gbsg2,
                 distribution = "coxph",  # 使用 Cox 似然
                 n.trees = 100,           # 树的数量
                 interaction.depth = 3,   # 树的深度
                 shrinkage = 0.1,         # 学习率
                 bag.fraction = 0.5,      # 随机采样的比例
                 verbose = FALSE)

summary(gbm_model)

# 7.4.5 变量重要性评估-------------------------------------------------------------------------

# 计算变量相对影响
rel_inf <- summary(gbm_model)$rel.inf

rel_inf_df <- data.frame(variable = names(gbsg2)[1:8], 
                         importance = rel_inf)

rel_inf_df

# 绘制变量相对影响条形图
ggplot(rel_inf_df, 
       aes(x = reorder(variable, importance), 
           y = importance)) +
  geom_bar(stat = "identity", fill = "#A0D8EF") +
  coord_flip() +
  labs(title = "变量相对影响（GBM-Survival）", x = "变量", y = "相对影响") +
  theme_minimal()

# 7.4.6 计算风险分数-------------------------------------------------------------------------

gbsg2$risk_score_gbm <- predict(gbm_model, 
                                newdata = gbsg2, 
                                n.trees = 100)

summary(gbsg2$risk_score_gbm)

# 风险分数的直方图
gbsg2 %>% 
  ggplot(aes(risk_score_gbm)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "GBM-Survival 风险分数的分布", x = "风险分数", y = "频数") +
  theme_minimal()

# 7.4.7 单变量部份依赖图 ----------------------------------------------------------

library(pdp)
pd_1d <- partial(gbm_model, pred.var = "pnodes", 
                 n.trees = 100, train = gbsg2)

autoplot(pd_1d, rug = TRUE, train = gbsg2) +
  geom_line(color="navy") +
  geom_smooth(color = "lightblue",
              se = FALSE)

# 7.4.8 多变量部分依赖图 (2D PDP)-------------------------------------------------------------------------
pd_2d <- partial(gbm_model, 
                 pred.var = c("pnodes", "age"), 
                 n.trees = 100, chull = TRUE)

autoplot(pd_2d, contour = TRUE, legend.title = "Log Hazard")

# 7.4.9 个体条件期望曲线ICE Curves-------------------------------------------------------------------------

ice_curves <- partial(gbm_model, 
                      pred.var = "pnodes", 
                      ice = TRUE, 
                      n.trees = 100)

plotPartial(ice_curves, alpha = 0.1, rug = TRUE, train = gbsg2)

autoplot(ice_curves, center = FALSE, alpha = 0.1, pdp.color = "red")

# 7.4.10 中心化个体条件期望曲线 c-ICE Curves-------------------------------------------------------------------------
cice_curves <- partial(gbm_model, 
                       pred.var = "pnodes", 
                       ice = TRUE, 
                       center = TRUE, 
                       n.trees = 100)
autoplot(cice_curves, alpha = 0.1, pdp.color = "blue")
