# 准备：安装和加载包 ---------------------------------------------------------------

# 需要的包列表
pkgs <- c("survival", "ggsurvfit", "survminer", "rms", 
          "tidyverse", "broom", "glue", "gtsummary", 
          "flextable", "officer", "showtext", "RColorBrewer",
          "tidycmprsk", "MASS")

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




# 随机森林 --------------------------------------------------------------------

# 加载必要的包
library(randomForestSRC)  # 随机生存森林
library(TH.data)         # GBSG2 数据集
library(tidyverse)       # 数据处理和可视化
library(survival)        # 生存分析工具
library(pec)             # 生存概率预测

# 加载数据
data(GBSG2, package = "TH.data")

# 数据预处理：确保数据格式正确
gbsg2 <- GBSG2 %>%
  mutate(
    horTh = as.factor(horTh),      # 转换为因子
    menostat = as.factor(menostat),
    tgrade = as.factor(tgrade)
  )

# 拟合随机生存森林模型
rsf_model <- rfsrc(Surv(time, cens) ~ horTh + age + menostat + tsize + 
                     tgrade + pnodes + progrec + estrec,
                   data = GBSG2,
                   ntree = 500,              # 树的数量
                   mtry = 3,                 # 每次分裂随机选择的特征数
                   nodesize = 15,            # 终端节点最小样本数
                   importance = TRUE)        # 计算变量重要性

# 查看模型摘要
print(rsf_model)

# 提取变量重要性
var_importance <- data.frame(
  Variable = names(rsf_model$importance),
  Importance = rsf_model$importance
) %>%
  arrange(desc(Importance))

var_importance

# 可视化变量重要性
ggplot(var_importance, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +   # 修正stat参数
  coord_flip() +
  labs(title = "变量重要性（随机生存森林）", x = "变量", y = "重要性") +
  theme_minimal()

# 预测风险分数（ensemble mortality）
risk_scores <- predict(rsf_model, newdata = gbsg2)$predicted
risk_scores


# 可视化风险分数的分布
ggplot(data.frame(risk = risk_scores), aes(x = risk)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
  labs(title = "随机生存森林预测的风险分数分布", x = "风险分数", y = "频数") +
  theme_minimal()

# 风险分数分层
library(dplyr)
gbsg2$risk_group <- cut(risk_scores, 
                        breaks = quantile(risk_scores, 
                                          probs = c(0, 0.33, 0.66, 1)), 
                        labels = c("Low", "Medium", "High"))
gbsg2$risk_scores <- risk_scores
table(gbsg2$risk_group)

# 低中高风险组分组箱线图
library(ggplot2)
ggplot(subset(gbsg2, !is.na(risk_group)), 
       aes(x = risk_group, y = risk_scores)) +
  geom_boxplot(aes(color = risk_group)) + 
  labs(x = "Risk Group", y = "Risk Score")

gbsg2 %>% 
  filter(!is.na(risk_group)) %>%
  ggplot(aes(x = risk_group, y = risk_scores, color = risk_group)) +
  geom_boxplot() + 
  labs(x = "Risk Group", y = "Risk Score")

# 低中高风险组分组KM曲线

# 拟合 Kaplan-Meier 模型
surv_fit <- survfit(Surv(time, cens) ~ risk_group, data = gbsg2)

# 绘制 Kaplan-Meier 曲线
ggsurvplot(surv_fit, 
           data = gbsg2,
           palette = "set1",  # 低、中、高风险颜色
           ggtheme = theme_minimal(),
           legend.title = "Risk Group",
           legend.labs = c("Low", "Medium", "High"),
           xlab = "Time (Days)",
           ylab = "Survival Probability",
           conf.int = TRUE)  # 添加置信区间

#添加 p 值（log-rank 检验）
survdiff(Surv(time, cens) ~ risk_group, data = gbsg2)


# 预测患者的生存曲线 ---------------------------------------------------------------

library(reshape2)
library(ggplot2)


# 假设 rsf_model 已经训练好，gbsg2[1:5, ] 是5个新样本
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


# ggplot画多患者曲线
library(ggplot2)
ggplot(df_long, aes(x = time, y = surv, color = patient)) +
  geom_line(size = 1.2) +
  labs(title = "前5个患者的预测生存曲线", 
       x = "时间（天）", y = "生存概率", 
       color = "患者") +
  theme_minimal(base_size = 15) +
  scale_color_brewer(palette = "Set1")















# Melanoma ----------------------------------------------------------------

# 查看Melanoma数据集
data(Melanoma, package = "MASS")
head(Melanoma)

# 检查数据结构
str(Melanoma)

# 将status变量转换为生存分析需要的格式
# status: 1=死亡, 2=死亡, 3=生存
# 通常我们需要: event=1（死亡），censor=0（生存）
Melanoma$event <- as.numeric(Melanoma$status != 3)

# 拟合随机生存森林模型
rsf_model <- rfsrc(Surv(time, event) ~ ., data = Melanoma, ntree = 1000, 
                   importance = TRUE, na.action = "na.impute")
print(rsf_model)

Melanoma$risk_scores <- predict(rsf_model, newdata = Melanoma)$predicted
Melanoma$risk_group <- cut(risk_scores, 
                           breaks = quantile(risk_scores, probs = c(0, 0.33, 0.66, 1)),
                           labels = c("Low", "Medium", "High"),
                           include.lowest = TRUE)

surv_fit <- survfit(Surv(time, event) ~ risk_group, data = Melanoma)
ggsurvplot(surv_fit, data = Melanoma, 
           palette ="Set2",
           ggtheme = theme_minimal(), legend.title = "Risk Group",
           legend.labs = c("Low", "Medium", "High"),
           xlab = "Time (Days)", ylab = "Survival Probability",
           conf.int = TRUE)
survdiff(Surv(time, status) ~ risk_group, data = Melanoma)
