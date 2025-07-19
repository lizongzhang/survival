
install.packages("survivalmodels")
install.packages("gbm")
# 加载必要的包
library(survival)        # 生存分析工具
library(TH.data)         # GBSG2 数据集
library(tidyverse)       # 数据处理和可视化
#library(survivalmodels)  # GBM-Survival 模型
library(gbm)             # Generalized Boosted Regression Modeling
library(pdp)

# 加载数据
data(GBSG2, package = "TH.data")

# 数据预处理
gbsg2 <- GBSG2 %>%
  mutate(
    horTh = as.factor(horTh),
    menostat = as.factor(menostat),
    tgrade = as.factor(tgrade)
  )

# 拟合 GBM-Survival 模型
gbm_model <- gbm(Surv(time, cens) ~ horTh + age + menostat + tsize + 
                   tgrade + pnodes + progrec + estrec,
                 data = gbsg2,
                 distribution = "coxph",  # 使用 Cox 似然
                 n.trees = 100,           # 树的数量
                 interaction.depth = 3,   # 树的深度
                 shrinkage = 0.1,         # 学习率
                 bag.fraction = 0.5,      # 随机采样的比例
                 verbose = FALSE)

# 查看模型摘要
summary(gbm_model)

# 提取变量相对影响
rel_inf <- summary(gbm_model)$rel.inf

rel_inf_df <- data.frame(variable = names(gbsg2)[1:8], 
                         importance = rel_inf)

ggplot(rel_inf_df, 
       aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "salmon") +
  coord_flip() +
  labs(title = "变量相对影响（GBM-Survival）", x = "变量", y = "相对影响") +
  theme_minimal()

# 预测风险分数
risk_scores <- predict(gbm_model, newdata = gbsg2, n.trees = 100)

# 可视化风险分数的分布
ggplot(data.frame(risk = risk_scores), aes(x = risk)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "GBM-Survival 风险分数的分布", x = "风险分数", y = "频数") +
  theme_minimal()



# 部分依赖图（示例：pnodes）
partial_dep <- plot(gbm_model, i.var = "pnodes", n.trees = 100)

str(partial_dep)

# 提取 x 和 y 数据
x_values <- partial_dep$panel.args[[1]]$x  # 存储了 pnodes 的网格值,对 pnodes 变量进行均匀采样或插值生成的一系列点
y_values <- partial_dep$panel.args[[1]]$y  # 存储了在每个 x_values（pnodes 网格点）上的部分依赖值，这些值表示 pnodes 对响应变量（通常是生存模型的预测输出）的平均影响，对数累积危险函数（log cumulative hazard），y_values 的变化反映 pnodes 对风险的贡献

# 创建数据框
y_values <- partial_dep$panel.args[[1]]$y
H_values <- exp(y_values)
pd_data <- data.frame(pnodes = partial_dep$panel.args[[1]]$x, H = H_values)
ggplot(pd_data, aes(x = pnodes, y = H)) +
  geom_line(color = "purple") +
  labs(title = "pnodes 的累积危险", x = "淋巴结数量", y = "累积危险") +
  theme_minimal()

S_values <- exp(-H_values)
pd_data$S <- S_values
ggplot(pd_data, aes(x = pnodes, y = S)) +
  geom_smooth(color = "purple") +
  labs(title = "pnodes 的生存概率", x = "淋巴结数量", y = "生存概率") +
  theme_minimal()


# partial dependence plots
library(pdp)
pdp_obj <- partial(gbm_model, 
                   pred.var = "pnodes",
                   n.trees = 100)
# Y轴代表预测的“对数风险比”（log hazard ratio）
autoplot(pdp_obj, rug = TRUE, train = gbsg2)


# 单变量部分依赖图 (1D PDP)
# 表示不同水平的淋巴结数量，反映癌症扩散程度。
# 表示从时间 0 到 $ t $ 的总风险水平，值越大表示风险越高。
# 若曲线呈上升趋势，表明 pnodes 增加与风险升高相关；若有拐点或平台，提示非线性效应（如低 pnodes 时风险低，超过某阈值后急升）。
# rug 标记：横轴上的小标记表示训练数据 pnodes 的分布，注意曲线在数据稀疏区域（标记少）的外推风险。
pd_1d <- partial(gbm_model, pred.var = "pnodes", 
                 n.trees = 100, train = GBSG2)
autoplot(pd_1d, rug = TRUE, train = gbsg2) +
  geom_line(color="navy") +
  geom_smooth(color = "lightblue",
              se = FALSE)





## 多变量部分依赖图 (2D PDP)
# 描述：展示两个预测变量（如 pnodes 和 age）的联合边际效应，通常以伪彩色级别图或等高线图形式。
# 特点：chull = TRUE 限制范围在训练数据凸包内，contour = TRUE 添加等高线。

# 颜色/伪彩色：表示对数累积危险（log hazard）的强度，图例（"Log Hazard"）显示颜色与值映射（例如深色为低风险，浅色为高风险）。
# 等高线：contour = TRUE 添加的线表示 log hazard 的等值线，密集线段提示风险梯度大。
# 含义：显示 pnodes 和 age 的联合边际效应，凸显交互效应或独立效应。
# 如何读懂：
# 颜色梯度：观察颜色从深到浅的变化，浅色区域（如 pnodes > 10 和 age > 60）表示高风险。
# 等高线解读：等高线密集处（如 pnodes = 5-10）表明风险快速变化。
# 凸包限制：chull = TRUE 确保仅显示训练数据范围，避免外推。
# 交互效应：若颜色变化主要随 pnodes 而非 age，说明 pnodes 占主导；若有交叉模式，提示两变量交互。
# 临床意义：例如，若 age < 50 时 pnodes 影响小，但 age > 70 时风险随 pnodes 激增，可能反映老年患者更脆弱。
pd_2d <- partial(gbm_model, 
                 pred.var = c("pnodes", "age"), 
                 n.trees = 100, chull = TRUE)
plotPartial(pd_2d, levelplot = TRUE, 
            contour = TRUE, contour.color = "white")
# 或者
autoplot(pd_2d, contour = TRUE, legend.title = "Log Hazard")


# 个体条件期望曲线 (ICE Curves)
# 描述：为每个观测值绘制 pnodes 对预测的个体效应曲线，揭示异质性。
# 方法：
# 使用 partial() 设置 ice = TRUE。
# 使用 plotPartial() 或 autoplot() 绘制。

ice_curves <- partial(gbm_model, 
                      pred.var = "pnodes", 
                      ice = TRUE, 
                      n.trees = 100)
plotPartial(ice_curves, alpha = 0.1, rug = TRUE, train = gbsg2)
# 或者
autoplot(ice_curves, center = FALSE, alpha = 0.1, pdp.color = "red")


# 中心化个体条件期望曲线 (c-ICE Curves)
# 特点：中心化后更易观察个体效应差异。
# 描述：基于 ICE 曲线，中心化后突出个体偏差。
# 方法：
# 使用 partial() 设置 ice = TRUE 和 center = TRUE。
# 使用 plotPartial() 或 autoplot() 绘制。

cice_curves <- partial(gbm_model, pred.var = "pnodes", ice = TRUE, 
                       center = TRUE, n.trees = 100)
autoplot(cice_curves, alpha = 0.1, pdp.color = "blue")


# 带平滑的局部依赖图
# 描述：在部分依赖图上叠加 LOESS 平滑曲线，增强趋势可视化。
# 方法：
# 使用 partial() 计算数据。
# 使用 plotPartial() 或 autoplot() 设置 smooth = TRUE。

pd_smooth <- partial(gbm_model, pred.var = "pnodes", 
                    n.trees = 100, train = gbsg2)
autoplot(pd_smooth, smooth = TRUE)
autoplot(pd_smooth, smooth = TRUE) + 
  geom_line(color="navy") +
  geom_smooth(color = "lightblue",
              se = FALSE)


# pdp package -------------------------------------------------------------

library(pdp)
library(gbm)
library(ggplot2)

# 假设 gbm_model 已训练
data(GBSG2, package = "TH.data")  # 使用 GBSG2 数据
gbm_model <- gbm.fit(x = GBSG2[, c("age", "horTh", "pnodes")], 
                     y = Surv(GBSG2$time, GBSG2$cens), 
                     distribution = "coxph", n.trees = 100, verbose = FALSE)

# 1D PDP
pd_1d <- partial(gbm_model, pred.var = "pnodes", n.trees = 100, train = GBSG2)
autoplot(pd_1d, pdp.color = "purple", rug = TRUE)

# 2D PDP
pd_2d <- partial(gbm_model, pred.var = c("pnodes", "age"), n.trees = 100, chull = TRUE, train = GBSG2)
autoplot(pd_2d, contour = TRUE, legend.title = "Log Hazard")

# ICE Curves
ice <- partial(gbm_model, pred.var = "pnodes", ice = TRUE, n.trees = 100, train = GBSG2)
autoplot(ice, alpha = 0.1, pdp.color = "red")

# c-ICE Curves
cice <- partial(gbm_model, pred.var = "pnodes", ice = TRUE, center = TRUE, n.trees = 100, train = GBSG2)
autoplot(cice, alpha = 0.1, pdp.color = "blue")
