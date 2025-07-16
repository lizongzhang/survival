
# 准备：安装和加载包 ---------------------------------------------------------------


# 需要的包列表
pkgs <- c("survival", "survminer", "tidyverse",
          "broom", "glue", "gtsummary")

# 检查并安装未安装的包
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# 加载所有包
lapply(pkgs, library, character.only = TRUE)




# 3.2 一个简单的例子 -------------------------------------------------------------

# 创建数据框
df <- data.frame(
  time = c(5, 6, 2, 4, 4),
  event = c(1, 0, 0, 1, 1)
)

# 创建生存对象
Surv(time, event)

# 拟合 KM 曲线，并指定数据
km <- survfit(Surv(time, event) ~ 1, data = df)

# 指定要显示的时刻
times_to_show <- 1:6

# 提取指定时刻的n.risk, n.event, survival prob
km_sum <- summary(km, times = times_to_show)

# 整理为数据框
km_table <- data.frame(
  time     = km_sum$time,
  n.risk   = km_sum$n.risk,
  n.event  = km_sum$n.event,
  survival = km_sum$surv,
  std.err  = km_sum$std.err,
  lower95  = km_sum$lower,
  upper95  = km_sum$upper
)

# 显示生存概率
km_table


# 绘制Kaplan-Meier曲线
ggsurvplot(km, data = df, 
           palette = "blue",
           conf.int = TRUE, 
           risk.table = TRUE)



# 3.3 绘制Kaplan-Meier曲线方法一：survminer::ggsurvplot()  ------------------------


## 数据准备 --------------------------------------------------------------------

## 加载示例数据, cancer列表中的lung数据框
data(cancer, package="survival")

## 定性变量的编码 
library(tidyverse)

lung <- lung %>%
  mutate(
    status = recode(status, `1` = 0, `2` = 1),    # 事件状态重编码：1=生存(0), 2=死亡(1)
    female = recode(sex, `1` = 0, `2` = 1)        # 性别重编码：1=男(0), 2=女(1)
  )

lung <- lung %>% 
  mutate(
    female = factor(female, levels = c(0, 1), 
                    labels = c("Male", "Female"))  # 性别变量转为因子factor
  )

## 计算生存概率-----------------------------------------------------------------

# 样本总数、事件数、删失数
table(lung$status)

# 计算生存概率
km <- survfit(Surv(time, status) ~ 1, data = lung)

# 显示对象的内部结构
str(km)

# 对生存分析结果km进行数据整洁化，并显示前6行
km %>%         # Kaplan-Meier生存分析对象
  tidy() %>%   # 转换为整洁数据框格式（需要broom包）
  head()       # 显示前6行结果

km %>% 
  tidy() %>% 
  tail() # 显示后6行结果


## 预测1年期生存概率 ------------------------------------------------------------

# 1年 = 365.25天
summary(survfit(Surv(time, status) ~ 1, data = lung), times = 365.25)

# 表格化输出结果
survfit(Surv(time, status) ~ 1, data = lung) %>% 
  gtsummary::tbl_survfit(
    times = 365.25,
    label_header = "**1-year survival (95% CI)**"
  )

## 绘制KM曲线-----------------------------------------------------------------

library(survminer)

ggsurvplot(
  km,                # KM拟合对象
  data = lung,       # 数据集
  linetype = 1,      # 实线
  palette = "jco",   # 颜色方案
  conf.int = FALSE,  # 不显示置信区间
  # risk.table = TRUE,     # 风险表（如需可取消注释）
  cumevents = TRUE,       # 显示累计事件
  # cumcensor = TRUE,      # 累计删失（如需可取消注释）
  tables.height = 0.2     # 表格高度
)

## 计算中位数生存时间 ------------------------------------------------------------

survfit(Surv(time, status) ~ 1, data = lung)

ggsurvplot(
  km,                # KM拟合对象
  data = lung,       # 数据集
  linetype = 1,      # 实线
  palette = "jco",   # 颜色方案
  conf.int = FALSE,  # 不显示置信区间
  surv.median.line = "hv", # 显示中位生存水平/垂直线
  # risk.table = TRUE,     # 风险表（如需可取消注释）
  cumevents = TRUE,       # 显示累计事件
  # cumcensor = TRUE,      # 累计删失（如需可取消注释）
  tables.height = 0.2     # 表格高度
)


## 不考虑删失和考虑删失对KM曲线的影响 ------------------------------------------------------



# Estimate the survivor function pretending that all censored observations are actual observations.
km_wrong <- survfit(Surv(time) ~ 1, lung)

# Estimate the survivor function from this dataset via kaplan-meier.
km <- survfit(Surv(time, status) ~ 1, lung)

# Plot the two and compare
ggsurvplot_combine(list(correct = km, wrong = km_wrong))


# 3.4 绘制Kaplan-Meier曲线方法二: ggsurvfit::ggsurvfit()  ---------------------


## 绘制全样本的KM曲线 --------------------------------------------------------------

library(ggsurvfit)

survfit2(Surv(time, status) ~ 1, data = lung) %>% 
  ggsurvfit(color = "navy") +
  theme_bw() +
  add_risktable(size = 3)

## 绘制不同组别的KM曲线 --------------------------------------------------------------

## 绘制男性/女性的KM曲线
survfit2(Surv(time, status) ~ female, data = lung) %>% 
  ggsurvfit() +
  theme_bw() +
  add_risktable(size = 3)


# ph.ecog 是患者的ECOG体能状态评分
# 0代表完全正常
# 1代表体力活动受限但能行走和从事轻工作
# 2代表能自理，但不能工作，白天一半以上时间能活动
# 3仅能部分自理，白天一半以上时间卧床或坐轮椅
# 4完全丧失生活自理能力，完全卧床或坐轮椅

# 查看ph.ecog的分布
table(lung$ph.ecog)


## 绘制ph.ecog = 0/1/2三个组别的KM曲线
lung %>%
  filter(ph.ecog %in% c(0, 1, 2)) %>%             # 只保留ph.ecog为0、1、2的数据
  survfit2(Surv(time, status) ~ ph.ecog, data = .) %>%
  ggsurvfit() +
  theme_bw()




## 自定义线型和颜色 --------------------------------------------------------------

survfit2(Surv(time, status) ~ female, data = lung) %>% 
  ggsurvfit(linetype_aes = TRUE) +                            # 启用分组线型
  scale_color_manual(values = c("Male" = "#b8c6e6",              # 男性淡蓝紫
                                "Female" =  "#ffb6b9")) +        # 女性淡粉
  scale_linetype_manual(values = c("Male" = "solid",          # 男性实线
                                   "Female" = "dashed")) +    # 女性虚线
  theme_bw() +                                                # 黑白主题
  theme(legend.position = "bottom")                           # 图例置底


## 使用期刊绘图模版 --------------------------------------------------------------

library(ggsci)

survfit2(Surv(time, status) ~ female, data = lung) %>%
  ggsurvfit(linewidth = 1) +                        # 绘制KM生存曲线，线宽为1
  add_confidence_interval() +                       # 添加置信区间
  add_risktable() +                                 # 添加风险表（在险人数）
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) + # 添加中位生存线
  scale_color_jama() +                              # JAMA期刊风格线条配色
  scale_fill_jama() +                               # JAMA期刊风格填充色（置信区间）
  theme_pubr() +                                      # pubr主题
  theme(legend.position = "bottom")                 # 图例置于底部


# ggsci的期刊模版

# "npg"      # Nature Publishing Group（自然出版集团）
# "aaas"     # American Association for the Advancement of Science / Science（美国科学促进会/《科学》杂志）
# "nejm"     # New England Journal of Medicine（新英格兰医学杂志）
# "lancet"   # The Lancet（柳叶刀杂志）
# "jama"     # Journal of the American Medical Association（美国医学会杂志）
# "jco"      # Journal of Clinical Oncology（临床肿瘤学杂志）
# "ucscgb"   # UCSC Genome Bioinformatics（加州大学圣克鲁兹基因组生物信息学）
# "d3"       # D3.js配色（著名数据可视化库）
# "igv"      # Integrative Genomics Viewer配色（综合基因组浏览器）
# "uchicago" # University of Chicago配色（芝加哥大学）
# "startrek" # Star Trek配色（星际迷航主题配色）

# 用法只需替换scale_color_xxx()和scale_fill_xxx()即可。




# 3.5 log-rank检验 --------------------------------------------------------------

# log-rank 检验两组或多组人的生存曲线是否存在显著差异。
# 如果 P 值 < 0.05，说明两组生存曲线差异显著
# 如果 P 值 > 0.05，说明差异步显著


## 两个组别的log-rank检验 ---------------------------------------------------------

survdiff(Surv(time, status) ~ female, data = lung_sub)


## 多个组别的log-rank检验 ---------------------------------------------------------


library(survival)

# 先筛选出ph.ecog为0/1/2的数据
lung_sub <- subset(lung, ph.ecog %in% c(0, 1, 2))

# log-rank 检验
survdiff(Surv(time, status) ~ ph.ecog, data = lung_sub)



## 在KM曲线图像中添加log-rank的结果 ---------------------------------------------------


## 绘制男性/女性的KM曲线并添加log-rank检验的结果

sf <- survfit2(Surv(time, status) ~ female, data = lung)

sf %>%
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  scale_ggsurvfit() +
  labs(caption = glue::glue("Log-rank {survfit2_p(sf)}"))

## 绘制ph.ecog = 0/1/2 的KM曲线并添加log-rank检验的结果

sf <- survfit2(Surv(time, status) ~ ph.ecog, data = lung_sub)

sf %>%
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  scale_ggsurvfit() +
  labs(caption = glue::glue("Log-rank {survfit2_p(sf)}"))

## 绘制age的四个等分组别的KM曲线

lung <- lung %>%
  mutate(age_group = cut_number(age, 4))

sf <- survfit2(Surv(time, status) ~ age_group, data = lung)

sf %>%
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  theme_minimal() +
  scale_ggsurvfit() +
  labs(caption = glue::glue("Log-rank {survfit2_p(sf)}"))


# 模拟数据
time <- c(5, 6, 2, 4, 4)          # 生存时间
event <- c(1, 0, 0, 1, 1)         # 1 = 事件发生, 0 = 删失

n <- length(time)
y_pos <- rev(1:n)                 # 从上到下排列个体编号（1在最上）

# 创建空白图形
plot(NULL, xlim = c(0, max(time)), ylim = c(0.5, n + 0.5),
     xlab = "t", ylab = "individual", yaxt = "n", 
     col.axis = "blue", col.lab = "blue")

axis(2, at = y_pos, labels = 1:n, col.axis = "blue")  # 添加个体编号

# 设置边框颜色为蓝色
box(col = "blue")

# 绘制生存线段
segments(x0 = 0, y0 = y_pos, x1 = time, y1 = y_pos, lwd = 2,
         col = "blue")

# 添加标记：事件为×，删失为●
for (i in 1:n) {
  if (event[i] == 1) {
    points(time[i], y_pos[i], pch = 4, cex = 1.5, lwd = 2, col = "blue")  # ×
  } else {
    points(time[i], y_pos[i], pch = 16, cex = 1.2, col = "blue")          # ●
  }
}


# 创建数据框
df <- data.frame(
  time = c(5, 6, 2, 4, 4),
  event = c(1, 0, 0, 1, 1)
)

# 拟合 KM 曲线，并指定数据
km <- survfit(Surv(time, event) ~ 1, data = df)

fit_sum <- summary(km)

df_surv <- data.frame(
  time = fit_sum$time,
  n.risk = fit_sum$n.risk,
  n.event = fit_sum$n.event,
  survival = fit_sum$surv,
  std.err = fit_sum$std.err,
  lower = fit_sum$lower,
  upper = fit_sum$upper
)

View(df_surv)


# 绘图
ggsurvplot(km, data = df, 
           palette = "blue",
           conf.int = TRUE, 
           risk.table = TRUE)

ggsurvplot(km, data = df, 
           linetype = 2, 
           palette = "blue",
           conf.int = FALSE, 
           surv.median.line = "hv",
           risk.table = TRUE, 
           cumevents = TRUE, 
           cumcensor = TRUE,
           tables.height  = 0.2)


data(cancer, package="survival")

km <- survfit(Surv(time, status) ~ 1, data = lung)

ggsurvplot(km, data = lung, 
           linetype = 2, 
           palette = "jco",
           conf.int = FALSE, 
           surv.median.line = "hv",
           #risk.table = TRUE, 
           cumevents = TRUE, 
           #cumcensor = TRUE,
           tables.height  = 0.2)
























# example: lung -----------------------------------------------------------

data(cancer, package="survival")


# 定性变量的编码 -----------------------------------------------------------------
library(tidyverse)

lung <- lung %>%
  mutate(
    status = recode(status, `1` = 0, `2` = 1),
    female = recode(sex, `1` = 0, `2` = 1)
  )

lung <- lung %>% 
  mutate(
    female = factor(female, levels = c(0, 1), 
                    labels = c("Male", "Female")),
   )
  

survfit2(Surv(time, status) ~ 1, data = lung) %>% 
  ggsurvfit() +
  theme_bw() +
  theme(legend.position = "bottom") +
  add_risktable(size = 3)

survfit2(Surv(time, status) ~ female, data = lung) %>% 
  ggsurvfit() +
  theme_bw() +
  theme(legend.position = "bottom") +
  add_risktable(size = 3)

survfit2(Surv(time, status) ~ female, data = lung) %>% 
  ggsurvfit(, linetype_aes = TRUE) +
  scale_color_manual(values = c("Male" = "red", "Female" = "blue")) +
  scale_linetype_manual(values = c("Male" = "solid", 
                                   "Female" = "dashed")) +
  theme_bw() +
  theme(legend.position = "bottom") 

survfit2(Surv(time, status) ~ female, data = lung) %>% 
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  add_risktable() +
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_ggsurvfit() 

library(survival)
library(ggsurvfit)
library(ggsci)      # 期刊配色
library(ggpubr)     # 主题

survfit2(Surv(time, status) ~ female, data = lung) %>%
  ggsurvfit(linewidth = 1) +
  scale_color_npg() +                  # Nature配色
  theme_pubr() +                       # 论文风格主题
  theme(legend.position = "bottom")

library(ggsurvfit)
library(ggsci)

survfit2(Surv(time, status) ~ female, data = lung) %>%
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  add_risktable() +
  add_quantile(y_value = 0.5, color = "gray50", linewidth = 0.75) +
  scale_color_jama() +    
  scale_fill_jama() +     
  theme_bw() +
  theme(legend.position = "bottom")
  

sf <- survfit2(Surv(time, status) ~ female, data = lung)

sf %>%
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  scale_ggsurvfit() +
  labs(caption = glue::glue("Log-rank {survfit2_p(sf)}"))






survfit2(Surv(time, status) ~ female, data = lung) %>%
  ggsurvfit(type = "cloglog") +
  scale_x_continuous(transform = "log")



# wrong_km ----------------------------------------------------------------


# Create dancedat data
dancedat <- data.frame(
  name = c("Chris", "Martin", "Conny", "Desi", "Reni", "Phil", 
           "Flo", "Andrea", "Isaac", "Dayra", "Caspar"),
  time = c(20, 2, 14, 22, 3, 7, 4, 15, 25, 17, 12),
  obs_end = c(1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0))

# Estimate the survivor function pretending that all censored observations are actual observations.
km_wrong <- survfit(Surv(time) ~ 1, data = dancedat)

# Estimate the survivor function from this dataset via kaplan-meier.
km <- survfit(Surv(time, obs_end) ~ 1, data = dancedat)

# Plot the two and compare
ggsurvplot_combine(list(correct = km, wrong = km_wrong))

# Estimate the survivor function pretending that all censored observations are actual observations.
km_wrong <- survfit(Surv(time) ~ 1, lung)

# Estimate the survivor function from this dataset via kaplan-meier.
km <- survfit(Surv(time, status) ~ 1, lung)

# Plot the two and compare
ggsurvplot_combine(list(correct = km, wrong = km_wrong))




# breast cancer -----------------------------------------------------------

# Kaplan-Meier estimate
km <- survfit(Surv(time, cens) ~ 1, data = GBSG2)

# Plot of the Kaplan-Meier estimate
ggsurvplot(km)

# Add the risk table to plot
ggsurvplot(km, risk.table = TRUE)

# Kaplan-Meier estimate
km <- survfit(Surv(time, cens) ~ 1, data = GBSG2)

# Plot of the Kaplan-Meier estimate
ggsurvplot(km)

# Add the risk table to plot
ggsurvplot(km, risk.table = TRUE)

# Add a line showing the median survival time
ggsurvplot(km, surv.median.line = "hv", type = "hv")


# Weibull model
wb <- survreg(Surv(time, cens) ~ 1, data = GBSG2)

# Compute the median survival from the model
predict(wb, type = "quantile", p = 0.5, newdata = data.frame(1))

# 70 Percent of patients survive beyond time point...
predict(wb, type = "quantile", p = 0.3, newdata = data.frame(1))


# Weibull model
wb <- survreg(Surv(time, cens) ~ 1, data = GBSG2)

# Retrieve survival curve from model probabilities 
surv <- seq(.99, .01, by = -.01)

# Get time for each probability
t <- predict(wb, type = "quantile", p = 1-surv, newdata = data.frame(1))

# Create data frame with the information
surv_wb <- data.frame(time = t, surv = surv)

# Look at first few lines of the result
head(surv_wb)

# Create data frame with the information needed for ggsurvplot_df
surv_wb <- data.frame(time = t, surv = surv, 
                      upper = NA, lower = NA, std.err = NA)

# Plot
ggsurvplot_df(fit = surv_wb, surv.geom = geom_line)


# Weibull model
wbmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2)
coef(wbmod)

# Retrieve survival curve from model
surv <- seq(.99, .01, by = -.01)
t_yes <- predict(wbmod, type = "quantile", p = 1 - surv,
                 newdata = data.frame(horTh = "yes"))

# Take a look at survival curve
str(t_yes)


# Weibull model
wbmod <- survreg(Surv(time, cens) ~ horTh + tsize, data = GBSG2)

# Imaginary patients
newdat <- expand.grid(
  horTh = levels(GBSG2$horTh),
  tsize = quantile(GBSG2$tsize, probs = c(0.25, 0.5, 0.75)))
newdat

# Compute survival curves
surv <- seq(.99, .01, by = -.01)
t <- predict(wbmod, type = "quantile", p = 1 - surv,
             newdata = newdat)

# How many rows and columns does t have?
dim(t)

# Use cbind() to combine the information in newdat with t
surv_wbmod_wide <- cbind(newdat, t)

# Use melt() to bring the data.frame to long format
surv_wbmod <- melt(surv_wbmod_wide, id.vars = c("horTh", "tsize"), 
                   variable.name = "surv_id", value.name = "time")


# Use surv_wbmod$surv_id to add the correct survival probabilities surv
surv_wbmod$surv <- surv[as.numeric(surv_wbmod$surv_id)]

# Add columns upper, lower, std.err, and strata to the data.frame
surv_wbmod[, c("upper", "lower", "std.err", "strata")] <- NA

# Take a look at the structure of the object
str(surv_wbmod)

# Plot the survival curves


# Plot the survival curves
ggsurvplot_df(surv_wbmod, surv.geom = geom_line,
              linetype = "horTh", color = "tsize", legend.title = NULL)



# Weibull model
wbmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2)

# Log-Normal model
lnmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2, dist = "lognormal")

# Newdata
newdat <- data.frame(horTh = levels(GBSG2$horTh))

# Surv
surv <- seq(.99, .01, by = -.01)

# Survival curve from Weibull model and log-normal model
wbt <- predict(wbmod, type = "quantile", p = 1- surv, newdata = newdat)
lnt <- predict(lnmod, type = "quantile", p = 1- surv, newdata = newdata)

# Melt the data.frame into long format.
surv_long <- melt(surv_wide, id.vars = c("horTh", "dist"),  variable.name = "surv_id", value.name = "time")

# Add column for the survival probabilities
surv_long$surv <- surv[as.numeric(surv_long$surv_id)]

# Add columns upper, lower, std.err, and strata contianing NA values
surv_long[, c("upper", "lower", "std.err", "strata")] <- NA

# Plot the survival curves
ggsurvplot_df(surv_long, surv.geom = geom_line,
              linetype = "horTh", color = "dist" , legend.title = NULL)



# Weibull model 228 cases sex ---------------------------------------------



# Look at the data set
str(dat)

# Estimate a Weibull model
wbmod <- survreg(Surv(time, status) ~ sex, data = dat)
coef(wbmod)



# coxph -------------------------------------------------------------------

# Compute Cox model
cxmod <- coxph(Surv(time, status) ~ performance, data = dat)

# Show model coefficient
coef(cxmod)


# Cox model
cxmod <- coxph(Surv(time, cens) ~ horTh + tsize, data = GBSG2)

# Imaginary patients
newdat <- expand.grid(
  horTh = levels(GBSG2$horTh),
  tsize = quantile(GBSG2$tsize, probs = c(0.25, 0.5, 0.75)))
rownames(newdat) <- letters[1:6]

# Inspect newdat
newdat

# Compute survival curves
cxsf <- survfit(cxmod, data = GBSG2, newdata = newdat, conf.type = "none")

# Look at first 6 rows of cxsf$surv and time points
head(cxsf$surv)
head(cxsf$time)


# Compute data.frame needed for plotting
surv_cxmod0 <- surv_summary(cxsf)

# Look at the first few lines
head(surv_cxmod0)

# Compute data.frame needed for plotting
surv_cxmod0 <- surv_summary(cxsf)

# Look at the first few lines
head(surv_cxmod0)

# Get a character vector of patient letters (patient IDs)
pid <- as.character(surv_cxmod0$strata)

# Multiple of the rows in newdat so that it fits with surv_cxmod0
m_newdat <- newdat[pid, ]

# Add patient info to data.frame
surv_cxmod <- cbind(surv_cxmod0, m_newdat)
head(surv_cxmod)


# Plot
# Plot
ggsurvplot_df(surv_cxmod, linetype = "horTh", color = "tsize",
              legend.title = NULL, censor = FALSE)

# Compute Cox model and survival curves
cxmod <- coxph(Surv(time, status) ~ performance, data = lung)
new_lung <- data.frame(performance = c(60, 70, 80, 90))
cxsf <- survfit(cxmod, data = lung, newdata = new_lung, conf.type = "none")

# Use the summary of cxsf to take a vector of patient IDs
surv_cxmod0 <- surv_summary(cxsf)
pid <- as.character(surv_cxmod0$strata)

# Duplicate rows in newdat to fit with surv_cxmod0 and add them in
m_newdat <- new_lung[pid, , drop = FALSE]
surv_cxmod <- cbind(surv_cxmod0, m_newdat)

# Plot
ggsurvplot_df(surv_cxmod, color = "performance", legend.title = NULL, censor = FALSE)


# Compute Kaplan-Meier curve
km <- survfit(Surv(time, status) ~ 1, data = lung)

# Compute Cox model
cxmod <- coxph(Surv(time, status) ~ performance, data = lung)


# Compute Cox model survival curves
new_lung <- data.frame(performance = c(60, 70, 80, 90))
cxsf <- survfit(cxmod, data = lung, newdata = new_lung, conf.type = "none")

# Plot Kaplan-Meier curve
ggsurvplot(km, conf.int = FALSE)

# Plot Cox model survival curves
ggsurvplot(cxsf, censor = FALSE)