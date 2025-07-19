install.packages("survival")
install.packages("ggsurvfit")
install.packages("gtsummary")
install.packages("tidyverse")

# 准备：安装和加载包 ---------------------------------------------------------------

# 需要的包列表
pkgs <- c("survival", "ggsurvfit", "survminer", "rms", 
          "tidyverse", "broom", "glue", "gtsummary", 
          "flextable", "officer", "showtext", "RColorBrewer")

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


# 4.4. 2 数据准备 --------------------------------------------------------------------

## 加载示例数据, cancer列表中的lung数据框
data(cancer, package="survival")

## 定性变量的编码 

lung <- lung %>%
  mutate(
    status = recode(status, `1` = 0, `2` = 1),    # 事件状态重编码：1=生存(0), 2=死亡(1)
    female = recode(sex, `1` = 0, `2` = 1)        # 性别重编码：1=男(0), 2=女(1)
  )


# 4.4.3 数据检查--------------------------------------------------------------------

vars <- c("female", "age", "wt.loss", "meal.cal", "ph.ecog", "ph.karno", "pat.karno")

# 样本总数、事件数、删失数
table(lung$status)

# 缺失值统计
sapply(lung[, vars], \(x) sum(is.na(x)))

# Event per Variable 计算
event_n <- sum(lung$status == 1)
epv <- event_n / length(vars)
epv


# 4.4.4 cox模型的估计 ----------------------------------------------------------------

fit <- coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.ecog, 
      data = lung)

summary(fit)

# 4.4.5 检验风险比例假定 ----------------------------------------------------------------

fit <- coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.karno, 
             data = lung)

# Schoenfeld残差检
cox.zph(fit)

# Schoenfeld残差检验的可视化

fit %>% 
  cox.zph() %>% 
  plot()

fit %>% 
  cox.zph() %>% 
  plot(var = "female")

fit %>% 
  cox.zph() %>% 
  plot(var = "ph.ecog")






# 4.4.6 表格化输出----------------------------------------------------------------

model1 <- coxph(Surv(time, status) ~ female + age + wt.loss, data = lung)


model2 <- coxph(Surv(time, status) ~ female + age + wt.loss +ph.ecog + ph.karno, 
                data = lung)

tbl1 <- tbl_regression(
  model1,
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3)
)

tbl2 <- tbl_regression(
  model2,
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3)
)

tbl_merge(
  tbls = list(tbl1, tbl2),
  tab_spanner = c("**Model 1**", "**Model 2**")
)

tbl_combined <- tbl_merge(
  tbls = list(tbl1, tbl2),
  tab_spanner = c("**Model 1**", "**Model 2**")
)


# 4.4.7 导出表格到WORD -------------------------------------------------------------

library(flextable)
library(officer)


read_docx() %>%
  body_add_par("Table 2. Results of Cox Proportional Hazards Regression Models") %>%
  body_add_flextable(as_flex_table(tbl_combined)) %>%
  print(target = "cox_results.docx")



# 4.5.1 绘制基线生存曲线 ----------------------------------------------------------------

library(survival)
library(ggsurvfit)

# 拟合Cox模型
fit <- coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.karno, 
             data = lung)

# 拟合基线生存曲线
sf <- survfit(fit)

# 用 ggsurvfit 绘图
ggsurvfit(sf) +
  labs(
    title = "Cox模型基线生存曲线",
    x = "随访时间（天）",
    y = "生存概率"
  ) +
  theme_minimal()



# 4.5.2 绘制不同人群特征组合下的预测生存曲线 ------------------------------------------------


newdata <- data.frame(
  female   = c(0, 1, 1),
  age      = c(60, 60, 70),
  wt.loss  = c(0, 0, 5),
  ph.ecog  = c(1, 1, 2),
  ph.karno = c(80, 80, 60)
)

sf <- survfit(fit, newdata = newdata)

library(survminer)
ggsurvplot(
  sf,
  data = lung,
  conf.int = FALSE,
  legend.labs = c("男性, 60岁, wt.loss =0, ph.ecog=1, ph.karno=80",
                  "女性, 60岁, wt.loss =0, ph.ecog=1, ph.karno=80",
                  "女性, 70岁, wt.loss =5, ECOG=2,ph.ecog=2,ph.karno=80" ),
  palette = "Set2", # 使用RColorBrewer的 Set2 调色板
  title = "不同变量组合下的预测生存曲线",
  xlab = "时间（天）",
  ylab = "生存概率"
)


# 4.5.3 Cox回归系数森林图 --------------------------------------------------------

coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.karno, 
      data = lung) %>% 
  ggforest(main = "Cox回归系数森林图", fontsize = 1.2,
           cpositions = c(0.1, 0.22, 0.4), refLabel = "对照", 
           noDigits = 3, data = lung)


# 4.5.4 列线图（Nomogram） -----------------------------------------------------------


#install.packages("rms")
library(rms)
dd <- datadist(lung)
options(datadist="dd")

# 用cph拟合
fit <- cph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.karno, 
           data = lung, x = TRUE, y = TRUE, surv = TRUE)

surv <- Survival(fit)


nom <- nomogram(fit, fun=list(function(x) surv(365, x)), funlabel="1年生存概率")

plot(nom)






# 4.5.5 高低风险组生存曲线----------------------------------------------------------------



library(survival)
library(survminer)

# 步骤1: 计算风险评分
risk_score <- predict(fit, type = "lp")  # "lp" = linear predictor

# 步骤2: 分组（以中位数为界）
lung$risk_group <- ifelse(risk_score > median(risk_score, na.rm = TRUE), 
                          "High risk", "Low risk")

# 步骤3: 拟合分组生存曲线
fit_group <- survfit(Surv(time, status) ~ risk_group, 
                     data = lung)

# 步骤4: 绘制生存曲线
ggsurvplot(
  fit_group,
  data = lung,
  pval = TRUE,                  # 显示log-rank检验p值
  legend.title = "Risk Group",
  legend.labs = c("High risk", "Low risk"),
  palette = c("#E41A1C", "#377EB8"),  # 红、蓝
  xlab = "时间（天）",
  ylab = "生存概率",
  title = "高低风险组生存曲线"
)




