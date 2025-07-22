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


# 4.3.2 变量编码 --------------------------------------------------------------------

## 加载数据, cancer列表list中的lung数据框
data(cancer, package="survival")
lung <- lung %>%
  mutate(
    status = recode(status, "1" = 0, "2" = 1),   # 删失 = 0, 死亡 = 1
    female = recode(sex, "1" = 0, "2" = 1))  

# 4.3.3 数据检查--------------------------------------------------------------------

vars <- c("female", "age", "wt.loss", 
          "meal.cal", "ph.ecog", "ph.karno", "pat.karno")

# 样本总数、事件数、删失数
table(lung$status)

# 缺失值统计
sapply(lung[, vars], \(x) sum(is.na(x)))

# Event per Variable 
event_n <- sum(lung$status == 1)
epv <- event_n / length(vars)
epv

# 4.3.4 Cox模型的估计 ----------------------------------------------------------------

fit <- coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.ecog, 
             data = lung)

summary(fit)

coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.karno, 
      data = lung) %>%
  summary()

# 4.3.5 检验风险比例假定 ----------------------------------------------------------------

fit <- coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.karno, 
             data = lung)

# Schoenfeld残差检
cox.zph(fit)

# Schoenfeld残差检验的可视化

# 平滑曲线接近水平，无明显上升或下降趋势, 风险比例假设成立。

fit %>% 
  cox.zph() %>% 
  plot()

fit %>% 
  cox.zph() %>% 
  plot(var = "female")

fit %>% 
  cox.zph() %>% 
  plot(var = "ph.ecog")


# 4.3.6 表格化输出----------------------------------------------------------------

model1 <- coxph(Surv(time, status) ~ female + age + wt.loss, data = lung)


model2 <- coxph(Surv(time, status) ~ female + age + wt.loss +ph.ecog + ph.karno, 
                data = lung)


library(gtsummary)
tbl1 <- tbl_regression(
  model1,
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3)
)
tbl1

tbl2 <- tbl_regression(
  model2,
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3)
)
tbl2

tbl_merge(
  tbls = list(tbl1, tbl2),
  tab_spanner = c("**Model 1**", "**Model 2**")
)

tbl_combined <- tbl_merge(
  tbls = list(tbl1, tbl2),
  tab_spanner = c("**Model 1**", "**Model 2**")
)


# 4.3.7 导出表格到WORD -------------------------------------------------------------

library(flextable)
library(officer)


read_docx() %>%
  body_add_par("Table 2. Results of Cox Proportional Hazards Regression Models") %>%
  body_add_flextable(as_flex_table(tbl_combined)) %>%
  print(target = "cox_results.docx")



# 4.4.1 绘制生存曲线 ----------------------------------------------------------------

library(survival)
library(ggsurvfit)

# Cox 模型拟合
fit <- coxph(Surv(time, status) ~ female + age + wt.loss + 
               ph.ecog + ph.karno, 
             data = lung)



## 绘制基线生存曲线 ----------------------------------------------------------------

# 定义基线个体（所有covariate协/解释变量设为 0 ）
newdf <- data.frame(
  female = 0,
  age = 0,  
  wt.loss = 0,
  ph.ecog = 0,
  ph.karno = 0
)

# 获取基线生存曲线
sf_base <- survfit(fit, newdata = newdf)

str(sf_base)


# 绘制基线生存曲线
library(ggsurvfit)
ggsurvfit(sf_base,
          col = "blue") +
  labs(
    title = "Cox模型基线生存曲线",
    x = "时间（天）",
    y = "生存概率"
  ) +
  theme_minimal()



## 绘制平均个体生存曲线 ----------------------------------------------------------------

library(survival)
library(ggsurvfit)

# average person: female = 0, mean(age),  mean(wt.loss) ,mean(ph.ecog), mean(ph.karno)
# survfit(fit) 默认生成的生存曲线是基于协变量的平均值（对于连续变量）或参考水平（对于分类变量）计算的

sf_ave_ind <- survfit(fit)

# 用 ggsurvfit 绘图
ggsurvfit::ggsurvfit(sf_ave_ind, 
          col = "red") +
  labs(title = "Cox模型基线生存曲线",
       x = "时间（天）",
       y = "生存概率") +
  theme_minimal()

## 合并两条生存曲线ggsurvplot_combine() ----------------------------------------------------------------
library(survminer)
ggsurvplot_combine(
  list(average = sf_ave_ind, baseline = sf_base),
  data = lung,
  palette = c("red", "blue"),  # 为两条曲线指定不同颜色
  linetype = c("solid", "dashed"),  # 为两条曲线指定不同线型
  legend.title = "曲线类型",
  legend.labs = c("平均个体生存曲线", "基线生存曲线"),
  title = "Cox模型生存曲线比较",
  xlab = "时间（天）",
  ylab = "生存概率"
) 


## 特定个体的生存曲线 ----------------------------------------------------------------

newdata <- data.frame(
  female   = c(0, 1, 1),
  age      = c(60, 60, 70),
  wt.loss  = c(0, 0, 5),
  ph.ecog  = c(1, 1, 2),
  ph.karno = c(80, 80, 60)
)

sf_3ind <- survfit(fit, newdata = newdata)

library(survminer)
survminer::ggsurvplot(
  sf_3ind,
  data = lung,
  conf.int = FALSE,
  legend.labs = c("男性, 60岁, wt.loss =0, ph.ecog=1, ph.karno=80",
                  "女性, 60岁, wt.loss =0, ph.ecog=1, ph.karno=80",
                  "女性, 70岁, wt.loss =5, ph.ecog=2,ph.karno=80" ),
  palette = "Set2", # 使用RColorBrewer的 Set2 调色板
  title = "不同变量组合下的预测生存曲线",
  xlab = "时间（天）",
  ylab = "生存概率"
)

# 男性和女性生存曲线对比

sf <- survfit(fit, newdata = data.frame(
  female = c(0,1),
  age = mean(lung$age, na.rm = TRUE),
  wt.loss = mean(lung$wt.loss, na.rm = TRUE),
  ph.ecog = mean(lung$ph.ecog, na.rm = TRUE),
  ph.karno = mean(lung$ph.karno, na.rm = TRUE)
))

sf <- survfit(fit, newdata = data.frame(
  female = c(0,0),
  age = mean(lung$age, na.rm = TRUE),
  wt.loss = mean(lung$wt.loss, na.rm = TRUE),
  ph.ecog = c(0,1,2,3),
  ph.karno = mean(lung$ph.karno, na.rm = TRUE)
))

ggsurvplot(sf, 
           data = lung, 
           legend.labs = c("0", "1", "2", "3"),
           conf.int = FALSE
           )

## 高低风险组生存曲线 --------------------------------------------------------

library(survival)
library(survminer)

# 步骤1: 计算风险评分
fit <- coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog,
             data = lung,
             na.action = na.exclude)  # 保留原始行号结构

lung$risk_score <- predict(fit, type = "lp", na.action = na.exclude)

# 步骤2: 分组（以中位数为界）
lung$risk_group <- ifelse(lung$risk_score > median(lung$risk_score, 
                                                   na.rm = TRUE), 
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


## 调整后生存曲线 --------------------------------------------------------

survminer::ggadjustedcurves(fit, 
                 variable = "age", 
                 data = lung)

survminer::ggadjustedcurves(fit, 
                            variable = "female", 
                            data = lung)

survminer::ggadjustedcurves(fit, 
                            variable = "ph.ecog", 
                            data = lung)


survminer::ggadjustedcurves(fit, 
                            variable = "wt.loss", 
                            data = lung)



# 4.4.2 累积风险曲线 --------------------------------------------------------

fit <- coxph(Surv(time, status) ~ female + age + wt.loss + 
               ph.ecog + ph.karno, data = lung) 

survfit(fit) %>% 
  str()

survfit(fit)$cumhaz

survfit(fit) %>% 
  ggsurvfit(type = "cumhaz") 

survfit(fit) %>% 
  ggsurvfit() 

survfit(fit) %>% 
  ggsurvfit(type = "risk") 


# 4.4.3 Cox回归系数森林图 --------------------------------------------------------


# Cox 模型拟合
fit <- coxph(Surv(time, status) ~ female + age + wt.loss + 
               ph.ecog + ph.karno, 
             data = lung)


survminer::ggforest(fit, 
                    main = "Cox回归系数森林图", 
                    fontsize = 1.2,
                    cpositions = c(0.1, 0.22, 0.4), 
                    refLabel = "对照", 
                    noDigits = 3, data = lung)

# 4.4.4 互补对数-对数图cloglog plot --------------------------------------------------------
# 生存曲线拟合（Kaplan-Meier）
fit_km <- survfit(Surv(time, status) ~ sex, 
                   data = lung)

# 绘制 cloglog 图：用于检验比例风险假设（PH Assumption）
ggsurvfit(fit_km, type = "cloglog") +
  scale_x_continuous(trans = "log") +
  labs(
    title = "Cloglog Plot for Proportional Hazards Assumption",
    x = "log(Time)", 
    y = "log(-log(Survival Probability))"
  ) +
  theme_minimal()

# 4.4.5 列线图rms::nomogram() --------------------------------------------------------
library(rms)
dd <- datadist(lung) #创建描述数据分布的信息
options(datadist="dd") #告诉rms包刚刚创建的数据分布信息

# 用cph拟合
fit <- cph(Surv(time, status) ~ female + age + wt.loss +
             ph.ecog + ph.karno, data = lung, 
           x = TRUE, y = TRUE, surv = TRUE)

#创建了输入“时间”和“变量值”，输出“生存概率”的函数surv
surv <- rms::Survival(fit)

nom <- rms::nomogram(fit, 
                fun=list(function(x) surv(365, x)), 
                funlabel="1年生存概率")

plot(nom)

# 4.4.6 校准曲线rms::calibrate() --------------------------------------------------------
library(rms)

dd <- datadist(lung)  #创建数据分布对象以处理连续变量
options(datadist = 'dd')  # 设置默认数据分布

fit_rms <- cph(Surv(time, status) ~ female + age + wt.loss + 
                 ph.ecog + ph.karno, 
               data = lung, x=TRUE, y=TRUE, surv=TRUE)  # 拟合Cox模型

cal <- calibrate(fit_rms, 
                 cmethod="KM", 
                 method="boot", u=365, m=50, B=100)  # 校准模型，计算365天生存概率
plot(cal)  # 绘制校准曲线




