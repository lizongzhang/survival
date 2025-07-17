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

# 数据编码 ---------------------------------------------------------------

# 加载 MASS 包中的 Melanoma 数据集
data(Melanoma, package = "MASS")

# Melanoma 中的status：1 表示因黑色素瘤死亡，2 表示存活，3 表示因其他原因死亡
# 重编码 status（将2存活改为0，1因黑色素瘤死亡保持为1，3因其他原因死亡改为2），并转换为因子型
Melanoma <- 
  Melanoma %>% 
  mutate(
    status = as.factor(recode(status, `2` = 0, `1` = 1, `3` = 2))
  )


head(Melanoma) # 查看数据前6行，检查数据结构


# 事件发生堆积条形图等，展示各事件累积发生的比例。
library(survminer)
library(survival)
ggcompetingrisks(
  fit = survfit(Surv(time, status) ~ 1, data = Melanoma),
  data = Melanoma,
  palette = "jco"
)



# 计算无分组的竞争风险累积发生概率（CIF）
cif <- tidycmprsk::cuminc(Surv(time, status) ~ 1, data = Melanoma)
cif


# 绘制无分组的CIF曲线，添加置信区间和风险表

# 绘制 Failure type "1"（也就是 type 1，对应 status==1 的事件）的累积发生概率（CIF）曲线

tidycmprsk::cuminc(Surv(time, status) ~ 1, data = Melanoma) %>% 
  ggcuminc() + 
  labs(
    x = "Days"
  ) + 
  add_confidence_interval() +
  add_risktable()


# 绘制 Failure type "2"（也就是 type 2，对应 status==2 的事件）的累积发生概率（CIF）曲线
tidycmprsk::cuminc(Surv(time, status) ~ 1, data = Melanoma) %>% 
  ggcuminc(outcome = "2") + 
  labs(
    x = "Days"
  ) + 
  add_confidence_interval() +
  add_risktable()


# 同时展示多种结局（outcome=1,2）的CIF曲线，设置y轴范围
tidycmprsk::cuminc(Surv(time, status) ~ 1, data = Melanoma) %>% 
  ggcuminc(outcome = c("1", "2")) +
  ylim(c(0, 1)) + 
  labs(
    x = "Days"
  )

# 按ulcer分组，计算并展示5年（1826.25天）时点的CIF表，并进行组间比较
tidycmprsk::cuminc(Surv(time, status) ~ ulcer, data = Melanoma) %>% 
  tbl_cuminc(
    times = 1826.25, 
    label_header = "**{time/365.25}-year cuminc**") %>% 
  add_p()

# 按ulcer分组，绘制CIF曲线，添加置信区间和风险表
tidycmprsk::cuminc(Surv(time, status) ~ ulcer, data = Melanoma) %>% 
  ggcuminc() + 
  labs(
    x = "Days"
  ) + 
  add_confidence_interval() +
  add_risktable()


#  Gray 检验
#使用 tbl_cuminc 和 add_p 进行 Gray 检验
tidycmprsk::cuminc(Surv(time, status) ~ ulcer, data = Melanoma) %>% 
  tidycmprsk::tbl_cuminc(
    times = c(365.25, 1826.25, 3652.5),  # 1年、5年、10年
    label_header = "**{time/365.25}-year cuminc**"
  ) %>% 
  add_p(
    pvalue_fun = function(x) style_pvalue(x, digits = 3)  # 用 function 语法
  )



# 竞争风险回归，分析 sex 和 age 对结局的影响
tidycmprsk::crr(Surv(time, status) ~ sex + age, data = Melanoma)

# 竞争风险回归结果整理为回归表（exp=TRUE 输出HR及置信区间）
tidycmprsk::crr(Surv(time, status) ~ sex + age, data = Melanoma) %>% 
  tbl_regression(exp = TRUE)

# COX回归，仅以 status==1 为结局，分析 sex 和 age 的影响，并整理输出
coxph(
  Surv(time, ifelse(status == 1, 1, 0)) ~ sex + age, 
  data = Melanoma
) %>% 
  tbl_regression(exp = TRUE)



# 对比Fine-Gray模型和coxph ---------------------------------------------------



library(tidycmprsk)
library(gtsummary)
library(survival)

# 竞争风险回归（Fine-Gray模型）
fit_fg <- tidycmprsk::crr(Surv(time, status) ~ sex + age, data = Melanoma)
tbl_fg <- tbl_regression(fit_fg, exp = TRUE, label = list(sex = "Sex", age = "Age")) %>%
  modify_header(label = "**变量**") %>%
  modify_caption("**Fine-Gray 竞争风险回归**")

# Cox回归
fit_cox <- coxph(Surv(time, ifelse(status == 1, 1, 0)) ~ sex + age, data = Melanoma)
tbl_cox <- tbl_regression(fit_cox, exp = TRUE, label = list(sex = "Sex", age = "Age")) %>%
  modify_header(label = "**变量**") %>%
  modify_caption("**Cox 比例风险模型**")

tbl_merge(
  tbls = list(tbl_fg, tbl_cox),
  tab_spanner = c("**Fine-Gray 竞争风险回归**", "**Cox 比例风险模型**")
)



# 拟合 Fine-Gray 模型
crr_model <- tidycmprsk::crr(Surv(time, status) ~ sex + age, data = Melanoma)

# 查看模型摘要
crr_model %>%
  tidy(conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  print()

# 准备个体数据（示例：三个个体）
# 准备个体数据（包含 time 和 status 列）
newdata <- data.frame(
  sex = c("Male", "Female", "Male"),  # 性别
  age = c(50, 60, 70),               # 年龄
  time = c(0, 0, 0),                 # 占位符，值不影响预测
  status = c(0, 0, 0)                # 占位符，值不影响预测
)

# 指定预测时间点（1000 和 2000 天）
times <- seq(500, 3000, 100)

# 预测个体 CIF
cif_pred <- predict(crr_model, newdata = newdata, time = times)

# 调试：检查 cif_pred 结构
str(cif_pred)
print(cif_pred)

# 将列表格式的 cif_pred 转换为数据框
cif_data <- tibble(
  Individual = rep(paste("Individual", 1:3, "(", newdata$sex, ", Age ", newdata$age, ")", sep = ""), each = length(times)),
  time = rep(times, times = 3),
  CIF = unlist(cif_pred)
)

# 查看整理后的数据
print(cif_data)



# 可视化 CIF
ggplot(cif_data, aes(x = time, y = CIF, color = Individual)) +
  geom_line() +
  geom_point() +
  labs(title = "个体 CIF 预测（黑色素瘤死亡）",
       x = "时间（天）", y = "累积发生概率",
       color = "个体") +
  theme_minimal() +
  ylim(0, 0.4)  # 限制 y 轴范围以更好地显示




