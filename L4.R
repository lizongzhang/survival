install.packages("survival")
install.packages("ggsurvfit")
install.packages("gtsummary")
install.packages("tidyverse")

library(survival)
library(ggsurvfit)
library(gtsummary)
library(tidyverse)

data(cancer, package="survival")

lung


# 定性变量的编码 -----------------------------------------------------------------

lung <- lung %>%
  mutate(
    status = recode(status, `1` = 0, `2` = 1),
    female = recode(sex, `1` = 0, `2` = 1)
  )


# 数据检查--------------------------------------------------------------------

vars <- c("female", "age", "wt.loss", "meal.cal", "ph.ecog", "ph.karno", "pat.karno")

# 样本总数、事件数、删失数
table(lung$status)

# 缺失值统计
sapply(lung[, vars], \(x) sum(is.na(x)))

# Event per Variable 计算
event_n <- sum(lung$status == 1)
epv <- event_n / length(vars)
epv


# cox模型的估计 ----------------------------------------------------------------

library(gtsummary)

coxph(Surv(time, status) ~ female, data = lung)

coxph(Surv(time, status) ~ female, data = lung) %>% 
  summary()

coxph(Surv(time, status) ~ female, data = lung) %>% 
  tbl_regression(exp = TRUE)

coxph(Surv(time, status) ~ female, data = lung) %>% 
  tbl_regression(exp = FALSE)


coxph(Surv(time, status) ~ female + age + wt.loss, data = lung) %>% 
  summary()

coxph(Surv(time, status) ~ female + age + wt.loss, data = lung) %>% 
  tbl_regression(exp = TRUE)

coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.karno + 
        pat.karno, data = lung) %>% 
  summary()

# 表格化输出----------------------------------------------------------------

# 拟合模型
model1 <- coxph(Surv(time, status) ~ female + age + wt.loss, data = lung)


model2 <- coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.karno + 
                  pat.karno, data = lung)

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


# 导出表格到WORD -------------------------------------------------------------

library(flextable)
library(officer)


read_docx() %>%
  body_add_par("Table 2. Results of Cox Proportional Hazards Regression Models") %>%
  body_add_flextable(as_flex_table(tbl_combined)) %>%
  print(target = "cox_results.docx")
