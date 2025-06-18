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

# status: 0 = censored, 1 = dead

lung <- 
  lung %>% 
  mutate(
    status = recode(status, `1` = 0, `2` = 1))


lung <- 
  lung %>% 
  mutate(
    gender = recode(sex, `1` = "male", `2` = "female"))

lung <- 
  lung %>% 
  mutate(
    female = recode(sex, `1` = 0, `2` = 1))

lung <- lung %>%
  mutate(gender = factor(gender, levels = c("male", "female")))


head(lung$status)

# 字符串转成日期
date_ex <- 
  tibble(
    sx_date = c("2007-06-22", 
                "2004-02-13",
                "2010-10-27"),
    last_fup_date = c("2017-05-15",
                      "2018-07-04",
                      "2016-10-31")
  )

date_ex

date_ex <- 
  date_ex %>% 
  mutate(
    sx_date = ymd(sx_date),
    last_fup_date = ymd(last_fup_date)
  )


# 天数转换成年
date_ex <- 
  date_ex %>% 
  mutate(
    os_yrs = as.duration(sx_date %--% last_fup_date) / dyears(1) )

date_ex

# create surve object

Surv(lung$time, lung$status)

s1 <- survfit(Surv(time, status) ~ 1, data = lung )
str(s1)

survfit2(Surv(time, status) ~ 1, data = lung) %>% 
  ggsurvfit() +
  labs(
    x = "Days", 
    y = "Overall Survival Probability"
  )

survfit2(Surv(time, status) ~ 1, data = lung) %>% 
  ggsurvfit() +
  labs(
    x = "Days", 
    y = "Overall Survival Probability"
  ) +
  add_confidence_interval() +
  add_risktable()


survfit2(Surv(time, status) ~ 1, data = lung) %>% 
  ggsurvfit() +
  labs(
    x = "Days", 
    y = "Overall Survival Probability"
  ) +
  add_confidence_interval(fill = "cyan", alpha = 0.3) +
  add_risktable() +
  scale_color_manual(values = c("cyan4")) 


summary(survfit(Surv(time, status) ~ 1, data = lung),
        times = 365.25)

survfit(Surv(time, status) ~ 1, data = lung) %>% 
  tbl_survfit(
    times = 365.25, 
    label_header = "1-year Survival (95% CI)"
  )


# 不考虑删失数据时的存活时长的中位数226
# 远远小于考虑删失数据的存活时长的中位数310


survfit(Surv(time, status) ~ 1, data = lung)

lung %>% 
  filter( status == 1) %>% 
  summarize(median_surv = median(time))
        
survfit(Surv(time, status) ~ 1, data = lung) %>% 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Median survival (95% CI)**"
  )


# sex role?

# Male=1 Female=2
survdiff(formula = Surv(time, status) ~ sex, data = lung)

survdiff(formula = Surv(time, status) ~ female, data = lung)


survdiff(formula = Surv(time, status) ~ gender, data = lung)

# Male=1 Female=2
coxph(Surv(time, status) ~ sex, data = lung) %>% 
  tbl_regression(exp = TRUE)

coxph(Surv(time, status) ~ sex, data = lung) %>% 
  summary()

coxph(Surv(time, status) ~ female, data = lung) %>% 
  tbl_regression(exp = TRUE)

coxph(Surv(time, status) ~ gender, data = lung) %>% 
  tbl_regression(exp = TRUE,
                 omit_reference_rows = TRUE
                 ) 

coxph(Surv(time, status) ~  female + age + ph.ecog + ph.karno + pat.karno + 
        meal.cal + wt.loss, data = lung) %>% 
  tbl_regression(exp = TRUE)

coxph(Surv(time, status) ~ age + sex + wt.loss, data=lung) %>% 
  summary()


lung.cox <- coxph(Surv(time, status) ~ age + sex + wt.loss, data=lung) 

summary(lung.cox)

library(broom)
lung.cox.tab <- tidy(lung.cox, exponentiate = TRUE, conf.int = TRUE) 

surv.at.means <- survfit(lung.cox) %>%  
  tidy() 

plotdata <- data.frame(age = mean(lung$age), 
                       sex = 1:2,
                       wt.loss = mean(lung$wt.loss, na.rm = T))

surv.by.sex <- survfit(lung.cox, newdata = plotdata)

tidy(surv.by.sex)

plot(surv.by.sex, 
     xlab = "Days", 
     ylab = "Survival Probability",
     col = c("blue", "red"),
     lwd = 2)

library(survminer)
ggsurvplot(surv.by.sex, data = plotdata, censor = F,
           legend.labs = c("male", "female"))




# cox regression estimates ------------------------------------------------
data(cancer, package="survival")

lung

# status: 0 = censored, 1 = dead

lung <- 
  lung %>% 
  mutate(
    status = recode(status, `1` = 0, `2` = 1))

lung <- 
  lung %>% 
  mutate(
    female = recode(sex, `1` = 0, `2` = 1))

model1 <- coxph(Surv(time, status) ~ female + age + wt.loss, 
                data = lung)
model1

model2 <- coxph(Surv(time, status) ~ female + age + wt.loss +
                  ph.ecog + pat.karno, 
                data = lung)
model2

# 拟合模型
model1 <- coxph(Surv(time, status) ~ female + age + wt.loss, data = lung)
model2 <- coxph(Surv(time, status) ~ female + age + wt.loss + ph.ecog + ph.karno + 
                  pat.karno, data = lung)

# 回归表（保留 HR 三位小数）
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

# 合并两张表
tbl_merge(
  tbls = list(tbl1, tbl2),
  tab_spanner = c("**Model 1**", "**Model 2**")
)


# export to word ----------------------------------------------------------

tbl_combined <- tbl_merge(
  tbls = list(tbl1, tbl2),
  tab_spanner = c("**Model 1**", "**Model 2**")
)


library(flextable)
library(officer)


read_docx() %>%
  body_add_par("Table 2. Results of Cox Proportional Hazards Regression Models") %>%
  body_add_flextable(as_flex_table(tbl_combined)) %>%
  print(target = "cox_results.docx")



# check data --------------------------------------------------------------

library(survival)
library(tidyverse)

data(lung, package = "survival")

lung <- lung %>%
  mutate(
    status = recode(status, `1` = 0, `2` = 1),
    female = recode(sex, `1` = 0, `2` = 1)
  )

vars <- c("female", "age", "wt.loss", "meal.cal", "ph.ecog", "ph.karno", "pat.karno")

# 样本总数、事件数、删失数
table(lung$status)

# 缺失值统计
sapply(lung[, vars], \(x) sum(is.na(x)))

# Event per Variable 计算
event_n <- sum(lung$status == 1)
epv <- event_n / length(vars)
epv





