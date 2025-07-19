# 准备：安装和加载包 ---------------------------------------------------------------

# 需要的包列表
pkgs <- c("survival", "ggsurvfit", "survminer", "rms", 
          "tidycmprsk", "cmprsk", 
          "tidyverse", "broom", "glue", "gtsummary", 
          "flextable", "officer", "showtext", 
          "RColorBrewer", "MASS")

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

# 6.3.1 数据背景 ---------------------------------------------------------------

# 加载 MASS 包中的 Melanoma 数据集
data(Melanoma, package = "MASS")

?MASS::Melanoma


# 6.3.2 数据准备 ---------------------------------------------------------------

# Melanoma 中的status：1 表示因黑色素瘤死亡，2 表示存活，3 表示因其他原因死亡
# 重编码 status（将2存活改为0，1因黑色素瘤死亡保持为1，3因其他原因死亡改为2）
# 并转换为因子型
Melanoma <- 
  Melanoma %>% 
  mutate(status = as.factor(recode(status, 
                                   `2` = 0, 
                                   `1` = 1, 
                                   `3` = 2)
                            )
  )

head(Melanoma) # 查看数据前6行，检查数据结构



# ### 6.3.3 计算CIF ---------------------------------------------------------

# 计算目标事件/竞争事件的累积发生概率（CIF）
cif <- tidycmprsk::cuminc(Surv(time, status) ~ 1, data = Melanoma)
cif

# 按ulcer分组，计算并展示5年（1826.25天）时点的CIF表，并进行组间比较
tidycmprsk::cuminc(Surv(time, status) ~ ulcer, data = Melanoma) %>% 
  tidycmprsk::tbl_cuminc(
    times = 1826.25, 
    label_header = "**{time/365.25}-year cuminc**") %>% 
  add_p()


# 6.3.4 CIF曲线 ---------------------------------------------------------

## survminer::ggcompetingrisks -----------------------------------------------------


# 比较各种结局的累积发生率随时间的变化
survminer::ggcompetingrisks(
  fit = survfit(Surv(time, status) ~ 1, data = Melanoma),
  data = Melanoma,
  palette = "jco"
)

# 查看所有调色板名称及类型
library(RColorBrewer)
display.brewer.all()

# 使用palette = "Set3"
ggcompetingrisks(
  fit = survfit(Surv(time, status) ~ 1, data = Melanoma),
  data = Melanoma,
  palette = "Set3"
)

ggcompetingrisks(
  fit = survfit(Surv(time, status) ~ 1, data = Melanoma),
  data = Melanoma,
  palette = "Paired"
)


## ggsurvfit::ggcuminc -----------------------------------------------------

#绘制 Failure type "1"（也就是 type 1，对应 status==1 的事件）的累积发生概率
tidycmprsk::cuminc(Surv(time, status) ~ 1, data = Melanoma) %>% 
  ggsurvfit::ggcuminc(col = "#B3D4FC") + 
  labs(x = "Days") + 
  add_confidence_interval() + # 添加置信区间
  add_risktable() # 风险表

# 绘制 Failure type "2"的累积发生概率（CIF）曲线
tidycmprsk::cuminc(Surv(time, status) ~ 1, data = Melanoma) %>% 
  ggcuminc(outcome = "2", color = "#D7B9D5") + 
  labs(x = "Days") + 
  #add_confidence_interval() +
  add_risktable()

# 多种结局（outcome=1,2）的CIF曲线
tidycmprsk::cuminc(Surv(time, status) ~ 1, data = Melanoma) %>% 
  ggcuminc(outcome = c("1", "2")) +
  ylim(c(0, 1)) + 
  labs(x = "Days")

tidycmprsk::cuminc(Surv(time, status) ~ 1, data = Melanoma) %>% 
  ggcuminc(outcome = c("1", "2")) +
  aes(color = outcome) +       # 强制美学映射
  scale_color_manual(values = c("1" = "#B3D4FC", 
                                "2" = "#F4C2C2"),
                     labels = c("death form Melanoma", 
                                "death from other cause")) +
  ylim(0, 1) +
  labs(x = "Days") +
  add_risktable() +
  scale_ggsurvfit()

# 按ulcer分组，绘制多结局CIF曲线，添加置信区间和风险表
tidycmprsk::cuminc(Surv(time, status) ~ ulcer, data = Melanoma) %>% 
  ggcuminc(outcome = c("1", "2")) + 
  labs(x = "Days") + 
  add_confidence_interval() +
  add_risktable()


# 6.3.5 灰色检验Gray Test -----------------------------------------------------

# Gray检验：比较两个或多个组在特定结局（某一竞争风险事件）上的累积发生率是否存在显著差异。
# p<0.05, 存在显著差异

tidycmprsk::cuminc(Surv(time, status) ~ ulcer, data = Melanoma) %>% 
  tidycmprsk::tbl_cuminc(
    times = c(365.25, 1826.25, 3652.5),  # 1年、5年、10年
    label_header = "**{time/365.25}-year cuminc**"
  ) %>% 
  add_p(
    pvalue_fun = function(x) style_pvalue(x, digits = 3)
  )



# 6.3.6 估计Fine-Gray 模型 ----------------------------------------------------


## 估计Fine-Gray模型 -----------------------------------------------------------

# 竞争风险回归，分析 sex 和 age 对结局的影响
tidycmprsk::crr(Surv(time, status) ~ sex + age, data = Melanoma)

## 表格化输出 -----------------------------------------------------------

tidycmprsk::crr(Surv(time, status) ~ sex + age, data = Melanoma) %>% 
  tbl_regression(exp = TRUE) # exp=TRUE 输出HR

## 多个模型表格化输出 -----------------------------------------------------------


fit1 <- tidycmprsk::crr(Surv(time, status) ~ sex + age, data = Melanoma)
fit2 <- tidycmprsk::crr(Surv(time, status) ~ sex + age  + year + thickness, data = Melanoma)
fit3 <- tidycmprsk::crr(Surv(time, status) ~ sex + age  + year + thickness + ulcer, data = Melanoma)

library(gtsummary)

tbl1 <- tbl_regression(fit1, exponentiate = TRUE)
tbl2 <- tbl_regression(fit2, exponentiate = TRUE)
tbl3 <- tbl_regression(fit3, exponentiate = TRUE)

# 合并报告
merged_tbl <- tbl_merge(
  tbls = list(tbl1, tbl2, tbl3),
  tab_spanner = c("**Model 1**", "**Model 2**", "**Model 3**")
)

merged_tbl


## 将表格输出到WORD文档 ------------------------------------------------------------

library(flextable)
library(officer)

# 将表格转换为 flextable 对象
ft <- as_flex_table(merged_tbl)

# 创建Word文档并写入表格
doc <- read_docx() %>%
  body_add_flextable(ft)

# 保存Word文件
print(doc, target = "FineGray_models.docx")



# 6.3.7 系数森林图 -------------------------------------------------------------

fit <- tidycmprsk::crr(Surv(time, status) ~ sex + age + 
                         year + thickness + ulcer, data = Melanoma)
tbl <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
names(tbl)

ggplot(tbl, aes(y = term, x = estimate)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  xlab("Subdistribution Hazard Ratio (HR)") +
  ylab("") +
  theme_minimal()


# 6.3.8 对比竞争风险模型和Cox模型 ----------------------------------------------------

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
