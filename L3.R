

# 准备：安装和加载包 ---------------------------------------------------------------

# 需要的包列表
pkgs <- c("survival", "ggsurvfit", "survminer", 
          "rms", "survRM2",
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


# 3.2 一个简单的例子 -------------------------------------------------------------

# 创建数据框
df <- data.frame(
  time = c(5, 6, 2, 4, 4),
  event = c(1, 0, 0, 1, 1)
)

survival::Surv(df$time, df$event)

# 拟合 KM 曲线，并指定数据
km <- survival::survfit(Surv(time, event) ~ 1, data = df)

# list列表
km

km$surv

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
survminer::ggsurvplot(km, data = df, 
           palette = "blue",
           conf.int = TRUE, 
           risk.table = TRUE)


# 3.3.1 绘制KM曲线 ------------------------------------------------------------

## 变量重新编码-----------------------------------------------------------------------

## 加载数据, cancer列表list中的lung数据框
library(tidyverse)

data(cancer, package="survival")

?lung

lung <- lung %>%
  mutate(
    status = recode(status, "1" = 0, "2" = 1),   # 删失 = 0, 死亡 = 1
    female = recode(sex, "1" = 0, "2" = 1))    #   男 = 0, 女 = 1


# 第2步：计算生存概率 --------------------------------------------------------------

# 波浪号（~）后的 1 表示拟合一个不含分层变量的Kaplan-Meier模型

km <- survfit(Surv(time, status) ~ 1, 
              data = lung)

# 显示对象的内部结构
str(km)

# 对生存分析结果km进行数据整洁化，并显示前6行
km %>%         # Kaplan-Meier生存分析对象
  broom::tidy() %>%   # 转换为整洁数据框格式（需要broom包）
  head()       # 显示前6行结果

km %>% 
  tidy() %>% 
  tail() # 显示后6行结果


# 第3步：绘制KM曲线  -----------------------------------------------------------------------

# 绘制KM曲线
library(survminer)

survminer::ggsurvplot(
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

# 第4步：预测1年期生存概率  -----------------------------------------------------------------------
# 1年 = 365天 
summary(survfit(Surv(time, status) ~ 1, data = lung), 
        times = 365)

# 1年 = 365.25天 考虑闰年效应
summary(survfit(Surv(time, status) ~ 1, data = lung), 
        times = 365.25)


#  不考虑删失和考虑删失对KM曲线的影响 -----------------------------------------------------------------------

# 不考虑删失的KM曲线
km_wrong <- survfit(Surv(time) ~ 1, lung)

# 考虑删失的KM曲线
km <- survfit(Surv(time, status) ~ 1, lung)

# Plot the two and compare
ggsurvplot_combine(list(correct = km, wrong = km_wrong))

# 3.3.2 生存概率的表格输出-----------------------------------------------------------------------
# 表格化输出结果
survfit(Surv(time, status) ~ 1, data = lung) %>% 
  gtsummary::tbl_survfit(
    times = 365.25,
    label_header = "**1-year survival (95% CI)**"
  )

survfit(Surv(time, status) ~ sex, data = lung) %>% 
  gtsummary::tbl_survfit(times = c(180, 365, 730),  # 半年、1年、2年
                         label_header = "**Survival at {time} days (95% CI)**")


#   3.3.3 处理删失数据的可视化-----------------------------------------------------------------------
survminer::ggsurvplot(km, data = lung, 
                      censor = TRUE, 
                      censor.shape = "|", 
                      censor.size = 4)

survminer::ggsurvplot(km, data = lung, 
                      censor = TRUE, 
                      censor.shape = "$", 
                      censor.size = 6)


#   3.3.4 分面（faceting）展示多组KM曲线-----------------------------------------------------------------------

# 使用lung数据集按性别进行Log-rank检验
km_by_sex <- survfit(Surv(time, status) ~ female, 
                     data = lung)

# 示例：按性别和治疗类型分面
survminer::ggsurvplot_facet(
  km_by_sex, 
  data = lung, 
  facet.by = "ph.ecog",  # 按ECOG状态分面
  palette = "jco",
  conf.int = TRUE
)


# 使用lung数据集按性别进行Log-rank检验
km_by_sex <- survfit(Surv(time, status) ~ female, data = lung)

# 示例：按性别和治疗类型分面
survminer::ggsurvplot_facet(
  km_by_sex, 
  data = lung, 
  facet.by = "ph.ecog",  # 按ECOG状态分面
  palette = "jco",
  conf.int = TRUE
)




#   3.3.5 交互式KM曲线（plotly）-----------------------------------------------------------------------

library(plotly)
km <- survfit(Surv(time, status) ~ 1, data = lung)
p <- ggsurvplot(km, data = lung)
plotly::ggplotly(p$plot)

#   3.3.6 patchwork：用于组合多个KM曲线图-----------------------------------------------------------------------
library(patchwork)
p1 <- ggsurvplot(survfit(Surv(time, status) ~ sex, 
                         data = lung), risk.table = TRUE)

# 先筛选出ph.ecog为0/1/2的数据
lung_sub <- subset(lung, ph.ecog %in% c(0, 1, 2))

p2 <- ggsurvplot(survfit(Surv(time, status) ~ ph.ecog, 
                         data = lung_sub), 
                 risk.table = TRUE)

p1$plot + p2$plot  # 使用patchwork组合图形


# 3.3.7 Log-rank检验 ------------------------------------------------------------

# 两组log-rank 检验
survival::survdiff(Surv(time, status) ~ female, 
                   data = lung)


# 先筛选出ph.ecog为0/1/2的数据
lung_sub <- subset(lung, ph.ecog %in% c(0, 1, 2))

# 多组log-rank 检验
survival::survdiff(Surv(time, status) ~ ph.ecog, 
                   data = lung_sub)

# 使用lung数据集按性别进行Log-rank检验
km_by_sex <- survfit(Surv(time, status) ~ sex, data = lung)

# 使用survminer绘制KM曲线并添加Log-rank p值
survminer::ggsurvplot(km_by_sex, 
                      data = lung, 
                      pval = TRUE,  # 显示Log-rank p值
                      pval.method = TRUE,  # 显示检验方法
                      palette = c("blue", "red"),
                      risk.table = TRUE)



# 3.4 `ggsurvfit`-----------------------------------------------------------------------


# 全样本KM曲线 -----------------------------------------------------------------

ggsurvfit::survfit2(Surv(time, status) ~ 1, data = lung) %>% 
  ggsurvfit::ggsurvfit(color = "navy") +
  theme_bw() +
  add_risktable(size = 3)


#   绘制男性/女性的KM曲线-----------------------------------------------------------------------
survfit2(Surv(time, status) ~ female, data = lung) %>% 
  ggsurvfit() +
  theme_bw() +
  add_risktable(size = 3)

#   按ph.ecog分组的KM曲线-----------------------------------------------------------------------
# ph.ecog 是患者的ECOG体能状态评分
# 0代表完全正常
# 1代表体力活动受限但能行走和从事轻工作
# 2代表能自理，但不能工作，白天一半以上时间能活动
# 3仅能部分自理，白天一半以上时间卧床或坐轮椅
# 4完全丧失生活自理能力，完全卧床或坐轮椅

# 查看ph.ecog的分布
table(lung$ph.ecog)


## 绘制ph.ecog = 0/1/2三个组别的KM曲线
# . place holder占位符
lung %>%
  filter(ph.ecog %in% c(0, 1, 2)) %>%             # 只保留ph.ecog为0、1、2的数据
  survfit2(Surv(time, status) ~ ph.ecog, data = .) %>%
  ggsurvfit() +
  theme_bw()

# 自定义线型和颜色-----------------------------------------------------------------------
# ggplot2 package
survfit2(Surv(time, status) ~ female, data = lung) %>% 
  ggsurvfit(linetype_aes = TRUE) +       #aesthetic                     # 启用分组线型
  scale_color_manual(values = c("0" = "#b8c6e6",              # 男性淡蓝紫
                                "1" =  "#ffb6b9"),
                     labels = c("0" = "Male",              # 男性淡蓝紫
                                "1" =  "Female")) +        # 女性淡粉
  scale_linetype_manual(values = c("0" = "solid",          # 男性实线
                                   "1" = "dashed"),
                        labels = c("0" = "Male",              # 男性淡蓝紫
                                   "1" =  "Female")) +    # 女性虚线
  theme_bw() +                                                # 黑白主题
  theme(legend.position = "bottom") 


#   使用期刊绘图模版-----------------------------------------------------------------------

library(ggsci)

survfit2(Surv(time, status) ~ female, data = lung) %>%
  ggsurvfit(linewidth = 1) +                        # 绘制KM生存曲线，线宽为1
  add_confidence_interval() +                       # 添加置信区间
  add_risktable() +                                 # 添加风险表（在险人数）
  add_quantile(y_value = 0.5, 
               color = "gray50", 
               linewidth = 0.75) + # 添加中位生存线
  scale_color_jama() +                              # JAMA期刊风格线条配色
  scale_fill_jama() +                               # JAMA期刊风格填充色（置信区间）
  theme_pubr() +                                      # pubr主题
  theme(legend.position = "bottom")                 # 图例置于底部


#   3.5 `survRM2`-----------------------------------------------------------------------

library(survRM2)
survRM2::rmst2(time = lung$time, 
               status = lung$status, 
               arm = lung$female)

library(survRM2)
survRM2::rmst2(time = lung$time, 
               status = lung$status, 
               arm = lung$female,
               tau = 365)

