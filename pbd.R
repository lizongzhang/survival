# 加载包
library(randomForestSRC)
library(survival)
library(survminer)

# 1. 查看并准备数据
data(pbc, package = "randomForestSRC")
str(pbc)
head(pbc)

psych::describe(pbc)
pbc_noNA <- na.omit(pbc)

# 去除缺失的结局变量
pbc2 <- na.omit(pbc[, c("time", "status", "age", "edema", "bili", "albumin", 
                        "alk.phos", "ast", "trt")])


# 检查以下变量名都在pbc中
vars <- c("time", "status", "age", "edema", "bili", "albumin", "alk.phos", "ast", "trt")
vars[!vars %in% colnames(pbc)]
# 这行会显示哪些变量不存在

# 选存在的列
pbc2 <- pbc[, vars[vars %in% colnames(pbc)]]

# 去除缺失
pbc2 <- na.omit(pbc2)

# 2. 构建随机生存森林模型
set.seed(123)
rsf_fit <- rfsrc(Surv(days, status) ~ ., data = pbc, ntree = 100, importance = TRUE)

library(randomForestSRC)
rsf_fit <- rfsrc(Surv(days, status) ~ ., data = pbc, ntree = 100, 
                 na.action = "na.impute",
                 importance = TRUE)

# 3. 变量重要性图
var_imp <- sort(rsf_fit$importance, decreasing = TRUE)
barplot(var_imp, main = "变量重要性", las = 2, col = "steelblue")

# 4. 计算风险分数并分组
risk_score <- rsf_fit$predicted
pbc$risk_group <- cut(risk_score, breaks = quantile(risk_score, c(0, 0.5, 1)),
                       labels = c("Low", "High"), include.lowest = TRUE)

# 5. 绘制分组生存曲线
fit <- survfit(Surv(days, status) ~ risk_group, data = pbc)
library(survminer)
ggsurvplot(fit, data = pbc, pval = TRUE, risk.table = TRUE,
           palette = c("#00BFC4", "#F8766D"),
           title = "随机生存森林分组的PBC患者生存曲线")
days()