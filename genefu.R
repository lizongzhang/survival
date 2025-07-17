# 基因组学真实数据+随机生存森林分析（gbsg2数据集，来自genefu包）

# 1. 安装&加载必要R包
if(!require(genefu)) install.packages("genefu")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genefu")

if(!require(randomForestSRC)) install.packages("randomForestSRC")
if(!require(survival)) install.packages("survival")
if(!require(survminer)) install.packages("survminer")

library(genefu)
library(randomForestSRC)
library(survival)
library(survminer)

# 2. 加载gbsg2数据集并查看结构
data(package = "genefu")
# 2. 加载数据
data(data.nkis)


# 3. 整理表达矩阵和生存信息
# 基因表达数据是 data.nkis$ge，行为基因，列为样本
# 生存信息 data.nkis$surv[, "t.rfs"] = 无复发生存时间，data.nkis$surv[, "e.rfs"] = 无复发生存结局
expr <- t(data.nkis$ge) # 转为样本为行
surv_info <- data.nkis$surv
# 合并成数据框
dat <- data.frame(time = surv_info[, "t.rfs"],
                  status = surv_info[, "e.rfs"],
                  expr)
# 去除缺失
dat <- na.omit(dat)


# 4. 随机生存森林建模（只选前100个基因做示例，真实分析可选全部）
set.seed(42)
sel_genes <- colnames(expr)[1:100]
dat_rsf <- dat[, c("time", "status", sel_genes)]
rsf_fit <- rfsrc(Surv(time, status) ~ ., data = dat_rsf, ntree = 100, importance = TRUE)

# 5. 查看变量（基因）重要性
imp <- sort(rsf_fit$importance, decreasing = TRUE)
print(imp[1:10])
barplot(imp[1:20], las = 2, col = "skyblue", main = "Top 20基因重要性")

# 6. 风险分组与生存曲线
dat_rsf$risk_score <- rsf_fit$predicted
dat_rsf$risk_group <- cut(dat_rsf$risk_score, breaks = quantile(dat_rsf$risk_score, c(0, 0.5, 1)), 
                          labels = c("Low", "High"), include.lowest = TRUE)
fit <- survfit(Surv(time, status) ~ risk_group, data = dat_rsf)
ggsurvplot(fit, data = dat_rsf, pval = TRUE, risk.table = TRUE, 
           palette = c("#00BFC4", "#F8766D"),
           title = "RSF分组的乳腺癌患者复发生存曲线")