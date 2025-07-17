# 安装包（如未安装）
if(!require(randomForestSRC)) 
  install.packages("randomForestSRC")
if(!require(survival)) 
  install.packages("survival")
if(!require(TCGAbiolinks)) 
  BiocManager::install("TCGAbiolinks")
if(!require(SummarizedExperiment)) 
  BiocManager::install("SummarizedExperiment")

library(randomForestSRC)
library(survival)
library(TCGAbiolinks)
library(SummarizedExperiment)

# 查询并下载表达量数据与生存信息（示例代码，实际下载较慢）
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM")

GDCdownload(query)
expr_data <- GDCprepare(query)

# 提取表达矩阵
expr_matrix <- assay(expr_data)

# 获取生存信息
clinical <- colData(expr_data)
surv_dat <- data.frame(
  sample = clinical$barcode,
  time = clinical$days_to_death,
  status = ifelse(clinical$vital_status == "Dead", 1, 0)
)

# 数据清理：只保留有生存信息的样本，且表达矩阵样本名一致
expr_matrix <- expr_matrix[, colnames(expr_matrix) %in% surv_dat$sample]
surv_dat <- surv_dat[surv_dat$sample %in% colnames(expr_matrix), ]

# 统一顺序
expr_matrix <- expr_matrix[, order(colnames(expr_matrix))]
surv_dat <- surv_dat[order(surv_dat$sample), ]

# 合并数据
rsf_dat <- data.frame(time = surv_dat$time, status = surv_dat$status, t(expr_matrix))
rsf_dat <- rsf_dat[complete.cases(rsf_dat), ]
