# 查看R的版本
R.version.string

# 查看RStudio的版本
RStudio.Version()$version



# 安装单个包
install.packages("survival")
install.packages("survminer")
install.packages("flexsurv")
install.packages("cmprsk")

# 一次性安装全部
install.packages(c("survival", "survminer", "flexsurv", "cmprsk"))

# 查看包的最新版的版本号
options(repos = c(CRAN = "https://cloud.r-project.org"))  # 使用主镜像
available.packages()["survminer", "Version"]


install.packages("writexl")  # 如果未安装，先安装
library(writexl)
write_xlsx(lung, path = "lung.xlsx")

install.packages("TH.data")
library(TH.data)
data(GBSG2)
write_xlsx(GBSG2, path = "breast.xlsx")


