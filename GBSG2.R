
library(survival)
library(survminer)


install.packages("TH.data")
library(TH.data)

library(tidyverse)

# Check out the help page for this dataset
help(GBSG2, package = "TH.data")

# Load the data
data(GBSG2, package = "TH.data")

# Look at the summary of the dataset
summary(GBSG2)

table(GBSG2$cens) %>% barplot()



# Create Surv-Object
sobj <- Surv(GBSG2$time, GBSG2$cens)

# Look at 10 first elements
sobj[1:10]

# Look at summary
summary(sobj)

# Look at structure
str(sobj)


# Load the UnempDur data
data(UnempDur, package = "Ecdat")

# Count censored and uncensored data
cens_employ_ft <- table(UnempDur$censor1)
cens_employ_ft

# Create barplot of censored and uncensored data
barplot(cens_employ_ft)

# Create Surv-Object
sobj <- Surv(UnempDur$spell, UnempDur$censor1)

# Look at 10 first elements
sobj[1:10]


# KM ----------------------------------------------------------------------

# Create time and event data
time <- c(5, 6, 2, 4, 4)
event <- c(1, 0, 0, 1, 1)

# Compute Kaplan-Meier estimate
km <- survfit(Surv(time, event) ~ 1)
km

# Take a look at the structure
str(km)

# Create data.frame
data.frame(time = km$time, n.risk = km$n.risk, n.event = km$n.event,
           n.censor = km$n.censor, surv = km$surv)

ggsurvplot(km, conf.int = FALSE, 
           risk.table = "nrisk_cumevents", 
           legend = "none")




library(survival)
library(survminer)

# 创建数据框
df <- data.frame(
  time = c(5, 6, 2, 4, 4),
  event = c(1, 0, 0, 1, 1)
)

# 拟合 KM 曲线，并指定数据
km <- survfit(Surv(time, event) ~ 1, data = df)

# 绘图
ggsurvplot(km, data = df, 
           conf.int = TRUE, 
           risk.table = TRUE)

# breast cancer -----------------------------------------------------------

# Kaplan-Meier estimate
km <- survfit(Surv(time, cens) ~ 1, data = GBSG2)

# Plot of the Kaplan-Meier estimate
ggsurvplot(km)

# Add the risk table to plot
ggsurvplot(km, risk.table = TRUE)

# Kaplan-Meier estimate
km <- survfit(Surv(time, cens) ~ 1, data = GBSG2)

# Plot of the Kaplan-Meier estimate
ggsurvplot(km)

# Add the risk table to plot
ggsurvplot(km, risk.table = TRUE)

# Add a line showing the median survival time
ggsurvplot(km, surv.median.line = "hv", type = "hv")


# Weibull model
wb <- survreg(Surv(time, cens) ~ 1, data = GBSG2)

# Compute the median survival from the model
predict(wb, type = "quantile", p = 0.5, newdata = data.frame(1))

# 70 Percent of patients survive beyond time point...
predict(wb, type = "quantile", p = 0.3, newdata = data.frame(1))


# Weibull model
wb <- survreg(Surv(time, cens) ~ 1, data = GBSG2)

# Retrieve survival curve from model probabilities 
surv <- seq(.99, .01, by = -.01)

# Get time for each probability
t <- predict(wb, type = "quantile", p = 1-surv, newdata = data.frame(1))

# Create data frame with the information
surv_wb <- data.frame(time = t, surv = surv)

# Look at first few lines of the result
head(surv_wb)

# Create data frame with the information needed for ggsurvplot_df
surv_wb <- data.frame(time = t, surv = surv, 
                      upper = NA, lower = NA, std.err = NA)

# Plot
ggsurvplot_df(fit = surv_wb, surv.geom = geom_line)


# Weibull model
wbmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2)
coef(wbmod)

# Retrieve survival curve from model
surv <- seq(.99, .01, by = -.01)
t_yes <- predict(wbmod, type = "quantile", p = 1 - surv,
                 newdata = data.frame(horTh = "yes"))

# Take a look at survival curve
str(t_yes)


# Weibull model
wbmod <- survreg(Surv(time, cens) ~ horTh + tsize, data = GBSG2)

# Imaginary patients
newdat <- expand.grid(
  horTh = levels(GBSG2$horTh),
  tsize = quantile(GBSG2$tsize, probs = c(0.25, 0.5, 0.75)))
newdat

# Compute survival curves
surv <- seq(.99, .01, by = -.01)
t <- predict(wbmod, type = "quantile", p = 1 - surv,
             newdata = newdat)

# How many rows and columns does t have?
dim(t)

# Use cbind() to combine the information in newdat with t
surv_wbmod_wide <- cbind(newdat, t)

# Use melt() to bring the data.frame to long format
surv_wbmod <- melt(surv_wbmod_wide, id.vars = c("horTh", "tsize"), 
                   variable.name = "surv_id", value.name = "time")


# Use surv_wbmod$surv_id to add the correct survival probabilities surv
surv_wbmod$surv <- surv[as.numeric(surv_wbmod$surv_id)]

# Add columns upper, lower, std.err, and strata to the data.frame
surv_wbmod[, c("upper", "lower", "std.err", "strata")] <- NA

# Take a look at the structure of the object
str(surv_wbmod)

# Plot the survival curves


# Plot the survival curves
ggsurvplot_df(surv_wbmod, surv.geom = geom_line,
              linetype = "horTh", color = "tsize", legend.title = NULL)



# Weibull model
wbmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2)

# Log-Normal model
lnmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2, dist = "lognormal")

# Newdata
newdat <- data.frame(horTh = levels(GBSG2$horTh))

# Surv
surv <- seq(.99, .01, by = -.01)

# Survival curve from Weibull model and log-normal model
wbt <- predict(wbmod, type = "quantile", p = 1- surv, newdata = newdat)
lnt <- predict(lnmod, type = "quantile", p = 1- surv, newdata = newdata)

# Melt the data.frame into long format.
surv_long <- melt(surv_wide, id.vars = c("horTh", "dist"),  variable.name = "surv_id", value.name = "time")

# Add column for the survival probabilities
surv_long$surv <- surv[as.numeric(surv_long$surv_id)]

# Add columns upper, lower, std.err, and strata contianing NA values
surv_long[, c("upper", "lower", "std.err", "strata")] <- NA

# Plot the survival curves
ggsurvplot_df(surv_long, surv.geom = geom_line,
              linetype = "horTh", color = "dist" , legend.title = NULL)



# Weibull model 228 cases sex ---------------------------------------------



# Look at the data set
str(dat)

# Estimate a Weibull model
wbmod <- survreg(Surv(time, status) ~ sex, data = dat)
coef(wbmod)



# coxph -------------------------------------------------------------------

# Compute Cox model
cxmod <- coxph(Surv(time, status) ~ performance, data = dat)

# Show model coefficient
coef(cxmod)


# Cox model
cxmod <- coxph(Surv(time, cens) ~ horTh + tsize, data = GBSG2)

# Imaginary patients
newdat <- expand.grid(
  horTh = levels(GBSG2$horTh),
  tsize = quantile(GBSG2$tsize, probs = c(0.25, 0.5, 0.75)))
rownames(newdat) <- letters[1:6]

# Inspect newdat
newdat

# Compute survival curves
cxsf <- survfit(cxmod, data = GBSG2, newdata = newdat, conf.type = "none")

# Look at first 6 rows of cxsf$surv and time points
head(cxsf$surv)
head(cxsf$time)


# Compute data.frame needed for plotting
surv_cxmod0 <- surv_summary(cxsf)

# Look at the first few lines
head(surv_cxmod0)

# Compute data.frame needed for plotting
surv_cxmod0 <- surv_summary(cxsf)

# Look at the first few lines
head(surv_cxmod0)

# Get a character vector of patient letters (patient IDs)
pid <- as.character(surv_cxmod0$strata)

# Multiple of the rows in newdat so that it fits with surv_cxmod0
m_newdat <- newdat[pid, ]

# Add patient info to data.frame
surv_cxmod <- cbind(surv_cxmod0, m_newdat)
head(surv_cxmod)


# Plot
# Plot
ggsurvplot_df(surv_cxmod, linetype = "horTh", color = "tsize",
              legend.title = NULL, censor = FALSE)

# Compute Cox model and survival curves
cxmod <- coxph(Surv(time, status) ~ performance, data = lung)
new_lung <- data.frame(performance = c(60, 70, 80, 90))
cxsf <- survfit(cxmod, data = lung, newdata = new_lung, conf.type = "none")

# Use the summary of cxsf to take a vector of patient IDs
surv_cxmod0 <- surv_summary(cxsf)
pid <- as.character(surv_cxmod0$strata)

# Duplicate rows in newdat to fit with surv_cxmod0 and add them in
m_newdat <- new_lung[pid, , drop = FALSE]
surv_cxmod <- cbind(surv_cxmod0, m_newdat)

# Plot
ggsurvplot_df(surv_cxmod, color = "performance", legend.title = NULL, censor = FALSE)


# Compute Kaplan-Meier curve
km <- survfit(Surv(time, status) ~ 1, data = lung)

# Compute Cox model
cxmod <- coxph(Surv(time, status) ~ performance, data = lung)


# Compute Cox model survival curves
new_lung <- data.frame(performance = c(60, 70, 80, 90))
cxsf <- survfit(cxmod, data = lung, newdata = new_lung, conf.type = "none")

# Plot Kaplan-Meier curve
ggsurvplot(km, conf.int = FALSE)

# Plot Cox model survival curves
ggsurvplot(cxsf, censor = FALSE)

