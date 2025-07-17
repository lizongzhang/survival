# Weibull model 基线生存曲线 --------------------------------------------


library(TH.data)
data("GBSG2")

# Weibull model
wb <- survreg(Surv(time, cens) ~ 1, data = GBSG2)

# Compute the median survival from the model
predict(wb, type = "quantile", p = 0.5, newdata = data.frame(1))

# 70 Percent of patients survive beyond time point...
predict(wb, type = "quantile", p = 0.3, newdata = data.frame(1))


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
library(survminer)
ggsurvplot_df(fit = surv_wb, surv.geom = geom_line)


# # Weibull模型单个解释变量 ---------------------------------------------------------
# 
# # Weibull model
# wbmod <- survreg(Surv(time, cens) ~ horTh, data = GBSG2)
# coef(wbmod)
# 
# # Retrieve survival curve from model
# surv <- seq(.99, .01, by = -.01)
# t_yes <- predict(wbmod, type = "quantile", p = 1 - surv,
#                  newdata = data.frame(horTh = "yes"))
# 
# # Take a look at survival curve
# str(t_yes)

# Weibull模型：两个解释变量 ---------------------------------------------------------

wbmod <- survreg(Surv(time, cens) ~ horTh + tsize, data = GBSG2)

summary(wbmod)

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
library(reshape2)
surv_wbmod <- melt(surv_wbmod_wide, id.vars = c("horTh", "tsize"), 
                   variable.name = "surv_id", value.name = "time")


# Use surv_wbmod$surv_id to add the correct survival probabilities surv
surv_wbmod$surv <- surv[as.numeric(surv_wbmod$surv_id)]

# Add columns upper, lower, std.err, and strata to the data.frame
surv_wbmod[, c("upper", "lower", "std.err", "strata")] <- NA

# Take a look at the structure of the object
str(surv_wbmod)

# Plot the survival curves
ggsurvplot_df(surv_wbmod, surv.geom = geom_line,
              linetype = "horTh", color = "tsize", 
              legend.title = NULL)