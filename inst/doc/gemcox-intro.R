## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 3.5)

## ----setup--------------------------------------------------------------------
library(gemcox)

## ----data---------------------------------------------------------------------
d <- gemcox_simulate(n = 400, p = 6, mu_sep = 0.5, beta_sep = 2, seed = 42)
dim(d$X)
mean(d$status)          # proportion of subjects with an event

## ----fit----------------------------------------------------------------------
fit <- gemcox(d$X, time = d$time, status = d$status, K = 2)
fit

## ----tau----------------------------------------------------------------------
head(round(fit$tau, 3))
hist(apply(fit$tau, 1, max), breaks = 20, xlim = c(0.5, 1),
     main = "Largest membership weight per subject", xlab = "max tau")

## ----coef---------------------------------------------------------------------
round(coef(fit), 3)
contrast_hat  <- coef(fit)[, 1] - coef(fit)[, 2]
contrast_true <- d$beta_list[[1]] - d$beta_list[[2]]
## direction of the estimated contrast vs the truth (sign is arbitrary)
abs(sum(contrast_hat * contrast_true)) /
  sqrt(sum(contrast_hat^2) * sum(contrast_true^2))

## ----summary------------------------------------------------------------------
summary(fit)

## ----test---------------------------------------------------------------------
ht <- gemcox_heterogeneity_test(fit, B = 19, seed = 1)
ht

## ----predict------------------------------------------------------------------
new <- gemcox_simulate(n = 3, p = 6, mu_sep = 0.5, beta_sep = 2, seed = 7)
predict(fit, new$X, type = "tau")
predict(fit, new$X, type = "survival", times = c(50, 150, 300))

