## Feature names must survive unchanged (known defect 3: make.names turned
## "IgG3 AMA-1 MFI | V20" into "IgG3.AMA.1.MFI...V20").

nm <- c("IgG3 AMA-1 MFI | V20", "IgG1 CSP-C (delta)", "Age at V5", "x-y|z w")
d <- gemcox_simulate(n = 200, p = 4, mu_sep = 1, beta_sep = 1, seed = 81)
X <- d$X
colnames(X) <- nm

test_that("names with spaces, hyphens and '|' survive a matrix fit", {
  fit <- suppressWarnings(gemcox(X, time = d$time, status = d$status))
  expect_identical(rownames(coef(fit)), nm)
  expect_identical(rownames(fit$mu), nm)
  expect_identical(colnames(fit$Sigma[[1]]), nm)
  expect_identical(names(fit$gmm_scaler$center), nm)
  expect_identical(colnames(fit$data$X_cox), nm)
  expect_identical(rownames(summary(fit)$coefficients), nm)
})

test_that("names survive a data frame built with check.names = FALSE", {
  df <- data.frame(X, check.names = FALSE)
  expect_identical(names(df), nm)
  fit <- suppressWarnings(gemcox(df, time = d$time, status = d$status))
  expect_identical(rownames(coef(fit)), nm)
  expect_equal(predict(fit, df[, rev(nm)]), predict(fit, X))
})

test_that("separate X_gmm and X_cox keep their own names", {
  Xc <- cbind(X, `CSP | V5` = stats::rnorm(200))
  fit <- suppressWarnings(gemcox(X[, 1:2], Xc, time = d$time, status = d$status))
  expect_identical(rownames(fit$mu), nm[1:2])
  expect_identical(rownames(coef(fit)), c(nm, "CSP | V5"))
})
