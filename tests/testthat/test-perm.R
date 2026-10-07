## Inference: the cross-validated partial likelihood, the K = 2 vs K = 1
## test, and its two null distributions.

test_that("cox_partial_loglik equals coxph's Breslow log-likelihood, with ties", {
  d <- gemcox_simulate(n = 200, p = 3, mu_sep = 0, beta_sep = 1, seed = 101)
  tt <- round(d$time)   # many ties
  eta <- as.numeric(d$X %*% c(0.3, -0.2, 0.5))
  ref <- coxph(Surv(tt, d$status) ~ offset(eta), ties = "breslow")$loglik
  expect_equal(cox_partial_loglik(tt, d$status, eta), ref, tolerance = 1e-10)
})

test_that("cross-validated partial likelihood follows Verweij & van Houwelingen", {
  d <- gemcox_simulate(n = 200, p = 3, mu_sep = 1, beta_sep = 1, seed = 102)
  folds <- stratified_folds(d$status, 4, seed = 1)
  cv <- suppressWarnings(cv_loglik_folds(d$X, d$X, d$time, d$status, folds, K = 1,
                                         criterion = "partial", args = list()))
  manual <- 0
  for (f in 1:4) {
    tr <- folds != f
    fit <- suppressWarnings(gemcox(d$X[tr, ], time = d$time[tr], status = d$status[tr], K = 1))
    eta <- as.numeric(d$X %*% coef(fit))   # the centring constant cancels at K = 1
    manual <- manual +
      coxph(Surv(d$time, d$status) ~ offset(eta), ties = "breslow")$loglik -
      coxph(Surv(d$time[tr], d$status[tr]) ~ offset(eta[tr]), ties = "breslow")$loglik
  }
  expect_equal(as.numeric(cv), manual, tolerance = 1e-8)
})

test_that("stratified folds balance events and leave the caller's RNG alone", {
  set.seed(5)
  before <- .Random.seed
  st <- rep(0:1, c(60, 40))
  f <- stratified_folds(st, 4, seed = 9)
  expect_identical(.Random.seed, before)
  expect_identical(f, stratified_folds(st, 4, seed = 9))
  expect_true(all(table(f[st == 1]) == 10))
})

test_that("the test returns a valid htest with the stated p-value formula", {
  d <- gemcox_simulate(n = 150, p = 3, mu_sep = 0, beta_sep = 2, seed = 103)
  set.seed(1)
  before <- .Random.seed
  for (nl in c("permutation", "bootstrap")) {
    pt <- suppressWarnings(gemcox_heterogeneity_test(
      X_gmm = d$X, time = d$time, status = d$status, B = 5, nfolds = 3, statistic = "cv",
      null = nl, seed = 7))
    expect_s3_class(pt, "htest")
    expect_identical(pt$null, nl)
    ok <- is.finite(pt$null.distribution)
    expect_equal(pt$n_valid, sum(ok))
    expect_equal(pt$p.value,
                 (1 + sum(pt$null.distribution[ok] >= pt$statistic)) / (1 + sum(ok)))
    expect_equal(unname(pt$statistic), unname(pt$cv_loglik["K2"] - pt$cv_loglik["K1"]))
    pt2 <- suppressWarnings(gemcox_heterogeneity_test(
      X_gmm = d$X, time = d$time, status = d$status, B = 5, nfolds = 3, statistic = "cv",
      null = nl, seed = 7))
    expect_identical(pt2$null.distribution, pt$null.distribution)
  }
  expect_identical(.Random.seed, before)
})

test_that("a fitted object can be passed directly, reusing its settings", {
  d <- gemcox_simulate(n = 150, p = 3, mu_sep = 0, beta_sep = 2, seed = 104)
  fit <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, lambda = 0.1))
  a <- suppressWarnings(gemcox_heterogeneity_test(fit, statistic = "cv", B = 3, nfolds = 3,
                                                  seed = 2))
  b <- suppressWarnings(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, statistic = "cv",
                                                status = d$status, lambda = 0.1,
                                                B = 3, nfolds = 3, seed = 2))
  expect_equal(a$statistic, b$statistic)
  expect_equal(a$null.distribution, b$null.distribution)
})

test_that("the parametric bootstrap sampler reproduces the fitted K = 1 model", {
  d <- gemcox_simulate(n = 1500, p = 3, mu_sep = 0, beta_sep = 0, K_true = 1, seed = 105)
  f1 <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 1, lambda = 0))
  draws <- with_local_seed(3, {
    s <- k1_bootstrap_sampler(f1, d$X, d$time, d$status)
    replicate(15, s(), simplify = FALSE)
  })
  rate <- mean(vapply(draws, function(y) mean(y$status), 0))
  expect_lt(abs(rate - mean(d$status)), 0.03)
  b <- rowMeans(vapply(draws, function(y) {
    as.numeric(coef(coxph(Surv(y$time, y$status) ~ d$X, ties = "breslow")))
  }, numeric(3)))
  expect_lt(max(abs(b - coef(f1))), 0.05)
  expect_true(all(vapply(draws, function(y) max(y$time) <= max(d$time), NA)))
})

## Calibration SMOKE test only: 50 small K = 1 datasets with a nonzero
## shared coefficient vector. The real calibration (>= 500 datasets, both
## nulls, type I error and power side by side) is experiment E5a.
## With 19 null replicates the exact rejection rate at 0.05 is 1/20, so at
## most 7 rejections of 50 are accepted (P(Bin(50, 0.05) >= 8) ~ 0.012).
test_that("both nulls are roughly calibrated on a nonzero-beta K = 1 null (smoke)", {
  skip_on_cran()
  for (nl in c("permutation", "bootstrap")) {
    p <- vapply(1:50, function(i) {
      d <- gemcox_simulate(n = 200, p = 4, mu_sep = 0, beta_sep = 0, K_true = 1,
                           seed = 5000 + i)
      suppressWarnings(gemcox_heterogeneity_test(
        X_gmm = d$X, time = d$time, status = d$status, B = 19, nfolds = 3, statistic = "cv",
        null = nl, seed = i))$p.value
    }, 0)
    expect_lte(sum(p <= 0.05), 7, label = sprintf("%s null: rejections of 50", nl))
  }
})

test_that("bootstrap is the default null, and `observed` reuses the observed statistic", {
  d <- gemcox_simulate(n = 150, p = 3, mu_sep = 0, beta_sep = 2, seed = 106)
  run <- function(...) {
    suppressWarnings(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status,
                                               B = 4, seed = 7, statistic = "cv", ...))
  }
  bs <- run(nfolds = 3)
  expect_identical(bs$null, "bootstrap")
  pm_full <- run(nfolds = 3, null = "permutation")
  pm_reuse <- run(nfolds = 3, null = "permutation", observed = bs)
  expect_identical(pm_reuse$statistic, bs$statistic)
  expect_identical(pm_reuse$cv_loglik, pm_full$cv_loglik)
  expect_identical(pm_reuse$null.distribution, pm_full$null.distribution)
  expect_identical(pm_reuse$p.value, pm_full$p.value)
  expect_error(run(nfolds = 4, null = "permutation", observed = bs), "same data")
})

test_that("the CV criterion exposes partial and GMM components that add up", {
  d <- gemcox_simulate(n = 200, p = 3, mu_sep = 2, beta_sep = 1, seed = 107)
  j <- suppressWarnings(gemcox_cv_loglik(d$X, time = d$time, status = d$status,
                                         nfolds = 3, seed = 2, criterion = "joint"))
  p <- suppressWarnings(gemcox_cv_loglik(d$X, time = d$time, status = d$status,
                                         nfolds = 3, seed = 2, criterion = "partial"))
  cj <- attr(j, "components")
  expect_equal(as.numeric(j), as.numeric(colSums(cj)), tolerance = 1e-10)
  expect_equal(as.numeric(p), as.numeric(cj["partial", ]), tolerance = 1e-12)
  expect_identical(attr(p, "components"), cj)
})

test_that("the in-sample LRT equals 2 (l2 - l1) from direct fits and is calibrated by the nulls", {
  d <- gemcox_simulate(n = 150, p = 3, mu_sep = 0, beta_sep = 2, seed = 108)
  run <- function(...) {
    suppressWarnings(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status,
                                               B = 4, seed = 9, statistic = "lrt", ...))
  }
  bs <- run()
  f1 <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 1))
  f2 <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 2))
  expect_equal(unname(bs$statistic), 2 * (f2$final_score - f1$final_score), tolerance = 1e-10)
  expect_identical(bs$statistic_type, "lrt")
  expect_equal(bs$p.value,
               (1 + sum(bs$null.distribution >= bs$statistic, na.rm = TRUE)) / (1 + bs$n_valid))
  pm_full <- run(null = "permutation")
  pm_reuse <- run(null = "permutation", observed = bs)
  expect_identical(pm_reuse$null.distribution, pm_full$null.distribution)
  expect_identical(pm_reuse$statistic, pm_full$statistic)
  expect_error(run(gamma = 0), "gamma = 1 and temp = 1")
  cv <- suppressWarnings(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status,
                                                   statistic = "cv", B = 2, nfolds = 3, seed = 9))
  expect_error(run(observed = cv), "same statistic")
})

test_that("the default test is the LRT with the bootstrap null", {
  d <- gemcox_simulate(n = 150, p = 3, mu_sep = 0, beta_sep = 2, seed = 109)
  ht <- suppressWarnings(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time,
                                                   status = d$status, B = 3, seed = 4))
  expect_identical(ht$statistic_type, "lrt")
  expect_identical(ht$null, "bootstrap")
  expect_match(ht$method, "parametric bootstrap null, in-sample likelihood ratio")
  expect_error(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status,
                                         B = 0), "number of null replicates")
})

test_that("gemcox_permutation_test() is a deprecated alias with its old defaults", {
  d <- gemcox_simulate(n = 150, p = 3, mu_sep = 0, beta_sep = 2, seed = 110)
  expect_warning(old <- gemcox_permutation_test(X_gmm = d$X, time = d$time, status = d$status,
                                                n_perm = 3, nfolds = 3, seed = 5),
                 class = "deprecatedWarning")
  new <- suppressWarnings(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time,
                                                    status = d$status, statistic = "cv",
                                                    B = 3, nfolds = 3, seed = 5))
  expect_identical(old$statistic_type, "cv")
  expect_identical(old$null.distribution, new$null.distribution)
  expect_identical(old$p.value, new$p.value)
})

test_that("K0 = 2 tests three clusters against two with a K = 2 bootstrap null", {
  d <- gemcox_simulate(n = 200, p = 8, mu_sep = 3, beta_sep = 2, K_true = 3, seed = 111)
  ht <- suppressWarnings(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time,
                                                   status = d$status, K0 = 2, B = 3, seed = 6))
  f2 <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 2))
  f3 <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 3))
  expect_equal(unname(ht$statistic), 2 * (f3$final_score - f2$final_score), tolerance = 1e-10)
  expect_identical(ht$null_K, 2)
  expect_match(ht$method, "K = 3 vs K = 2")
  expect_match(names(ht$statistic), "l3 - l2")
  expect_error(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status,
                                         K0 = 2, statistic = "cv"), "K0 >= 2 needs")
  expect_error(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status,
                                         K0 = 2, null = "permutation"), "K0 >= 2 needs")
})

test_that("the K0 >= 2 bootstrap sampler reproduces the fitted model's event rate", {
  d <- gemcox_simulate(n = 1500, p = 8, mu_sep = 3, beta_sep = 1, K_true = 2, seed = 112)
  f2 <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status, K = 2, lambda = 0))
  draws <- with_local_seed(4, {
    s <- kK_bootstrap_sampler(f2, d$X, d$X, d$time, d$status)
    replicate(10, s(), simplify = FALSE)
  })
  rate <- mean(vapply(draws, function(y) mean(y$status), 0))
  expect_lt(abs(rate - mean(d$status)), 0.03)
  expect_true(all(vapply(draws, function(y) max(y$time) <= max(d$time), NA)))
})

test_that("the LRT counts fits that stop at max_iter (added after the frozen study)", {
  d <- gemcox_simulate(n = 150, p = 3, mu_sep = 0, beta_sep = 2, seed = 7)
  args <- list(X_gmm = d$X, time = d$time, status = d$status, B = 3, seed = 2)
  ## the count adds nothing to the test itself
  a <- suppressWarnings(do.call(gemcox_heterogeneity_test, args))
  b <- suppressWarnings(do.call(gemcox_heterogeneity_test, c(args, list(max_iter = 1))))
  ## null-model fit + 2 fits on the observed data + 2 per null replicate
  expect_equal(a$n_fits, 1 + 2 + 2 * 3)
  expect_equal(b$n_fits, a$n_fits)
  ## with max_iter = 1 every K = 2 fit stops at the cap
  expect_gte(b$n_not_converged, 1 + 3)
  expect_lte(a$n_not_converged, b$n_not_converged)
  cv <- suppressWarnings(gemcox_heterogeneity_test(X_gmm = d$X, time = d$time, status = d$status,
                                                   statistic = "cv", B = 2, nfolds = 3, seed = 2))
  expect_true(is.na(cv$n_not_converged))
})
