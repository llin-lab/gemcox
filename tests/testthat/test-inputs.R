d <- gemcox_simulate(n = 100, p = 3, mu_sep = 1, beta_sep = 1, seed = 91)

test_that("mismatched dimensions give informative errors", {
  expect_error(gemcox(d$X, d$X[-1, ], time = d$time, status = d$status),
               "X_gmm has 100 rows but X_cox has 99")
  expect_error(gemcox(d$X, time = d$time[-1], status = d$status),
               "time must be a numeric vector of length 100")
  expect_error(gemcox(d$X, time = d$time, status = d$status[-1]),
               "status must be a vector of length 100")
})

test_that("non-binary status is rejected", {
  expect_error(gemcox(d$X, time = d$time, status = d$status + 1), "status must be binary")
  s <- d$status
  s[3] <- NA
  expect_error(gemcox(d$X, time = d$time, status = s), "status must be binary")
  expect_error(gemcox(d$X, time = d$time, status = rep(0, 100)), "no events")
})

test_that("negative or missing times are rejected", {
  tt <- d$time
  tt[5] <- -1
  expect_error(gemcox(d$X, time = tt, status = d$status), "time must be non-negative")
  tt[5] <- NA
  expect_error(gemcox(d$X, time = tt, status = d$status), "time contains NA")
})

test_that("NA and non-finite features are rejected, naming the column", {
  X <- d$X
  X[2, "f2"] <- NA
  expect_error(gemcox(X, time = d$time, status = d$status), "X_gmm contains NA.*f2")
  expect_error(gemcox(d$X, X, time = d$time, status = d$status), "X_cox contains NA")
  X[2, "f2"] <- Inf
  expect_error(gemcox(X, time = d$time, status = d$status), "non-finite")
})

test_that("non-numeric columns and duplicated names are rejected", {
  df <- data.frame(d$X, grp = "a")
  expect_error(gemcox(df, time = d$time, status = d$status), "non-numeric column.*grp")
  X <- d$X
  colnames(X) <- c("a", "a", "b")
  expect_error(gemcox(X, time = d$time, status = d$status), "duplicated column names")
})

test_that("invalid settings are rejected", {
  expect_error(gemcox(d$X, time = d$time, status = d$status, K = 0), "K must be")
  expect_error(gemcox(d$X, time = d$time, status = d$status, K = 1.5), "K must be")
  expect_error(gemcox(d$X, time = d$time, status = d$status, gamma = -1), "gamma must be")
  expect_error(gemcox(d$X, time = d$time, status = d$status, alpha = 2), "alpha must be")
  expect_error(gemcox(d$X, time = d$time, status = d$status, lambda = NULL), "lambda must be")
  expect_error(gemcox(d$X, time = d$time, status = d$status, init = "supervised"),
               "should be one of")
  expect_error(gemcox_control(max_jump = 0), "max_jump must be")
})

test_that("logical status is accepted", {
  f <- suppressWarnings(gemcox(d$X, time = d$time, status = d$status == 1, K = 1))
  expect_s3_class(f, "gemcox")
})
