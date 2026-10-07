## The GeMLR-style workflow: read_data -> fit_model -> runCV -> finalModel ->
## plot_beta_heatmap. The wrappers must reproduce the underlying functions.

test_that("read_data splits biomarkers, outcome and indicators from a data frame or a file", {
  x <- gemcox_example
  a <- read_data(x, time_col = "time", status_col = "status", Indi_col = "vaccine")
  expect_equal(a$dim, 6); expect_equal(a$numdata, 300)
  expect_identical(colnames(a$X), paste0("marker", 1:6))
  expect_identical(colnames(a$Indi), "vaccine")
  expect_equal(a$time, x$time); expect_equal(a$status, x$status)
  expect_equal(unname(apply(a$Xs, 2, sd)), rep(1, 6), tolerance = 1e-12)
  f <- tempfile(fileext = ".csv"); on.exit(unlink(f))
  utils::write.csv(x, f, row.names = FALSE)
  b <- read_data(f, time_col = "time", status_col = "status", Indi_col = "vaccine")
  expect_equal(b$X, a$X, tolerance = 1e-12)
  r <- tempfile(fileext = ".rds"); saveRDS(x, r)
  expect_equal(read_data(r, time_col = 8, status_col = 9)$dim, 7)   # vaccine becomes a biomarker
  expect_error(read_data(x, time_col = "time"), "status_col")
  y <- x; y$marker1[3] <- NA
  expect_error(read_data(y, time_col = "time", status_col = "status"), "Missing values in: marker1")
  y <- x; y$site <- "A"
  expect_error(read_data(y, time_col = "time", status_col = "status"), "non-numeric: site")
})

test_that("fit_model is gemcox() with indicators in the Cox models only", {
  a <- read_data(gemcox_example, time_col = "time", status_col = "status", Indi_col = "vaccine")
  f <- suppressWarnings(fit_model(a$X, a$time, a$status, Indi = a$Indi, K = 2, nseeds = 2))
  g <- suppressWarnings(gemcox(as.matrix(a$X), cbind(as.matrix(a$X), a$Indi), time = a$time,
                               status = a$status, K = 2, n_starts = 2, seed = 1))
  expect_identical(f$beta, g$beta)
  expect_identical(f$beta_sd, g$beta_scaled)
  expect_identical(rownames(f$beta_sd), c(paste0("marker", 1:6), "vaccine"))
  expect_identical(f$Indi_names, "vaccine")
  expect_true(f$metrics$cindex > 0.5 && f$metrics$cindex < 1)
  h <- suppressWarnings(fit_model(a$X, a$time, a$status, vargmm = 1:3, varcox = c("marker1", "marker2"),
                                  K = 2, nseeds = 1))
  expect_identical(h$vargmm, paste0("marker", 1:3))
  expect_identical(rownames(h$beta), c("marker1", "marker2"))
  expect_error(fit_model(a$X, a$time, a$status, vargmm = "nope"), "vargmm")
})

test_that("runCV gives per-fold values that sum to gemcox_cv_loglik, and finalModel picks the best K", {
  X <- as.matrix(gemcox_example[, paste0("marker", 1:6)])
  tt <- gemcox_example$time; st <- gemcox_example$status
  cv <- suppressWarnings(runCV(X, tt, st, ncmp = 1:2, k = 3, verbose = 0))
  expect_identical(dim(cv$cvLLfinal), c(3L, 2L))
  expect_identical(colnames(cv$cvLLfinal), c("cluster=1", "cluster=2"))
  ref <- suppressWarnings(gemcox_cv_loglik(X, time = tt, status = st, K = 1:2, nfolds = 3, seed = 1,
                                           lambda = 0.05, n_starts = 1))
  expect_equal(unname(colSums(cv$cvLLfinal)), unname(as.numeric(ref)), tolerance = 1e-8)
  expect_identical(cv$best_K, (1:2)[which.max(colMeans(cv$cvLLfinal))])
  fin <- suppressWarnings(finalModel(cv, X, tt, st, nseeds = 1))
  expect_identical(fin$K_selected, cv$best_K)
  expect_identical(ncol(fin$beta), as.integer(cv$best_K))
})

test_that("plot_beta_heatmap draws a fit or a matrix and can save to a file", {
  a <- read_data(gemcox_example, time_col = "time", status_col = "status")
  f <- suppressWarnings(fit_model(a$X, a$time, a$status, K = 2, nseeds = 1))
  png_file <- tempfile(fileext = ".png"); on.exit(unlink(png_file))
  m <- plot_beta_heatmap(f, output_file = png_file)
  expect_true(file.exists(png_file) && file.size(png_file) > 0)
  expect_identical(m, f$beta_sd)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_identical(plot_beta_heatmap(matrix(1:4, 2)), matrix(1:4, 2))
})
