## The single data-generating mechanism: mu_sep and beta_sep are controlled
## independently and the clusters have equal coefficient norms.

test_that("clusters have equal coefficient norms and the design directions are orthogonal", {
  for (p in c(4, 10, 80)) {
    d <- gemcox_simulate(n = 50, p = p, mu_sep = 1, beta_sep = 2, seed = 121)
    b1 <- d$beta_list[[1]]
    b2 <- d$beta_list[[2]]
    expect_equal(sqrt(sum(b1^2)), sqrt(sum(b2^2)), tolerance = 1e-12)
    v <- (b2 - b1) / sqrt(sum((b2 - b1)^2))
    base <- (b1 + b2) / 2
    expect_equal(sum(base * v), 0, tolerance = 1e-12)
    ## the profile direction u is orthogonal to v
    mu_diff <- colMeans(d$X[d$Z == 2, , drop = FALSE]) - colMeans(d$X[d$Z == 1, , drop = FALSE])
    expect_equal(length(mu_diff), p)
  }
})

test_that("u and v are orthonormal (checked on noise-free means)", {
  d <- gemcox_simulate(n = 20000, p = 6, mu_sep = 2, beta_sep = 2, seed = 122)
  mu_diff <- colMeans(d$X[d$Z == 2, ]) - colMeans(d$X[d$Z == 1, ])
  v <- d$beta_list[[2]] - d$beta_list[[1]]
  expect_equal(sqrt(sum(mu_diff^2)), 2, tolerance = 0.05)
  expect_lt(abs(sum(mu_diff * v)) / (sqrt(sum(mu_diff^2)) * sqrt(sum(v^2))), 0.03)
})

test_that("event rate is calibrated and seeds are reproducible without touching the RNG", {
  set.seed(8)
  before <- .Random.seed
  a <- gemcox_simulate(n = 2000, p = 5, mu_sep = 0, beta_sep = 1, seed = 123)
  expect_identical(.Random.seed, before)
  expect_identical(a, gemcox_simulate(n = 2000, p = 5, mu_sep = 0, beta_sep = 1, seed = 123))
  expect_lt(abs(mean(a$status) - 0.55), 0.03)
  expect_true(all(a$time >= 0.1))
})

test_that("K_true = 1 gives one cluster with the base coefficient", {
  d <- gemcox_simulate(n = 100, p = 5, mu_sep = 3, beta_sep = 3, K_true = 1, seed = 124)
  expect_true(all(d$Z == 1))
  expect_length(d$beta_list, 1)
  expect_gt(sqrt(sum(d$beta_list[[1]]^2)), 0.5)   # a nonzero shared effect
})

test_that("K_true = 3 uses an equilateral triangle geometry with equal coefficient norms", {
  d <- gemcox_simulate(n = 300, p = 10, mu_sep = 1, beta_sep = 2, K_true = 3, seed = 125)
  expect_setequal(unique(d$Z), 1:3)
  norms <- vapply(d$beta_list, function(b) sqrt(sum(b^2)), 0)
  expect_equal(max(norms) - min(norms), 0, tolerance = 1e-12)
  expect_equal(as.numeric(dist(do.call(rbind, d$beta_list))), rep(2, 3), tolerance = 1e-12)
  expect_error(gemcox_simulate(50, 6, 1, 1, K_true = 3, seed = 1), "p >= 8")
  expect_error(gemcox_simulate(50, 10, 1, 1, K_true = 4, seed = 1), "K_true must be")
})

test_that("K_true = 1 and 2 draws are unchanged by the K_true = 3 extension", {
  ## Checksums computed with the frozen sim-freeze-v1 version of
  ## gemcox_simulate() under the Mersenne-Twister generator.
  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  ck <- function(d) c(sum(d$X), sum(d$time), sum(d$status), sum(d$Z))
  frozen <- list(
    list(K_true = 2, rho = 0,   ck = c(0.5115530178, 7974.1771050688, 30, 88)),
    list(K_true = 1, rho = 0,   ck = c(3.5725587508, 7679.6624758164, 34, 60)),
    list(K_true = 2, rho = 0.5, ck = c(26.7871996394, 8521.9395764449, 34, 88)))
  for (f in frozen) {
    d <- gemcox_simulate(n = 60, p = 10, mu_sep = 1, beta_sep = 2, rho = f$rho,
                         K_true = f$K_true, seed = 126)
    expect_equal(ck(d), f$ck, tolerance = 1e-9)
  }
})

test_that("features = 'lognormal' standardises skewed noise and keeps the design", {
  ## Added after sim-freeze-v2 for experiment C5. The Gaussian default is
  ## covered by the checksum tests above and in test-lcc.R.
  g <- gemcox_simulate(n = 20000, p = 6, mu_sep = 1, beta_sep = 2, seed = 127)
  l <- gemcox_simulate(n = 20000, p = 6, mu_sep = 1, beta_sep = 2, features = "lognormal",
                       seed = 127)
  ## same draws: same labels and coefficients, and within a cluster each
  ## feature is an increasing function of its Gaussian counterpart
  expect_identical(l$Z, g$Z)
  expect_identical(l$beta_list, g$beta_list)
  for (j in 1:6) {
    expect_equal(cor(l$X[l$Z == 1, j], g$X[g$Z == 1, j], method = "spearman"), 1)
  }
  ## standardised and right-skewed about the cluster means
  z <- l$X[l$Z == 1, ] - matrix(colMeans(l$X[l$Z == 1, ]), sum(l$Z == 1), 6, byrow = TRUE)
  expect_equal(unname(apply(z, 2, sd)), rep(1, 6), tolerance = 0.06)
  skew <- apply(z, 2, function(x) mean(x^3) / sd(x)^3)
  expect_true(all(skew > 2))                     # theoretical 3.69 at sdlog = 0.8
  ## same separation: the class-mean difference matches the Gaussian design
  dl <- colMeans(l$X[l$Z == 2, ]) - colMeans(l$X[l$Z == 1, ])
  dg <- colMeans(g$X[g$Z == 2, ]) - colMeans(g$X[g$Z == 1, ])
  expect_lt(max(abs(dl - dg)), 0.06)
  expect_lt(abs(mean(l$status) - 0.55), 0.03)
  expect_error(gemcox_simulate(50, 6, 1, 1, rho = 0.5, features = "lognormal", seed = 1),
               "requires rho = 0")
})

test_that("administrative censoring: fixed follow-up, no other censoring (E6)", {
  d <- gemcox_simulate(n = 2000, p = 6, mu_sep = 1, beta_sep = 2, censoring = "administrative",
                       followup = 150, seed = 128)
  expect_true(all(d$time <= 150))
  expect_true(all(d$status[d$time < 150] == 1))      # censored only at the follow-up time
  expect_true(all(d$status[d$time == 150] == 0))
  r <- gemcox_simulate(n = 2000, p = 6, mu_sep = 1, beta_sep = 2, seed = 128)
  expect_identical(d$X, r$X)                          # same features and labels as random censoring
  expect_identical(d$Z, r$Z)
  expect_equal(length(d$mu_list), 2)
  u <- d$mu_list[[2]] - d$mu_list[[1]]
  expect_equal(sqrt(sum(u^2)), 1, tolerance = 1e-12)  # |mu_2 - mu_1| = mu_sep
  expect_error(gemcox_simulate(50, 6, 1, 1, censoring = "administrative", seed = 1), "followup")
})
