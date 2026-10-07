test_that("package namespace loads with its declared imports", {
  expect_true(isNamespaceLoaded("gemcox"))
  for (pkg in c("survival", "glmnet", "stats")) {
    expect_true(requireNamespace(pkg, quietly = TRUE), info = pkg)
  }
})
