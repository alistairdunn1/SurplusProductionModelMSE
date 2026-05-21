test_that("MSE package loads", {
  expect_true("MSE" %in% loadedNamespaces() ||
    requireNamespace("MSE", quietly = TRUE))
})
