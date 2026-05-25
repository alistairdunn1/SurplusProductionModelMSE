test_that("SurplusProductionModelMSE package loads", {
  expect_true("SurplusProductionModelMSE" %in% loadedNamespaces() ||
    requireNamespace("SurplusProductionModelMSE", quietly = TRUE))
})
