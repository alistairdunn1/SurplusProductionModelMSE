# Tests for reusable scenario-builder helpers

test_that("create_baseline_scenarios returns expected scenario names", {
  scenarios <- create_baseline_scenarios(
    f_target = 0.03,
    assessment_frequency = 2L,
    catch_allocation = c(A = 0.6, B = 0.4)
  )

  expect_length(scenarios, 4)
  expect_true(all(vapply(scenarios, inherits, logical(1), "mse_scenario")))
  expect_equal(
    vapply(scenarios, function(x) x$name, character(1)),
    c("no_catch", "constant_f_0.03", "hockey_20_50", "hockey_10_40")
  )
})

test_that("create_hcr_grid_scenarios returns expected count and names", {
  f_grid <- c(0.01, 0.02, 0.03)
  scenarios <- create_hcr_grid_scenarios(
    f_grid = f_grid,
    assessment_frequency = 1L,
    catch_allocation = c(A = 0.5, B = 0.5)
  )

  expect_length(scenarios, length(f_grid) * 3)
  expect_equal(
    vapply(scenarios, function(x) x$name, character(1)),
    c(
      "cf_0.010", "hs2050_0.010", "hs1040_0.010",
      "cf_0.020", "hs2050_0.020", "hs1040_0.020",
      "cf_0.030", "hs2050_0.030", "hs1040_0.030"
    )
  )
})

test_that("create_hcr_grid_scenarios include selection changes count", {
  scenarios <- create_hcr_grid_scenarios(
    f_grid = c(0.02, 0.03),
    assessment_frequency = 1L,
    include = c("constant_f", "hockey_10_40")
  )

  expect_length(scenarios, 4)
  expect_equal(
    vapply(scenarios, function(x) x$name, character(1)),
    c("cf_0.020", "hs1040_0.020", "cf_0.030", "hs1040_0.030")
  )
})

test_that("scenario builders validate assessment_frequency", {
  expect_error(
    create_baseline_scenarios(0.03, assessment_frequency = 0),
    "assessment_frequency"
  )
  expect_error(
    create_hcr_grid_scenarios(c(0.02, 0.03), assessment_frequency = 0),
    "assessment_frequency"
  )
})

test_that("scenario builders validate catch_allocation", {
  expect_error(
    create_baseline_scenarios(0.03, catch_allocation = c(A = -1, B = 2)),
    "catch_allocation"
  )
  expect_error(
    create_hcr_grid_scenarios(c(0.02, 0.03), catch_allocation = c(A = 0, B = 0)),
    "positive"
  )
})

test_that("create_hcr_grid_scenarios validates include values", {
  expect_error(
    create_hcr_grid_scenarios(c(0.02, 0.03), include = c("constant_f", "bad_rule")),
    "Unknown include option"
  )
})
