# Tests for MSE configuration S3 classes
# Tests constructors, validation, and print methods for:
#   om_config, em_config, impl_error, mse_scenario

# ── om_config ──────────────────────────────────────────────────────────────────

test_that("om_config creates valid single-area object", {
  om <- om_config(
    true_params = list(
      r = 0.3, K = 5000, m = 2,
      sigma_obs = 0.2, q = 1e-4, B0 = 4000
    )
  )
  expect_s3_class(om, "om_config")
  expect_equal(om$n_areas, 1L)
  expect_equal(om$movement_rate, 0)
  expect_null(om$movement_cost_matrix)
  expect_null(om$attractiveness)
  expect_equal(om$decay, 0)
  expect_equal(om$true_params$r, 0.3)
})

test_that("om_config creates valid multi-area object", {
  dm <- matrix(c(0, 100, 200, 100, 0, 100, 200, 100, 0), 3, 3)
  om <- om_config(
    n_areas = 3,
    movement_rate = 0.1,
    movement_cost_matrix = dm,
    attractiveness = c(1, 1.2, 0.8),
    decay = 0.01,
    true_params = list(
      r = 0.3, K = 5000, m = 2,
      sigma_obs = 0.2, q = rep(1e-4, 3),
      B0 = c(2000, 2500, 1500)
    )
  )
  expect_s3_class(om, "om_config")
  expect_equal(om$n_areas, 3L)
  expect_equal(om$movement_rate, 0.1)
  expect_equal(nrow(om$movement_cost_matrix), 3)
  expect_length(om$attractiveness, 3)
})

test_that("om_config allows NULL true_params", {
  om <- om_config(true_params = NULL)
  expect_s3_class(om, "om_config")
  expect_null(om$true_params)
})

test_that("om_config rejects invalid n_areas", {
  expect_error(om_config(n_areas = 0))
  expect_error(om_config(n_areas = -1))
  expect_error(om_config(n_areas = 1.5))
  expect_error(om_config(n_areas = "abc"))
})

test_that("om_config rejects invalid movement_rate", {
  expect_error(om_config(movement_rate = -0.1))
  expect_error(om_config(movement_rate = 1.5))
})

test_that("om_config allows asymmetric movement_cost_matrix", {
  bad_dm <- matrix(c(0, 100, 200, 50, 0, 100, 200, 100, 0), 3, 3)
  expect_s3_class(
    om_config(
      n_areas = 3, movement_cost_matrix = bad_dm,
      true_params = list(
        r = 0.3, K = 5000, m = 2,
        sigma_obs = 0.2, q = rep(1e-4, 3),
        B0 = rep(2000, 3)
      )
    ),
    "om_config"
  )
})

test_that("om_config rejects wrong-sized movement_cost_matrix", {
  dm2 <- matrix(c(0, 100, 100, 0), 2, 2)
  expect_error(
    om_config(
      n_areas = 3, movement_cost_matrix = dm2,
      true_params = list(
        r = 0.3, K = 5000, m = 2,
        sigma_obs = 0.2, q = rep(1e-4, 3),
        B0 = rep(2000, 3)
      )
    )
  )
})

test_that("om_config rejects movement_cost_matrix with non-zero diagonal", {
  bad_dm <- matrix(c(1, 100, 100, 0), 2, 2)
  expect_error(
    om_config(
      n_areas = 2, movement_cost_matrix = bad_dm,
      true_params = list(
        r = 0.3, K = 5000, m = 2,
        sigma_obs = 0.2, q = rep(1e-4, 2),
        B0 = rep(2000, 2)
      )
    )
  )
})

test_that("om_config rejects missing required true_params fields", {
  expect_error(
    om_config(true_params = list(r = 0.3, K = 5000))
  )
})

test_that("om_config rejects invalid parameter values", {
  base_params <- list(
    r = 0.3, K = 5000, m = 2,
    sigma_obs = 0.2, q = 1e-4, B0 = 4000
  )
  # r out of range
  expect_error(om_config(true_params = modifyList(base_params, list(r = 0))))
  expect_error(om_config(true_params = modifyList(base_params, list(r = 3))))
  # K must be positive

  expect_error(om_config(true_params = modifyList(base_params, list(K = 0))))
  # sigma_obs must be positive
  expect_error(om_config(true_params = modifyList(base_params, list(sigma_obs = 0))))
})

test_that("om_config print method works", {
  om <- om_config(
    true_params = list(
      r = 0.3, K = 5000, m = 2,
      sigma_obs = 0.2, q = 1e-4, B0 = 4000
    )
  )
  out <- capture.output(print(om))
  expect_true(any(grepl("Operating Model", out)))
  expect_true(any(grepl("Areas", out)))
  expect_true(any(grepl("r =", out)))
})

test_that("validate_om_config rejects wrong class", {
  expect_error(validate_om_config(list()), "not of class")
})


# ── em_config ──────────────────────────────────────────────────────────────────

test_that("em_config creates valid default object", {
  em <- em_config()
  expect_s3_class(em, "em_config")
  expect_equal(em$n_areas, 1L)
  expect_false(em$estimate_movement)
  expect_false(em$aggregate_areas)
  expect_null(em$fixed_params)
})

test_that("em_config creates valid configured object", {
  em <- em_config(
    n_areas = 1,
    aggregate_areas = TRUE,
    fixed_params = list(m = 2)
  )
  expect_s3_class(em, "em_config")
  expect_true(em$aggregate_areas)
  expect_equal(em$fixed_params$m, 2)
})

test_that("em_config rejects invalid n_areas", {
  expect_error(em_config(n_areas = 0))
  expect_error(em_config(n_areas = -1))
})

test_that("em_config rejects non-logical flags", {
  expect_error(em_config(estimate_movement = "yes"))
  expect_error(em_config(aggregate_areas = 1))
})

test_that("em_config rejects unnamed fixed_params", {
  expect_error(em_config(fixed_params = list(2, 0.3)))
})

test_that("em_config print method works", {
  em <- em_config(fixed_params = list(m = 2))
  out <- capture.output(print(em))
  expect_true(any(grepl("Estimation Model", out)))
  expect_true(any(grepl("Fixed parameters", out)))
  expect_true(any(grepl("m = 2", out)))
})

test_that("validate_em_config rejects wrong class", {
  expect_error(validate_em_config(list()), "not of class")
})


# ── impl_error ─────────────────────────────────────────────────────────────────

test_that("impl_error creates valid default object", {
  ie <- impl_error()
  expect_s3_class(ie, "impl_error")
  expect_equal(ie$cv, 0.1)
  expect_equal(ie$bias, 1)
  expect_equal(ie$autocorr, 0)
  expect_equal(ie$max_overage, 1.1)
})

test_that("impl_error creates valid custom object", {
  ie <- impl_error(cv = 0.2, bias = 1.05, autocorr = 0.5, max_overage = 1.2)
  expect_s3_class(ie, "impl_error")
  expect_equal(ie$cv, 0.2)
  expect_equal(ie$bias, 1.05)
  expect_equal(ie$autocorr, 0.5)
  expect_equal(ie$max_overage, 1.2)
})

test_that("impl_error rejects invalid cv", {
  expect_error(impl_error(cv = 0))
  expect_error(impl_error(cv = -0.1))
})

test_that("impl_error rejects invalid bias", {
  expect_error(impl_error(bias = -1))
})

test_that("impl_error rejects invalid autocorr", {
  expect_error(impl_error(autocorr = -0.1))
  expect_error(impl_error(autocorr = 1.1))
})

test_that("impl_error rejects invalid max_overage", {
  expect_error(impl_error(max_overage = 0.9))
})

test_that("impl_error print method works", {
  ie <- impl_error()
  out <- capture.output(print(ie))
  expect_true(any(grepl("Implementation Error", out)))
  expect_true(any(grepl("CV", out)))
  expect_true(any(grepl("unbiased", out)))
})

test_that("impl_error print shows over-catch for bias > 1", {
  ie <- impl_error(bias = 1.1)
  out <- capture.output(print(ie))
  expect_true(any(grepl("over-catch", out)))
})

test_that("validate_impl_error rejects wrong class", {
  expect_error(validate_impl_error(list()), "not of class")
})


# ── mse_scenario ───────────────────────────────────────────────────────────────

test_that("create_scenario creates valid object", {
  my_hcr <- function(biomass, reference_points) biomass * 0.1
  sc <- create_scenario("TestHCR", my_hcr)
  expect_s3_class(sc, "mse_scenario")
  expect_equal(sc$name, "TestHCR")
  expect_true(is.function(sc$harvest_control_rule))
  expect_null(sc$implementation_error)
  expect_equal(sc$assessment_frequency, 1L)
})

test_that("create_scenario accepts implementation_error", {
  my_hcr <- function(biomass, reference_points) biomass * 0.1
  ie <- impl_error(cv = 0.2)
  sc <- create_scenario("WithError", my_hcr,
    implementation_error = ie,
    assessment_frequency = 2
  )
  expect_s3_class(sc$implementation_error, "impl_error")
  expect_equal(sc$assessment_frequency, 2L)
})

test_that("create_scenario rejects non-character name", {
  my_hcr <- function(biomass, reference_points) 100
  expect_error(create_scenario(123, my_hcr))
})

test_that("create_scenario rejects non-function HCR", {
  expect_error(create_scenario("bad", "not_a_function"))
})

test_that("create_scenario rejects HCR with too few arguments", {
  bad_hcr <- function(x) x * 0.1
  expect_error(create_scenario("bad", bad_hcr), "at least 2 arguments")
})

test_that("create_scenario rejects invalid impl_error class", {
  my_hcr <- function(biomass, reference_points) 100
  expect_error(create_scenario("bad", my_hcr, implementation_error = list(cv = 0.1)))
})

test_that("create_scenario rejects invalid assessment_frequency", {
  my_hcr <- function(biomass, reference_points) 100
  expect_error(create_scenario("bad", my_hcr, assessment_frequency = 0))
})

test_that("create_scenario print method works", {
  my_hcr <- function(biomass, reference_points) biomass * 0.1
  sc <- create_scenario("TestHCR", my_hcr)
  out <- capture.output(print(sc))
  expect_true(any(grepl("MSE Scenario", out)))
  expect_true(any(grepl("TestHCR", out)))
  expect_true(any(grepl("none", out)))
})

test_that("create_scenario print shows implementation error info", {
  my_hcr <- function(biomass, reference_points) biomass * 0.1
  ie <- impl_error(cv = 0.15)
  sc <- create_scenario("WithErr", my_hcr, implementation_error = ie)
  out <- capture.output(print(sc))
  expect_true(any(grepl("0.15", out)))
})

test_that("validate_mse_scenario rejects wrong class", {
  expect_error(validate_mse_scenario(list()), "not of class")
})
