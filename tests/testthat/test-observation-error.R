# Tests for observation error module

# --- Helper: build a minimal mock fitted ProductionModel ---
mock_fitted_model <- function(residuals, cpue = NULL) {
  obj <- list(
    fitted = TRUE,
    results = list(residuals = residuals),
    data = list(
      years = seq_len(NROW(residuals)),
      cpue = if (!is.null(cpue)) cpue else rep(1, NROW(residuals))
    )
  )
  class(obj) <- "ProductionModel"
  obj
}

# ========================================================================
# calibrate_observation_error
# ========================================================================

test_that("calibrate_observation_error works with vector residuals", {
  set.seed(42)
  r <- rnorm(50, sd = 0.25)
  fit <- mock_fitted_model(r)
  params <- calibrate_observation_error(fit)

  expect_s3_class(params, "obs_error_params")
  expect_true(is.numeric(params$sigma))
  expect_true(is.numeric(params$rho))
  expect_equal(length(params$sigma), 1)
  expect_null(params$labels)
  # sigma should be close to 0.25 for 50 draws
  expect_true(abs(params$sigma - 0.25) < 0.1)
})

test_that("calibrate_observation_error works with matrix residuals", {
  set.seed(123)
  r <- matrix(rnorm(60, sd = 0.3), nrow = 20, ncol = 3)
  colnames(r) <- c("A1", "A2", "A3")
  fit <- mock_fitted_model(r)
  params <- calibrate_observation_error(fit)

  expect_s3_class(params, "obs_error_params")
  expect_equal(length(params$sigma), 1) # pooled across areas
  expect_null(params$labels)
  expect_true(abs(params$sigma - 0.3) < 0.1)
  expect_equal(nrow(params$by_series), 3)
})

test_that("calibrate_observation_error retains area-specific temporal sequences", {
  # Each area has perfect negative lag-one correlation. Flattening the matrix
  # would insert a non-temporal transition between the two areas.
  r <- cbind(A1 = c(-1, 1, -1), A2 = c(-1, 1, -1))
  params <- calibrate_observation_error(mock_fitted_model(r))

  expect_lt(params$rho, -0.99)
  expect_equal(params$by_series$n_pairs, c(2L, 2L))
})

test_that("calibrate_observation_error does not bridge missing temporal observations", {
  r <- c(-1, NA, 1, -1, 1)
  params <- calibrate_observation_error(mock_fitted_model(r))

  expect_equal(params$by_series$n_pairs, 2L)
})

test_that("calibrate_observation_error works with 3D array residuals (multi-index)", {
  set.seed(99)
  r <- array(rnorm(120, sd = 0.2),
    dim = c(20, 2, 3),
    dimnames = list(NULL, c("A1", "A2"), c("survey", "longline", "pots"))
  )
  fit <- mock_fitted_model(r)
  params <- calibrate_observation_error(fit)

  expect_s3_class(params, "obs_error_params")
  expect_equal(length(params$sigma), 3)
  expect_equal(params$labels, c("survey", "longline", "pots"))
  expect_equal(length(params$rho), 3)
  expect_equal(length(params$n_obs), 3)
  # Each label should have 20*2 = 40 obs
  expect_equal(as.integer(params$n_obs), rep(40L, 3))
})

test_that("calibrate_observation_error handles NAs in residuals", {
  set.seed(7)
  r <- rnorm(30, sd = 0.15)
  r[c(5, 10, 15, 20)] <- NA
  fit <- mock_fitted_model(r)
  params <- calibrate_observation_error(fit)

  expect_equal(params$n_obs, 26)
  expect_true(is.numeric(params$sigma) && is.finite(params$sigma))
})

test_that("calibrate_observation_error rejects non-ProductionModel", {
  expect_error(
    calibrate_observation_error(list(a = 1)),
    "ProductionModel"
  )
})

test_that("calibrate_observation_error rejects unfitted model", {
  fit <- mock_fitted_model(rnorm(10))
  fit$fitted <- FALSE
  expect_error(calibrate_observation_error(fit), "fitted")
})

test_that("calibrate_observation_error returns rho ~ 0 for iid residuals", {
  set.seed(55)
  r <- rnorm(200, sd = 0.1)
  fit <- mock_fitted_model(r)
  params <- calibrate_observation_error(fit)
  expect_true(abs(params$rho) < 0.2)
})

test_that("calibrate_observation_error detects autocorrelation", {
  set.seed(33)
  n <- 200
  eps <- numeric(n)
  rho_true <- 0.7
  for (t in 2:n) eps[t] <- rho_true * eps[t - 1] + rnorm(1, sd = 0.1)
  fit <- mock_fitted_model(eps)
  params <- calibrate_observation_error(fit)
  expect_true(params$rho > 0.4)
})

test_that("print.obs_error_params works for scalar", {
  params <- structure(
    list(sigma = 0.25, rho = 0.1, labels = NULL, n_obs = 20L),
    class = "obs_error_params"
  )
  out <- capture.output(print(params))
  expect_true(any(grepl("sigma", out)))
})

test_that("print.obs_error_params works for multi-index", {
  params <- structure(
    list(
      sigma = c(survey = 0.2, longline = 0.3),
      rho = c(survey = 0.1, longline = 0.05),
      labels = c("survey", "longline"),
      n_obs = c(survey = 40, longline = 40)
    ),
    class = "obs_error_params"
  )
  out <- capture.output(print(params))
  expect_true(any(grepl("survey", out)))
  expect_true(any(grepl("longline", out)))
})

# ========================================================================
# simulate_cpue
# ========================================================================

test_that("simulate_cpue returns correct dimensions for vector biomass", {
  B <- c(5000, 4800, 4600, 4400, 4200)
  result <- simulate_cpue(B, list(sigma = 0.2, rho = 0),
    q = 1e-4,
    years = 2020:2024, seed = 1
  )
  expect_s3_class(result, "data.frame")
  expect_true("year" %in% names(result))
  expect_true("cpue" %in% names(result))
  expect_equal(nrow(result), 5)
  expect_true(all(result$cpue > 0))
})

test_that("simulate_cpue is reproducible with seed", {
  B <- seq(5000, 3000, length.out = 10)
  r1 <- simulate_cpue(B, list(sigma = 0.3, rho = 0.2),
    q = 1e-4,
    years = 2020:2029, seed = 42
  )
  r2 <- simulate_cpue(B, list(sigma = 0.3, rho = 0.2),
    q = 1e-4,
    years = 2020:2029, seed = 42
  )
  expect_identical(r1, r2)
})

test_that("simulate_cpue gives different results with different seeds", {
  B <- rep(5000, 5)
  r1 <- simulate_cpue(B, list(sigma = 0.3, rho = 0),
    q = 1e-4,
    years = 1:5, seed = 1
  )
  r2 <- simulate_cpue(B, list(sigma = 0.3, rho = 0),
    q = 1e-4,
    years = 1:5, seed = 2
  )
  expect_false(identical(r1$cpue, r2$cpue))
})

test_that("simulate_cpue works with matrix biomass (multi-area)", {
  B <- matrix(c(5000, 4500, 4000, 3000, 2800, 2600), nrow = 3, ncol = 2)
  colnames(B) <- c("Ross", "Amundsen")
  result <- simulate_cpue(B, list(sigma = 0.15, rho = 0),
    q = 1e-4,
    years = 2020:2022, seed = 10
  )
  expect_true("area" %in% names(result))
  expect_true(all(c("Ross", "Amundsen") %in% result$area))
  expect_equal(nrow(result), 6) # 3 years x 2 areas
})

test_that("simulate_cpue works with multi-label", {
  B <- c(5000, 4800, 4600)
  result <- simulate_cpue(B, list(sigma = c(0.2, 0.3), rho = c(0, 0.1)),
    q = 1e-4, years = 2020:2022,
    labels = c("survey", "longline"), seed = 5
  )
  expect_true("label" %in% names(result))
  expect_true(all(c("survey", "longline") %in% result$label))
  expect_equal(nrow(result), 6) # 3 years x 2 labels
})

test_that("simulate_cpue applies missing data", {
  B <- rep(5000, 100)
  result <- simulate_cpue(B, list(sigma = 0.1, rho = 0),
    q = 1e-4,
    years = seq_len(100), missing_prob = 0.3, seed = 77
  )
  # Should have fewer than 100 rows (some missing)
  expect_true(nrow(result) < 100)
  expect_true(nrow(result) > 50) # but not nearly all missing
  expect_true(all(!is.na(result$cpue)))
})

test_that("simulate_cpue with zero sigma gives deterministic result", {
  B <- c(5000, 4000, 3000)
  result <- simulate_cpue(B, list(sigma = 0, rho = 0),
    q = 1e-4,
    years = 2020:2022, seed = 1
  )
  expected <- 1e-4 * c(5000, 4000, 3000) * exp(0 - 0)
  expect_equal(result$cpue, expected, tolerance = 1e-10)
})

test_that("simulate_cpue mean CPUE is approximately q*B (bias correction)", {
  # With many realisations, the mean CPUE should be close to q*B
  B <- rep(10000, 1)
  vals <- replicate(5000, {
    simulate_cpue(B, list(sigma = 0.3, rho = 0),
      q = 1e-4,
      years = 2020L, seed = NULL
    )$cpue
  })
  # E[q*B*exp(eps - sigma^2/2)] = q*B = 1.0
  expect_true(abs(mean(vals) - 1.0) < 0.05)
})

test_that("simulate_cpue accepts obs_error_params object", {
  params <- structure(
    list(sigma = 0.2, rho = 0.1, labels = NULL, n_obs = 30L),
    class = "obs_error_params"
  )
  B <- rep(5000, 5)
  result <- simulate_cpue(B, params, q = 1e-4, years = 2020:2024, seed = 99)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
})

test_that("simulate_cpue uses labels from obs_error_params when not explicit", {
  params <- structure(
    list(
      sigma = c(survey = 0.2, longline = 0.3),
      rho = c(survey = 0.0, longline = 0.1),
      labels = c("survey", "longline"),
      n_obs = c(survey = 40, longline = 40)
    ),
    class = "obs_error_params"
  )
  B <- rep(5000, 3)
  result <- simulate_cpue(B, params, q = 1e-4, years = 1:3, seed = 11)
  expect_true("label" %in% names(result))
  expect_true(all(c("survey", "longline") %in% result$label))
})

test_that("simulate_cpue recycles q across areas", {
  B <- matrix(c(5000, 4000, 3000, 6000, 5000, 4000), nrow = 3, ncol = 2)
  result <- simulate_cpue(B, list(sigma = 0.1, rho = 0),
    q = 2e-4,
    years = 1:3, seed = 3
  )
  expect_equal(nrow(result), 6)
})

test_that("simulate_cpue validates year length", {
  B <- c(5000, 4000)
  expect_error(
    simulate_cpue(B, list(sigma = 0.1, rho = 0), q = 1e-4, years = 1:5),
    "years"
  )
})
