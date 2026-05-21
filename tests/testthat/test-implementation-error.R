# Tests for implementation error module (Module 2.2)

# ========================================================================
# Input validation
# ========================================================================

test_that("apply_implementation_error rejects non-impl_error config", {
  expect_error(
    apply_implementation_error(1000, list(cv = 0.1)),
    "impl_error"
  )
})

test_that("apply_implementation_error rejects negative TAC", {
  cfg <- impl_error(cv = 0.1)
  expect_error(apply_implementation_error(-100, cfg), "tac")
})

# ========================================================================
# Perfect implementation (NULL config)
# ========================================================================

test_that("NULL impl_config gives perfect implementation", {
  result <- apply_implementation_error(1000, NULL)
  expect_equal(result$catch, 1000)
  expect_equal(result$eps, 0)
})

test_that("NULL impl_config works with vector TAC", {
  result <- apply_implementation_error(c(500, 600, 700), NULL)
  expect_equal(result$catch, c(500, 600, 700))
  expect_equal(result$eps, c(0, 0, 0))
})

# ========================================================================
# Lognormal distribution (test_implementation_error_distribution)
# ========================================================================

test_that("implementation error produces lognormal-distributed catch", {
  cfg <- impl_error(cv = 0.2, bias = 1, autocorr = 0, max_overage = Inf)
  tac <- 1000
  catches <- replicate(5000, {
    apply_implementation_error(tac, cfg)$catch
  })

  # Log-catches should be approximately normal
  log_catches <- log(catches)
  # Shapiro test on a subsample (max 5000)
  p_val <- shapiro.test(sample(log_catches, 500))$p.value
  expect_true(p_val > 0.01)

  # CV should be close to the configured value
  empirical_cv <- sd(catches) / mean(catches)
  expect_true(abs(empirical_cv - 0.2) < 0.03)
})

test_that("implementation error is reproducible with seed", {
  cfg <- impl_error(cv = 0.15)
  r1 <- apply_implementation_error(1000, cfg, seed = 42)
  r2 <- apply_implementation_error(1000, cfg, seed = 42)
  expect_identical(r1, r2)
})

test_that("different seeds give different results", {
  cfg <- impl_error(cv = 0.15)
  r1 <- apply_implementation_error(1000, cfg, seed = 1)
  r2 <- apply_implementation_error(1000, cfg, seed = 2)
  expect_false(identical(r1$catch, r2$catch))
})

# ========================================================================
# Catch-TAC ratio (test_catch_tac_ratio)
# ========================================================================

test_that("unbiased implementation has E[catch] ~ TAC", {
  cfg <- impl_error(cv = 0.1, bias = 1, autocorr = 0, max_overage = Inf)
  tac <- 1000
  catches <- replicate(10000, {
    apply_implementation_error(tac, cfg)$catch
  })
  # Mean catch should be close to TAC (bias-corrected)
  expect_true(abs(mean(catches) / tac - 1) < 0.02)
})

test_that("biased implementation shifts mean catch", {
  # bias = 1.2 means expected catch = 1.2 * TAC
  cfg <- impl_error(cv = 0.1, bias = 1.2, autocorr = 0, max_overage = Inf)
  tac <- 1000
  catches <- replicate(5000, {
    apply_implementation_error(tac, cfg)$catch
  })
  mean_ratio <- mean(catches) / tac
  expect_true(abs(mean_ratio - 1.2) < 0.05)
})

test_that("under-catch bias works correctly", {
  cfg <- impl_error(cv = 0.1, bias = 0.8, autocorr = 0, max_overage = Inf)
  tac <- 1000
  catches <- replicate(5000, {
    apply_implementation_error(tac, cfg)$catch
  })
  mean_ratio <- mean(catches) / tac
  expect_true(abs(mean_ratio - 0.8) < 0.05)
})

# ========================================================================
# Overage constraints (test_overage_constraints)
# ========================================================================

test_that("max_overage caps realized catch", {
  cfg <- impl_error(cv = 0.3, bias = 1, autocorr = 0, max_overage = 1.1)
  tac <- 1000
  catches <- replicate(5000, {
    apply_implementation_error(tac, cfg)$catch
  })
  # No catch should exceed 1.1 * TAC
  expect_true(all(catches <= 1.1 * tac + 1e-10))
  # Some catches should be at the cap (with CV=0.3, frequent exceedances)
  expect_true(sum(catches >= 1.1 * tac - 1) > 0)
})

test_that("max_overage = 1 allows no over-catch", {
  cfg <- impl_error(cv = 0.2, bias = 1, autocorr = 0, max_overage = 1.0)
  tac <- 1000
  catches <- replicate(2000, {
    apply_implementation_error(tac, cfg)$catch
  })
  expect_true(all(catches <= tac + 1e-10))
})

test_that("catch is always non-negative", {
  cfg <- impl_error(cv = 0.5, bias = 0.5, autocorr = 0, max_overage = 1.1)
  tac <- 100
  catches <- replicate(2000, {
    apply_implementation_error(tac, cfg)$catch
  })
  expect_true(all(catches >= 0))
})

# ========================================================================
# Autocorrelation (test_implementation_error_autocorrelation)
# ========================================================================

test_that("autocorrelation propagates through chained calls", {
  cfg <- impl_error(cv = 0.2, bias = 1, autocorr = 0.8, max_overage = Inf)
  tac <- 1000

  n_years <- 200
  n_reps <- 500
  lag1_cors <- numeric(n_reps)

  for (rep in seq_len(n_reps)) {
    eps_seq <- numeric(n_years)
    prev_eps <- 0
    for (t in seq_len(n_years)) {
      result <- apply_implementation_error(tac, cfg, previous_eps = prev_eps)
      eps_seq[t] <- result$eps
      prev_eps <- result$eps
    }
    lag1_cors[rep] <- cor(eps_seq[-n_years], eps_seq[-1])
  }
  mean_cor <- mean(lag1_cors)
  expect_true(mean_cor > 0.5, label = paste("mean lag-1 cor =", round(mean_cor, 3)))
})

test_that("zero autocorrelation gives uncorrelated errors", {
  cfg <- impl_error(cv = 0.2, bias = 1, autocorr = 0, max_overage = Inf)
  tac <- 1000

  n_years <- 300
  eps_seq <- numeric(n_years)
  prev_eps <- 0
  for (t in seq_len(n_years)) {
    result <- apply_implementation_error(tac, cfg, previous_eps = prev_eps)
    eps_seq[t] <- result$eps
    prev_eps <- result$eps
  }
  lag1_cor <- cor(eps_seq[-n_years], eps_seq[-1])
  expect_true(abs(lag1_cor) < 0.2)
})

# ========================================================================
# Zero-catch / fishery closure (test_zero_catch_scenario)
# ========================================================================

test_that("zero TAC gives zero catch", {
  cfg <- impl_error(cv = 0.2, bias = 1.2, autocorr = 0.5)
  result <- apply_implementation_error(0, cfg, previous_eps = 0.5, seed = 1)
  expect_equal(result$catch, 0)
  expect_equal(result$eps, 0)
})

test_that("zero TAC resets error state", {
  cfg <- impl_error(cv = 0.2, autocorr = 0.8)
  # Build up error state
  r1 <- apply_implementation_error(1000, cfg, previous_eps = 0, seed = 10)
  # Close fishery
  r2 <- apply_implementation_error(0, cfg, previous_eps = r1$eps)
  expect_equal(r2$eps, 0)
  # Reopen — should start fresh (no autocorrelation carried over)
  r3 <- apply_implementation_error(1000, cfg, previous_eps = r2$eps, seed = 10)
  expect_equal(r3$eps, apply_implementation_error(1000, cfg,
    previous_eps = 0,
    seed = 10
  )$eps)
})

# ========================================================================
# Multi-area support
# ========================================================================

test_that("apply_implementation_error works with vector TAC (multi-area)", {
  cfg <- impl_error(cv = 0.1)
  tac <- c(500, 600, 700)
  result <- apply_implementation_error(tac, cfg, seed = 5)
  expect_equal(length(result$catch), 3)
  expect_equal(length(result$eps), 3)
  expect_true(all(result$catch > 0))
})

test_that("multi-area previous_eps is used correctly", {
  cfg <- impl_error(cv = 0.15, autocorr = 0.7, max_overage = Inf)
  tac <- c(500, 600)
  prev <- c(0.5, -0.3)
  result <- apply_implementation_error(tac, cfg, previous_eps = prev, seed = 99)
  expect_equal(length(result$catch), 2)
  expect_equal(length(result$eps), 2)
})

test_that("multi-area with some zero TAC areas", {
  cfg <- impl_error(cv = 0.1, autocorr = 0.5)
  tac <- c(1000, 0, 800)
  result <- apply_implementation_error(tac, cfg,
    previous_eps = c(0.2, 0.3, 0.1),
    seed = 7
  )
  expect_equal(result$catch[2], 0)
  expect_equal(result$eps[2], 0)
  expect_true(result$catch[1] > 0)
  expect_true(result$catch[3] > 0)
})

# ========================================================================
# Edge cases
# ========================================================================

test_that("very small CV gives catch very close to TAC", {
  cfg <- impl_error(cv = 0.001, bias = 1, autocorr = 0, max_overage = Inf)
  tac <- 1000
  result <- apply_implementation_error(tac, cfg, seed = 1)
  expect_true(abs(result$catch - tac) < 5)
})

test_that("previous_eps scalar is recycled for vector TAC", {
  cfg <- impl_error(cv = 0.1, autocorr = 0.5, max_overage = Inf)
  tac <- c(500, 600, 700)
  result <- apply_implementation_error(tac, cfg, previous_eps = 0.2, seed = 3)
  expect_equal(length(result$catch), 3)
})
