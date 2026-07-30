# Tests for performance metrics module (Module 2.4)

# --- Helper: create a simple om_config for metrics ---
make_metrics_config <- function(r = 0.3, K = 5000, m = 2, B_initial = 5000,
                                n_areas = 1) {
  om_config(
    n_areas = n_areas,
    true_params = list(
      r = r, K = K, m = m, sigma_obs = 0.2,
      q = rep(1e-4, n_areas), B_initial = rep(B_initial / n_areas, n_areas)
    )
  )
}

# --- Helper: simulate biomass/catch trajectories ---
# Creates n_sims trajectories with deterministic decline
make_trajectories <- function(n_sims = 100, n_years = 20,
                              B_initial = 5000, decline_rate = 0.02,
                              catch_frac = 0.05, noise_sd = 0.1) {
  set.seed(42)
  biomass <- matrix(NA, n_sims, n_years)
  catch <- matrix(NA, n_sims, n_years)
  for (s in seq_len(n_sims)) {
    B <- B_initial
    for (t in seq_len(n_years)) {
      B <- B * (1 - decline_rate) * exp(rnorm(1, 0, noise_sd))
      B <- max(B, 1)
      biomass[s, t] <- B
      catch[s, t] <- catch_frac * B
    }
  }
  list(biomass = biomass, catch = catch)
}

# ========================================================================
# equilibrium_f
# ========================================================================

test_that("equilibrium_f returns correct Schaefer values", {
  # Standard PT Schaefer (m=2): F_eq = r/(m-1)*(1-(x*B_initial/K)^(m-1)) = r*(1-x)
  r <- 0.3
  K <- 5000
  m <- 2

  f50 <- equilibrium_f(0.5, r, K, m)
  expect_equal(unname(f50), r * (1 - 0.5)) # 0.15

  f20 <- equilibrium_f(0.2, r, K, m)
  expect_equal(unname(f20), r * (1 - 0.2)) # 0.24
})

test_that("equilibrium_f returns 0 when x*B_initial >= K", {
  expect_equal(unname(equilibrium_f(1.0, r = 0.3, K = 5000, m = 2)), 0)
})

test_that("equilibrium_f names output correctly", {
  f <- equilibrium_f(c(0.5, 0.2), r = 0.3, K = 5000, m = 2)
  expect_equal(names(f), c("F50%K", "F20%K"))
})

test_that("equilibrium_f handles B_initial != K", {
  # B_initial = 4000, K = 5000, x = 0.5 => B_target = 2000
  # F = r/(m-1)*(1-(2000/5000)^1) = 0.3*0.6 = 0.18
  f <- equilibrium_f(0.5, r = 0.3, K = 5000, m = 2, B_unfished = 4000)
  expect_equal(unname(f), 0.3 * (1 - (2000 / 5000)), tolerance = 1e-10)
})

test_that("equilibrium_f validates inputs", {
  expect_error(equilibrium_f(-0.1, 0.3, 5000, 2), "x")
  expect_error(equilibrium_f(1.5, 0.3, 5000, 2), "x")
})

# ========================================================================
# calculate_performance_metrics — basic structure
# ========================================================================

test_that("calculate_performance_metrics returns correct class", {
  cfg <- make_metrics_config()
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  expect_s3_class(perf, "mse_performance")
  expect_equal(perf$n_sims, 100)
  expect_equal(perf$n_years, 20)
  expect_equal(perf$n_areas, 1)
})

test_that("result contains all expected components", {
  cfg <- make_metrics_config()
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  expect_true("biomass_risk" %in% names(perf))
  expect_true("f_risk" %in% names(perf))
  expect_true("catch_stats" %in% names(perf))
  expect_true("biomass_ratios" %in% names(perf))
  expect_true("reference" %in% names(perf))
})

test_that("default thresholds produce 50% and 20% entries", {
  cfg <- make_metrics_config()
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  expect_true("50%K" %in% names(perf$biomass_risk))
  expect_true("20%K" %in% names(perf$biomass_risk))
  expect_true("F50%K" %in% names(perf$f_risk))
  expect_true("F20%K" %in% names(perf$f_risk))
})

# ========================================================================
# Biomass risk metrics
# ========================================================================

test_that("Pr(B < threshold) per_year has correct length", {
  cfg <- make_metrics_config()
  traj <- make_trajectories(n_years = 15)
  perf <- calculate_performance_metrics(traj, cfg)

  per_year <- perf$biomass_risk[["50%K"]]$aggregate$per_year
  expect_equal(length(per_year), 15)
})

test_that("biomass always above threshold gives zero risk", {
  cfg <- make_metrics_config(B_initial = 5000)
  # All biomass stays at B_initial = 5000, threshold at 50% = 2500
  n_sims <- 50
  n_years <- 10
  biomass <- matrix(5000, n_sims, n_years)
  catch <- matrix(0, n_sims, n_years)
  traj <- list(biomass = biomass, catch = catch)
  perf <- calculate_performance_metrics(traj, cfg)

  risk50 <- perf$biomass_risk[["50%K"]]$aggregate
  expect_equal(risk50$final_year, 0)
  expect_equal(risk50$ever, 0)
  expect_true(all(risk50$per_year == 0))
})

test_that("biomass always below threshold gives risk = 1", {
  cfg <- make_metrics_config(B_initial = 5000)
  n_sims <- 50
  n_years <- 10
  # Biomass at 10% of B_initial = 500, below both 50% and 20% thresholds
  biomass <- matrix(500, n_sims, n_years)
  catch <- matrix(10, n_sims, n_years)
  traj <- list(biomass = biomass, catch = catch)
  perf <- calculate_performance_metrics(traj, cfg)

  risk50 <- perf$biomass_risk[["50%K"]]$aggregate
  expect_equal(risk50$final_year, 1)
  expect_equal(risk50$ever, 1)

  risk20 <- perf$biomass_risk[["20%K"]]$aggregate
  expect_equal(risk20$final_year, 1)
  expect_equal(risk20$ever, 1)
})

test_that("'ever' probability >= final_year probability", {
  cfg <- make_metrics_config()
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  for (nm in names(perf$biomass_risk)) {
    agg <- perf$biomass_risk[[nm]]$aggregate
    expect_true(agg$ever >= agg$final_year - 1e-10)
  }
})

# ========================================================================
# F risk metrics
# ========================================================================

test_that("F risk uses correct equilibrium threshold", {
  cfg <- make_metrics_config(r = 0.3, K = 5000, m = 2)
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  # F50%K for standard PT Schaefer = 0.3*(1-0.5) = 0.15
  expect_equal(unname(perf$reference$f_thresholds["F50%K"]),
    0.15,
    tolerance = 1e-10
  )
  # F20%K = 0.3*(1-0.2) = 0.24
  expect_equal(unname(perf$reference$f_thresholds["F20%K"]),
    0.24,
    tolerance = 1e-10
  )
})

test_that("F risk structure matches biomass risk", {
  cfg <- make_metrics_config()
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  for (nm in names(perf$f_risk)) {
    agg <- perf$f_risk[[nm]]$aggregate
    expect_true("per_year" %in% names(agg))
    expect_true("final_year" %in% names(agg))
    expect_true("ever" %in% names(agg))
  }
})

test_that("zero harvest gives zero F risk", {
  cfg <- make_metrics_config()
  biomass <- matrix(3000, 50, 10)
  catch <- matrix(0, 50, 10)
  traj <- list(biomass = biomass, catch = catch)
  perf <- calculate_performance_metrics(traj, cfg)

  for (nm in names(perf$f_risk)) {
    expect_equal(perf$f_risk[[nm]]$aggregate$ever, 0)
  }
})

# ========================================================================
# Catch statistics
# ========================================================================

test_that("catch stats are computed correctly", {
  cfg <- make_metrics_config()
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  cs <- perf$catch_stats$aggregate
  expect_true(cs$mean_catch > 0)
  expect_true(cs$median_catch > 0)
  expect_true(cs$sd_catch > 0)
  expect_true(cs$aav >= 0)
})

test_that("AAV is zero for constant catch", {
  cfg <- make_metrics_config()
  biomass <- matrix(3000, 50, 10)
  catch <- matrix(100, 50, 10)
  traj <- list(biomass = biomass, catch = catch)
  perf <- calculate_performance_metrics(traj, cfg)

  expect_equal(perf$catch_stats$aggregate$aav, 0)
})

test_that("AAV increases with catch variability", {
  cfg <- make_metrics_config()
  set.seed(1)
  n_sims <- 50
  n_years <- 20

  # Low variability
  catch_low <- matrix(
    100 + rnorm(n_sims * n_years, 0, 1),
    n_sims, n_years
  )
  traj_low <- list(
    biomass = matrix(3000, n_sims, n_years),
    catch = pmax(catch_low, 0)
  )

  # High variability
  catch_high <- matrix(
    100 + rnorm(n_sims * n_years, 0, 50),
    n_sims, n_years
  )
  traj_high <- list(
    biomass = matrix(3000, n_sims, n_years),
    catch = pmax(catch_high, 1)
  )

  aav_low <- calculate_performance_metrics(traj_low, cfg)$catch_stats$aggregate$aav
  aav_high <- calculate_performance_metrics(traj_high, cfg)$catch_stats$aggregate$aav
  expect_true(aav_high > aav_low)
})

# ========================================================================
# Biomass ratios
# ========================================================================

test_that("biomass ratios are correct for known values", {
  cfg <- make_metrics_config(r = 0.3, K = 5000, m = 2, B_initial = 5000)
  # BMSY for Schaefer = K/2 = 2500
  biomass <- matrix(2500, 50, 10) # exactly at BMSY
  catch <- matrix(0, 50, 10)
  traj <- list(biomass = biomass, catch = catch)
  perf <- calculate_performance_metrics(traj, cfg)

  expect_equal(perf$biomass_ratios$aggregate$mean_depletion, 0.5,
    tolerance = 1e-10
  )
  expect_equal(perf$biomass_ratios$aggregate$mean_b_bmsy, 1.0,
    tolerance = 1e-10
  )
})

test_that("reference points are correct for Schaefer", {
  cfg <- make_metrics_config(r = 0.3, K = 5000, m = 2)
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  expect_equal(perf$reference$BMSY, 2500, tolerance = 1e-6)
  # Standard PT: FMSY = r/m = 0.3/2 = 0.15
  expect_equal(perf$reference$FMSY, 0.15, tolerance = 1e-6)
  # MSY = FMSY * BMSY = 0.15 * 2500 = 375
  expect_equal(perf$reference$MSY, 375, tolerance = 1e-6)
})

# ========================================================================
# Multi-area support
# ========================================================================

test_that("multi-area metrics are computed correctly", {
  cfg <- make_metrics_config(K = 6000, B_initial = 6000, n_areas = 3)
  # K_area = c(2000, 2000, 2000)
  n_sims <- 30
  n_years <- 10
  n_areas <- 3

  biomass <- array(NA, dim = c(n_sims, n_years, n_areas))
  catch <- array(NA, dim = c(n_sims, n_years, n_areas))
  for (a in 1:n_areas) {
    biomass[, , a] <- matrix(2000 - (a - 1) * 500, n_sims, n_years)
    catch[, , a] <- matrix(50, n_sims, n_years)
  }
  traj <- list(biomass = biomass, catch = catch)
  perf <- calculate_performance_metrics(traj, cfg)

  expect_equal(perf$n_areas, 3)

  # Check by_area exists for biomass risk
  for (nm in names(perf$biomass_risk)) {
    expect_equal(length(perf$biomass_risk[[nm]]$by_area), 3)
    expect_equal(
      names(perf$biomass_risk[[nm]]$by_area),
      c("A1", "A2", "A3")
    )
  }

  # Check by_area exists for F risk
  for (nm in names(perf$f_risk)) {
    expect_equal(length(perf$f_risk[[nm]]$by_area), 3)
  }

  # Check by_area catch stats
  expect_equal(length(perf$catch_stats$by_area), 3)
})

test_that("aggregate biomass sums across areas correctly", {
  cfg <- make_metrics_config(K = 6000, B_initial = 6000, n_areas = 2)
  # K_area = c(3000, 3000), threshold 50% = 3000 total
  n_sims <- 50
  n_years <- 5

  biomass <- array(NA, dim = c(n_sims, n_years, 2))
  biomass[, , 1] <- 1000 # area 1
  biomass[, , 2] <- 1500 # area 2
  # Total = 2500 < 50%*6000 = 3000 => risk should be 1
  catch <- array(0, dim = c(n_sims, n_years, 2))
  traj <- list(biomass = biomass, catch = catch)
  perf <- calculate_performance_metrics(traj, cfg)

  expect_equal(perf$biomass_risk[["50%K"]]$aggregate$final_year, 1)
})

# ========================================================================
# Custom thresholds
# ========================================================================

test_that("custom thresholds are respected", {
  cfg <- make_metrics_config()
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg, thresholds = c(0.75, 0.3, 0.1))

  expect_equal(length(perf$biomass_risk), 3)
  expect_true("75%K" %in% names(perf$biomass_risk))
  expect_true("30%K" %in% names(perf$biomass_risk))
  expect_true("10%K" %in% names(perf$biomass_risk))
})

# ========================================================================
# Input validation
# ========================================================================

test_that("rejects missing biomass", {
  cfg <- make_metrics_config()
  expect_error(
    calculate_performance_metrics(list(catch = matrix(1, 10, 5)), cfg),
    "biomass"
  )
})

test_that("rejects missing catch", {
  cfg <- make_metrics_config()
  expect_error(
    calculate_performance_metrics(list(biomass = matrix(1, 10, 5)), cfg),
    "catch"
  )
})

test_that("rejects non-om_config", {
  expect_error(
    calculate_performance_metrics(
      list(
        biomass = matrix(1, 10, 5),
        catch = matrix(1, 10, 5)
      ),
      list()
    ),
    "om_config"
  )
})

# ========================================================================
# print and summary methods
# ========================================================================

test_that("print.mse_performance produces output", {
  cfg <- make_metrics_config()
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  out <- capture.output(print(perf))
  expect_true(any(grepl("MSE Performance", out)))
  expect_true(any(grepl("Biomass Risk", out)))
  expect_true(any(grepl("Exploitation-rate Risk", out)))
  expect_true(any(grepl("Catch", out)))
})

test_that("summary.mse_performance returns a data frame", {
  cfg <- make_metrics_config()
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  s <- summary(perf)
  expect_s3_class(s, "data.frame")
  expect_true("metric" %in% names(s))
  expect_true("scope" %in% names(s))
  expect_true("value" %in% names(s))
  expect_true(nrow(s) > 0)
})

test_that("summary contains expected metric names", {
  cfg <- make_metrics_config()
  traj <- make_trajectories()
  perf <- calculate_performance_metrics(traj, cfg)

  s <- summary(perf)
  expect_true(any(grepl("Pr\\(B<50%K\\)", s$metric)))
  expect_true(any(grepl("Pr\\(B<20%K\\)", s$metric)))
  expect_true(any(grepl("Pr\\(U>U50%B0\\)", s$metric)))
  expect_true(any(grepl("mean_catch", s$metric)))
  expect_true(any(grepl("AAV", s$metric)))
  expect_true(any(grepl("mean_B_BMSY", s$metric)))
})

# ========================================================================
# Harvest rate provided vs computed
# ========================================================================

test_that("provided harvest_rate is used when available", {
  cfg <- make_metrics_config()
  biomass <- matrix(3000, 50, 10)
  catch <- matrix(100, 50, 10)
  # Provide explicit F = 0.5 (much higher than catch/biomass ~ 0.033)
  hr <- matrix(0.5, 50, 10)
  traj <- list(biomass = biomass, catch = catch, harvest_rate = hr)
  perf <- calculate_performance_metrics(traj, cfg)

  # F50%K = 0.075, so Pr(F > 0.075) should be 1 with F=0.5
  expect_equal(perf$f_risk[["F50%K"]]$aggregate$final_year, 1)
})
