# Tests for MSE Simulation Loop (Module 2.5)

# ========================================================================
# Helpers
# ========================================================================

make_simple_om <- function() {
  om_config(
    n_areas = 1L,
    true_params = list(
      r = 0.3, K = 5000, m = 2,
      sigma_obs = 0.2, q = 1e-4, B_initial = 5000
    )
  )
}

make_spatial_om <- function() {
  dm <- matrix(c(0, 1, 1, 0), 2, 2)
  om_config(
    n_areas = 2L,
    movement_rate = 0.1,
    movement_cost_matrix = dm,
    true_params = list(
      r = 0.3, K = 5000, m = 2,
      sigma_obs = 0.2, q = c(1e-4, 1e-4),
      B_initial = c(3000, 2000)
    )
  )
}

make_scenario <- function(name = "test", assess_freq = 1L, impl = NULL) {
  create_scenario(
    name = name,
    harvest_control_rule = hcr_constant_f(0.05),
    implementation_error = impl,
    assessment_frequency = assess_freq
  )
}

make_stable_em <- function() {
  em_config(
    fixed_params = list(m = 2, sigma_obs = 0.2),
    control = list(eval.max = 2000, iter.max = 1000),
    n_starts = 3L,
    calculate_se = FALSE,
    priors = list(
      r = list(dist = "lognormal", meanlog = log(0.3), sdlog = 0.5),
      K = list(dist = "lognormal", meanlog = log(5000), sdlog = 0.5)
    ),
    initial_depletion = 1
  )
}

# ========================================================================
# Input validation
# ========================================================================

test_that("mse_simulation rejects invalid operating_model", {
  expect_error(
    mse_simulation(initial_tac = 0, 
      operating_model = "not_om",
      scenarios = make_scenario()
    ),
    "om_config"
  )
})

test_that("mse_simulation rejects om_config without true_params", {
  om <- om_config(n_areas = 1L)
  expect_error(
    mse_simulation(initial_tac = 0, 
      operating_model = om,
      scenarios = make_scenario()
    ),
    "true_params"
  )
})

test_that("mse_simulation rejects invalid estimation_model", {
  expect_error(
    mse_simulation(initial_tac = 0, 
      operating_model = make_simple_om(),
      estimation_model = "bad",
      scenarios = make_scenario()
    ),
    "em_config"
  )
})

test_that("mse_simulation rejects invalid scenarios", {
  expect_error(
    mse_simulation(initial_tac = 0, 
      operating_model = make_simple_om(),
      scenarios = list("bad")
    ),
    "mse_scenario"
  )
})

test_that("mse_simulation rejects zero n_sims", {
  expect_error(
    mse_simulation(initial_tac = 0, 
      operating_model = make_simple_om(),
      scenarios = make_scenario(),
      n_sims = 0
    ),
    "n_sims"
  )
})

# ========================================================================
# Simulation loop runs without error (no EM fitting)
# ========================================================================

test_that("simulation loop completes with no assessments", {
  om <- make_simple_om()
  sc <- make_scenario()
  # min_assess_years > n_proj_years → no EM fit

  res <- mse_simulation(initial_tac = 0, 
    operating_model = om,
    scenarios = sc,
    n_sims = 5L,
    n_proj_years = 10L,
    min_assess_years = 99L,
    seed = 1
  )

  expect_s3_class(res, "mse_result")
  expect_equal(res$n_sims, 5L)
  expect_equal(res$n_proj_years, 10L)
  expect_length(res$results, 1L)
  expect_true("test" %in% names(res$results))
})

test_that("simulation loop accepts single scenario (not wrapped in list)", {
  om <- make_simple_om()
  sc <- make_scenario()

  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 3L, n_proj_years = 5L,
    min_assess_years = 99L, seed = 1
  )
  expect_s3_class(res, "mse_result")
})

test_that("simulation loop runs with multiple scenarios", {
  om <- make_simple_om()
  sc1 <- make_scenario("low_f")
  sc2 <- create_scenario("high_f", hcr_constant_f(0.1))

  res <- mse_simulation(initial_tac = 0, 
    om,
    scenarios = list(sc1, sc2),
    n_sims = 3L, n_proj_years = 8L,
    min_assess_years = 99L, seed = 1
  )

  expect_length(res$results, 2L)
  expect_true(all(c("low_f", "high_f") %in% names(res$results)))
})

# ========================================================================
# Trajectory storage
# ========================================================================

test_that("trajectory arrays have correct dimensions (single area)", {
  om <- make_simple_om()
  sc <- make_scenario()

  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 4L,
    n_proj_years = 8L, min_assess_years = 99L, seed = 1
  )

  traj <- res$results$test$trajectories
  expect_equal(dim(traj$biomass), c(4, 8, 1))
  expect_equal(dim(traj$catch), c(4, 8, 1))
  expect_equal(dim(traj$harvest_rate), c(4, 8, 1))
  expect_equal(dim(traj$tac), c(4, 8))
  expect_equal(dim(traj$estimated_biomass), c(4, 8))
})

test_that("trajectory arrays have correct dimensions (multi-area)", {
  om <- make_spatial_om()
  sc <- make_scenario()

  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 3L,
    n_proj_years = 6L, min_assess_years = 99L, seed = 1
  )

  traj <- res$results$test$trajectories
  expect_equal(dim(traj$biomass), c(3, 6, 2))
  expect_equal(dim(traj$catch), c(3, 6, 2))
  expect_equal(dim(traj$harvest_rate), c(3, 6, 2))
  expect_equal(dim(traj$tac), c(3, 6))
})

test_that("trajectory biomass starts near B_initial", {
  om <- make_simple_om()
  sc <- make_scenario()

  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 10L,
    n_proj_years = 5L, min_assess_years = 99L, seed = 1
  )

  traj <- res$results$test$trajectories
  # Year 1 biomass should equal B_initial (5000) for all sims

  expect_equal(as.numeric(traj$biomass[, 1, 1]), rep(5000, 10))
})

test_that("trajectory catch is positive", {
  om <- make_simple_om()
  sc <- make_scenario()

  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 5L,
    n_proj_years = 8L, min_assess_years = 99L, seed = 1
  )

  traj <- res$results$test$trajectories
  expect_true(all(traj$catch >= 0))
})

test_that("trajectory harvest_rate equals catch / biomass", {
  om <- make_simple_om()
  sc <- make_scenario()

  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 5L,
    n_proj_years = 8L, min_assess_years = 99L, seed = 1
  )

  traj <- res$results$test$trajectories
  expected_hr <- traj$catch / pmax(traj$biomass, 1e-8)
  expect_equal(traj$harvest_rate, expected_hr)
})

test_that("no estimated_biomass before first assessment", {
  om <- make_simple_om()
  sc <- make_scenario()

  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 3L,
    n_proj_years = 5L, min_assess_years = 99L, seed = 1
  )

  traj <- res$results$test$trajectories
  # No assessment → all estimated_biomass are NA
  expect_true(all(is.na(traj$estimated_biomass)))
})

# ========================================================================
# Initial TAC and TAC behaviour
# ========================================================================

test_that("custom initial_tac is used", {
  om <- make_simple_om()
  sc <- make_scenario()

  res <- mse_simulation(om,
    scenarios = sc, n_sims = 3L,
    n_proj_years = 5L, min_assess_years = 99L,
    initial_tac = 100, seed = 1
  )

  traj <- res$results$test$trajectories
  expect_equal(traj$tac[1, 1], 100)
})

test_that("TAC is constant between assessments (interim catch)", {
  om <- make_simple_om()
  # Assessment every 3 years, but first assessment needs 3+ years of data
  sc <- create_scenario("biennial", hcr_constant_f(0.05),
    assessment_frequency = 3L
  )

  res <- mse_simulation(om,
    estimation_model = make_stable_em(),
    scenarios = sc, n_sims = 3L,
    n_proj_years = 12L, min_assess_years = 3L,
    initial_tac = 100, seed = 42
  )

  traj <- res$results$biennial$trajectories
  # Years 1-2: initial_tac = 100 (yr 1,2 are before first assessment at yr=3)
  expect_equal(traj$tac[1, 1], 100)
  expect_equal(traj$tac[1, 2], 100)

  # After first assessment at yr=3, TAC changes, then stays constant until yr=6
  tac_yr3 <- traj$tac[1, 3]
  expect_equal(traj$tac[1, 4], tac_yr3)
  expect_equal(traj$tac[1, 5], tac_yr3)
})

# ========================================================================
# Assessment frequency
# ========================================================================

test_that("assessment_frequency = 1 updates TAC every year after min_assess", {
  om <- make_simple_om()
  sc <- create_scenario("annual", hcr_constant_f(0.05),
    assessment_frequency = 1L
  )

  res <- mse_simulation(om,
    estimation_model = make_stable_em(),
    scenarios = sc, n_sims = 2L,
    n_proj_years = 10L, min_assess_years = 3L,
    initial_tac = 100, seed = 42
  )

  traj <- res$results$annual$trajectories
  # Before min_assess: TAC = 100
  expect_equal(traj$tac[1, 1], 100)
  expect_equal(traj$tac[1, 2], 100)

  # After min_assess (yr >= 3 when assess_freq divides yr):
  # TAC should change from initial at some point
  # (EM fits at yr=3 with 2 data points... actually min_assess=3 means
  #  nrow(cpue_history) >= 3 at yr=4 → first assessment at yr=4)
  # Years 1-3: initial_tac; Year 4+: EM-driven TAC
  tac_later <- traj$tac[1, 4:10]
  # TAC should vary (EM updates each year from yr=4 on)
  # At minimum, it should differ from the initial_tac at some point
  expect_true(any(tac_later != 100))
})

test_that("assessment_frequency = 2 only assesses on even years", {
  om <- make_simple_om()
  sc <- create_scenario("biennial", hcr_constant_f(0.05),
    assessment_frequency = 2L
  )

  res <- mse_simulation(om,
    estimation_model = make_stable_em(),
    scenarios = sc, n_sims = 2L,
    n_proj_years = 10L, min_assess_years = 3L,
    initial_tac = 100, seed = 42
  )

  traj <- res$results$biennial$trajectories
  # 1st assess: yr divisible by 2 with >= 3 data points
  # yr=2: 1 data point (yr=1). No.
  # yr=4: 3 data points (yr 1,2,3). Yes!
  # yr=6: Yes. yr=8: Yes. yr=10: Yes.
  # So TAC changes at yr=4, stays same at yr=5, changes at yr=6, etc.
  tac_yr4 <- traj$tac[1, 4]
  tac_yr5 <- traj$tac[1, 5]
  expect_equal(tac_yr4, tac_yr5) # no assessment at yr=5

  tac_yr6 <- traj$tac[1, 6]
  tac_yr7 <- traj$tac[1, 7]
  expect_equal(tac_yr6, tac_yr7) # no assessment at yr=7
})

# ========================================================================
# Implementation error integration
# ========================================================================

test_that("perfect implementation: catch equals TAC", {
  om <- make_simple_om()
  sc <- create_scenario("perfect", hcr_constant_catch(150))

  res <- mse_simulation(om,
    scenarios = sc, n_sims = 5L,
    n_proj_years = 5L, min_assess_years = 99L,
    initial_tac = 150, seed = 1
  )

  traj <- res$results$perfect$trajectories
  catch_total <- apply(traj$catch, c(1, 2), sum)
  expect_equal(catch_total, traj$tac, tolerance = 1e-10)
})

test_that("implementation error makes catch differ from TAC", {
  om <- make_simple_om()
  impl <- impl_error(cv = 0.2, bias = 1, autocorr = 0)
  sc <- create_scenario("noisy", hcr_constant_catch(150),
    implementation_error = impl
  )

  res <- mse_simulation(om,
    scenarios = sc, n_sims = 10L,
    n_proj_years = 5L, min_assess_years = 99L,
    initial_tac = 150, seed = 1
  )

  traj <- res$results$noisy$trajectories
  catch_total <- apply(traj$catch, c(1, 2), sum)
  # Catch should not exactly equal TAC due to implementation error
  expect_false(all(catch_total == traj$tac))
})

# ========================================================================
# Convergence handling
# ========================================================================

test_that("an assessment is not applied before sufficient data are available", {
  om <- make_simple_om()
  # Use a very low F so biomass barely changes → EM may struggle
  sc <- create_scenario("tricky", hcr_constant_f(0.001),
    assessment_frequency = 1L
  )

  # This construction does not guarantee a failed fit. The hard-error policy
  # is tested directly through the missing-EM-reference test below.
  expect_no_error(
    mse_simulation(initial_tac = 0, om,
      estimation_model = make_stable_em(),
      scenarios = sc, n_sims = 2L,
      n_proj_years = 8L, min_assess_years = 3L,
      seed = 42
    )
  )
})

test_that("initial TAC is retained when no assessment is scheduled", {
  om <- make_simple_om()
  # If all EM fits fail, TAC uses true biomass + true ref points
  # We can't force failure easily, but we can verify the loop
  # completes with reasonable TAC values
  sc <- make_scenario()

  res <- mse_simulation(om,
    scenarios = sc, n_sims = 3L,
    n_proj_years = 5L, min_assess_years = 99L, initial_tac = 375, seed = 1
  )

  traj <- res$results$test$trajectories
  # No assessment ever occurs, so every TAC remains the supplied initial_tac.
  expect_equal(unique(as.vector(traj$tac)), 375)
})

# ========================================================================
# Seed reproducibility
# ========================================================================

test_that("same seed produces identical results", {
  om <- make_simple_om()
  sc <- make_scenario()

  res1 <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 5L,
    n_proj_years = 8L, min_assess_years = 99L, seed = 123
  )
  res2 <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 5L,
    n_proj_years = 8L, min_assess_years = 99L, seed = 123
  )

  expect_equal(
    res1$results$test$trajectories$biomass,
    res2$results$test$trajectories$biomass
  )
  expect_equal(
    res1$results$test$trajectories$catch,
    res2$results$test$trajectories$catch
  )
})

test_that("different seeds produce different results", {
  om <- make_simple_om()
  impl <- impl_error(cv = 0.2)
  sc <- create_scenario("test", hcr_constant_f(0.05),
    implementation_error = impl
  )

  res1 <- mse_simulation(initial_tac = 200, om,
    scenarios = sc, n_sims = 5L,
    n_proj_years = 8L, min_assess_years = 99L, seed = 1
  )
  res2 <- mse_simulation(initial_tac = 200, om,
    scenarios = sc, n_sims = 5L,
    n_proj_years = 8L, min_assess_years = 99L, seed = 2
  )

  expect_false(identical(
    res1$results$test$trajectories$catch,
    res2$results$test$trajectories$catch
  ))
})

# ========================================================================
# Performance metrics integration
# ========================================================================

test_that("performance metrics are calculated for each scenario", {
  om <- make_simple_om()
  sc1 <- make_scenario("s1")
  sc2 <- create_scenario("s2", hcr_constant_f(0.1))

  res <- mse_simulation(initial_tac = 0, om,
    scenarios = list(sc1, sc2), n_sims = 5L,
    n_proj_years = 8L, min_assess_years = 99L, seed = 1
  )

  expect_s3_class(res$results$s1$performance, "mse_performance")
  expect_s3_class(res$results$s2$performance, "mse_performance")
})

# ========================================================================
# Multi-area operation
# ========================================================================

test_that("multi-area simulation runs without error", {
  om <- make_spatial_om()
  sc <- make_scenario()

  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 3L,
    n_proj_years = 6L, min_assess_years = 99L, seed = 1
  )

  expect_s3_class(res, "mse_result")
  traj <- res$results$test$trajectories
  expect_equal(dim(traj$biomass)[3], 2L)
  # Per-area biomass at year 1 should match B_initial
  expect_equal(traj$biomass[1, 1, ], c(3000, 2000))
})

test_that("multi-area TAC is distributed proportional to biomass", {
  om <- make_spatial_om()
  sc <- create_scenario("const", hcr_constant_catch(100))

  res <- mse_simulation(om,
    scenarios = sc, n_sims = 3L,
    n_proj_years = 3L, min_assess_years = 99L,
    initial_tac = 100, seed = 1
  )

  traj <- res$results$const$trajectories
  # Year 1: B = (3000, 2000), total = 5000
  # TAC area = 100 * (3000/5000, 2000/5000) = (60, 40)
  # With perfect implementation:
  expect_equal(traj$catch[1, 1, 1], 60)
  expect_equal(traj$catch[1, 1, 2], 40)
})

# ========================================================================
# Observation error params
# ========================================================================

test_that("obs_error_params is constructed from true_params when NULL", {
  om <- make_simple_om()
  sc <- make_scenario()

  # Should not error — obs_error_params built internally
  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 3L,
    n_proj_years = 5L, min_assess_years = 99L, seed = 1
  )
  expect_s3_class(res, "mse_result")
})

test_that("custom obs_error_params is accepted", {
  om <- make_simple_om()
  sc <- make_scenario()

  oep <- structure(
    list(sigma = 0.1, rho = 0.5, labels = NULL, n_obs = NA_integer_),
    class = "obs_error_params"
  )

  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 3L,
    n_proj_years = 5L, min_assess_years = 99L,
    obs_error_params = oep, seed = 1
  )
  expect_s3_class(res, "mse_result")
})

# ========================================================================
# Print method
# ========================================================================

test_that("print.mse_result works", {
  om <- make_simple_om()
  sc <- make_scenario()
  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 3L,
    n_proj_years = 5L, min_assess_years = 99L, seed = 1
  )

  out <- capture.output(print(res))
  expect_true(any(grepl("MSE Simulation Result", out)))
  expect_true(any(grepl("3", out))) # n_sims
  expect_true(any(grepl("5", out))) # n_proj_years
})

# ========================================================================
# EM-fitting integration tests (small scale)
# ========================================================================

test_that("simulation with EM fitting completes", {
  skip_on_cran()
  om <- make_simple_om()
  sc <- create_scenario("with_em", hcr_constant_f(0.05),
    assessment_frequency = 1L
  )

  # Small sim with EM fitting enabled (min_assess = 5)
  expect_no_error(
    res <- mse_simulation(initial_tac = 0, 
      om,
      estimation_model = make_stable_em(),
      scenarios = sc, n_sims = 2L,
      n_proj_years = 10L, min_assess_years = 5L,
      seed = 42
    )
  )

  expect_s3_class(res, "mse_result")
  traj <- res$results$with_em$trajectories
  # Loop completes and produces finite biomass throughout
  expect_true(all(is.finite(traj$biomass)))
  expect_true(all(is.finite(traj$catch)))
})

test_that("EM updates TAC after first assessment", {
  skip_on_cran()
  om <- make_simple_om()
  sc <- create_scenario("em_update", hcr_constant_f(0.05),
    assessment_frequency = 1L
  )

  res <- mse_simulation(
    om,
    estimation_model = make_stable_em(),
    scenarios = sc, n_sims = 2L,
    n_proj_years = 10L, min_assess_years = 5L,
    initial_tac = 100, seed = 42
  )

  traj <- res$results$em_update$trajectories
  # Before first assessment, TAC = 100
  expect_equal(traj$tac[1, 1], 100)
  # After assessment, TAC should change
  # First assessment at yr=5 (when 4 data points available... no, at yr=6)
  # Actually: yr=1..5 data accumulated before yr=6 check
  # yr=6: nrow(cpue_history)=5 >= 5, yr%%1==0 → assess
  tac_after <- traj$tac[1, 7]
  expect_true(tac_after != 100)
})

# ========================================================================
# Stored call
# ========================================================================

test_that("mse_result stores the matched call", {
  om <- make_simple_om()
  sc <- make_scenario()
  res <- mse_simulation(initial_tac = 0, om,
    scenarios = sc, n_sims = 3L,
    n_proj_years = 5L, min_assess_years = 99L, seed = 1
  )
  expect_true(is.call(res$call))
})

# ========================================================================
# Internal helpers
# ========================================================================

test_that(".build_hcr_ref_points with EM result", {
  em_res <- list(K = 5000, msy = 200, bmsy = 2500, fmsy = 0.08)
  tp <- list(r = 0.3, K = 5000, m = 2)
  ref <- SurplusProductionModelMSE:::.build_hcr_ref_points(em_res, tp)
  expect_equal(ref$K, 5000)
  expect_equal(ref$MSY, 200)
})

test_that(".build_hcr_ref_points requires an EM result", {
  tp <- list(r = 0.3, K = 5000, m = 2)
  expect_error(
    SurplusProductionModelMSE:::.build_hcr_ref_points(NULL, tp),
    "successful estimation-model result"
  )
})

test_that("historical EM data require CPUE years contained in the catch history", {
  expect_error(
    SurplusProductionModelMSE:::.validate_historical_data(list(
      cpue_data = data.frame(year = 2018:2019, cpue = c(0.4, 0.3)),
      catch_data = data.frame(year = 2017:2018, catch = c(200, 220))
    )),
    "years must be contained"
  )
})

test_that("historical EM data initialise the first projected HCR decision", {
  om <- make_simple_om()
  historical_data <- list(
    cpue_data = data.frame(year = 2018:2022, cpue = c(0.45, 0.43, 0.41, 0.39, 0.37)),
    catch_data = data.frame(year = 2018:2022, catch = c(200, 220, 240, 260, 280))
  )
  res <- mse_simulation(
    om,
    estimation_model = make_stable_em(),
    scenarios = make_scenario(),
    n_sims = 1L,
    n_proj_years = 3L,
    initial_tac = 300,
    historical_data = historical_data,
    min_assess_years = 3L,
    seed = 2
  )

  expect_equal(res$projection_years, 2023:2025)
  expect_equal(res$results$test$trajectories$em_fit_success[1, 1], 1)
  expect_true(is.finite(res$results$test$trajectories$tac[1, 1]))
})

test_that("a converged EM fit supplies a valid sequential warm start", {
  skip_on_cran()
  historical_cpue <- data.frame(
    year = 2018:2022,
    cpue = c(0.45, 0.43, 0.41, 0.39, 0.37)
  )
  historical_catch <- data.frame(
    year = 2018:2022,
    catch = c(200, 220, 240, 260, 280)
  )
  em <- make_stable_em()
  first_fit <- SurplusProductionModelMSE:::.fit_em_safely(
    historical_cpue,
    historical_catch,
    em,
    make_simple_om()
  )

  expect_true(all(c(
    "log_r", "log_K.A1", "log_m", "log_q.A1", "log_d0"
  ) %in% names(first_fit$params_init)))
  expect_no_error(
    SurplusProductionModelMSE:::.fit_em_safely(
      rbind(historical_cpue, data.frame(year = 2023, cpue = 0.35)),
      rbind(historical_catch, data.frame(year = 2023, catch = 300)),
      em,
      make_simple_om(),
      params_init = first_fit$params_init
    )
  )
})

test_that("no exploitation cap is applied unless explicitly configured", {
  om <- make_simple_om()
  sc <- make_scenario()
  res <- mse_simulation(
    om,
    scenarios = sc,
    n_sims = 1L,
    n_proj_years = 1L,
    min_assess_years = 99L,
    initial_tac = 2 * sum(om$true_params$B_initial),
    seed = 1
  )
  realised_catch <- sum(res$results$test$trajectories$catch[1, 1, ])
  expect_equal(realised_catch, 2 * sum(om$true_params$B_initial))
})

test_that("an explicit maximum exploitation rate constrains and records catch", {
  om <- make_simple_om()
  om$max_harvest_rate <- 0.8
  sc <- make_scenario()
  res <- mse_simulation(
    om,
    scenarios = sc,
    n_sims = 1L,
    n_proj_years = 1L,
    min_assess_years = 99L,
    initial_tac = 2 * sum(om$true_params$B_initial),
    seed = 1
  )

  trajectories <- res$results$test$trajectories
  expected_maximum <- 0.8 * sum(om$true_params$B_initial)
  expect_equal(sum(trajectories$requested_catch[1, 1, ]),
    2 * sum(om$true_params$B_initial))
  expect_equal(sum(trajectories$catch[1, 1, ]), expected_maximum)
  expect_true(all(trajectories$catch_constrained[1, 1, ]))
})

test_that("catch-share index weights spatial biomass rather than total biomass", {
  index_mean <- SurplusProductionModelMSE:::.catch_share_index(
    biomass = c(100, 10),
    catch = c(1, 9),
    q = c(0.01, 0.02)
  )
  # 0.1 x 0.01 x 100 + 0.9 x 0.02 x 10 = 0.28.
  expect_equal(index_mean, 0.28)
})

test_that(".make_obs_params creates valid obs_error_params", {
  tp <- list(sigma_obs = 0.2)
  oep <- SurplusProductionModelMSE:::.make_obs_params(tp)
  expect_s3_class(oep, "obs_error_params")
  expect_equal(oep$sigma, 0.2)
  expect_equal(oep$rho, 0)
})
