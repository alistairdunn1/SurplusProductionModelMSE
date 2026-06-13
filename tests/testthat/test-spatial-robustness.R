# Tests for Spatial Robustness Testing (Module 2.4)

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

run_quick_mse <- function(om, scenarios, n_sims = 5L, n_proj = 8L,
                          seed = 1, ...) {
  mse_simulation(initial_tac = 0, 
    operating_model = om,
    scenarios = scenarios,
    n_sims = n_sims,
    n_proj_years = n_proj,
    min_assess_years = 99L,
    seed = seed,
    ...
  )
}

# ========================================================================
# compare_scenarios: input validation
# ========================================================================

test_that("compare_scenarios rejects non-mse_result input", {
  expect_error(
    compare_scenarios("bad", "mean_catch", "AAV"),
    "mse_result"
  )
})

test_that("compare_scenarios rejects unknown metric_x", {
  om <- make_simple_om()
  sc <- create_scenario("s1", hcr_constant_f(0.05))
  res <- run_quick_mse(om, sc)

  expect_error(
    compare_scenarios(res, "nonexistent_metric", "mean_catch"),
    "not found"
  )
})

test_that("compare_scenarios rejects unknown metric_y", {
  om <- make_simple_om()
  sc <- create_scenario("s1", hcr_constant_f(0.05))
  res <- run_quick_mse(om, sc)

  expect_error(
    compare_scenarios(res, "mean_catch", "nonexistent_metric"),
    "not found"
  )
})

# ========================================================================
# compare_scenarios: basic functionality
# ========================================================================

test_that("compare_scenarios returns a ggplot with single mse_result", {
  om <- make_simple_om()
  sc1 <- create_scenario("low_f", hcr_constant_f(0.03))
  sc2 <- create_scenario("high_f", hcr_constant_f(0.10))
  res <- run_quick_mse(om, list(sc1, sc2))

  p <- compare_scenarios(res, "mean_catch", "Pr(B<20%K)_final")
  expect_s3_class(p, "ggplot")
})

test_that("compare_scenarios returns ggplot with list of mse_results", {
  om <- make_simple_om()
  sc1 <- create_scenario("low", hcr_constant_f(0.03))
  sc2 <- create_scenario("high", hcr_constant_f(0.10))
  res1 <- run_quick_mse(om, sc1, seed = 1)
  res2 <- run_quick_mse(om, sc2, seed = 2)

  p <- compare_scenarios(list(res1, res2), "mean_catch", "AAV")
  expect_s3_class(p, "ggplot")
})

test_that("compare_scenarios plot contains correct number of points", {
  om <- make_simple_om()
  sc1 <- create_scenario("s1", hcr_constant_f(0.03))
  sc2 <- create_scenario("s2", hcr_constant_f(0.06))
  sc3 <- create_scenario("s3", hcr_constant_f(0.10))
  res <- run_quick_mse(om, list(sc1, sc2, sc3))

  p <- compare_scenarios(res, "mean_catch", "Pr(B<50%K)_final")
  plot_data <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(nrow(plot_data), 3L)
})

test_that("compare_scenarios works with biomass ratio metrics", {
  om <- make_simple_om()
  sc1 <- create_scenario("s1", hcr_constant_f(0.03))
  sc2 <- create_scenario("s2", hcr_constant_f(0.10))
  res <- run_quick_mse(om, list(sc1, sc2))

  p <- compare_scenarios(res, "mean_B_BMSY", "mean_catch")
  expect_s3_class(p, "ggplot")
})

# ========================================================================
# Pareto optimality
# ========================================================================

test_that("Pareto frontier is correct for simple dominance", {
  # A dominates B in mean_catch (higher better) and risk (lower better)
  # So A is Pareto optimal and B is not
  x <- c(200, 100) # mean_catch: higher better → dir = +1
  y <- c(0.1, 0.5) # Pr(...): lower better → dir = -1
  result <- SurplusProductionModelMSE:::.is_pareto_optimal(x, y, "mean_catch", "Pr(B<20%K)_final")
  expect_equal(result, c(TRUE, FALSE))
})

test_that("Pareto frontier with trade-off keeps both points", {
  # Neither dominates the other
  x <- c(200, 100) # A has higher catch
  y <- c(0.5, 0.1) # B has lower risk
  result <- SurplusProductionModelMSE:::.is_pareto_optimal(x, y, "mean_catch", "Pr(B<20%K)_final")
  expect_equal(result, c(TRUE, TRUE))
})

test_that("single scenario is always Pareto optimal", {
  result <- SurplusProductionModelMSE:::.is_pareto_optimal(100, 0.3, "mean_catch", "AAV")
  expect_true(result)
})

test_that("metric direction: risk metrics are lower-is-better", {
  expect_equal(SurplusProductionModelMSE:::.metric_direction("Pr(B<50%K)_final"), -1)
  expect_equal(SurplusProductionModelMSE:::.metric_direction("Pr(F>F50%K)_ever"), -1)
  expect_equal(SurplusProductionModelMSE:::.metric_direction("AAV"), -1)
})

test_that("metric direction: catch and biomass ratios are higher-is-better", {
  expect_equal(SurplusProductionModelMSE:::.metric_direction("mean_catch"), 1)
  expect_equal(SurplusProductionModelMSE:::.metric_direction("mean_B_BMSY"), 1)
  expect_equal(SurplusProductionModelMSE:::.metric_direction("final_B_K"), 1)
})

# ========================================================================
# run_self_test: input validation
# ========================================================================

test_that("run_self_test rejects non-om_config", {
  expect_error(
    run_self_test(initial_tac = 0, "bad", hcr_constant_f(0.05)),
    "om_config"
  )
})

test_that("run_self_test rejects non-function HCR", {
  om <- make_simple_om()
  expect_error(
    run_self_test(initial_tac = 0, om, "not_a_function"),
    "harvest_control_rule"
  )
})

# ========================================================================
# run_self_test: single area
# ========================================================================

test_that("self-test completes for single-area OM", {
  om <- make_simple_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.05),
    n_sims = 5L, n_proj_years = 8L, seed = 1,
    min_assess_years = 99L
  )
  expect_s3_class(st, "self_test_result")
})

test_that("self-test reports biomass conservation", {
  om <- make_simple_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.05),
    n_sims = 5L, n_proj_years = 8L, seed = 1,
    min_assess_years = 99L
  )
  expect_true(st$biomass_conserved)
})

test_that("self-test reports positive biomass", {
  om <- make_simple_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.05),
    n_sims = 5L, n_proj_years = 8L, seed = 1,
    min_assess_years = 99L
  )
  expect_true(st$positive_biomass)
})

test_that("self-test summary is a data frame", {
  om <- make_simple_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.05),
    n_sims = 5L, n_proj_years = 8L, seed = 1,
    min_assess_years = 99L
  )
  expect_s3_class(st$summary, "data.frame")
  expect_true("metric" %in% names(st$summary))
  expect_true("value" %in% names(st$summary))
})

test_that("self-test contains mse_result", {
  om <- make_simple_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.05),
    n_sims = 3L, n_proj_years = 5L, seed = 1,
    min_assess_years = 99L
  )
  expect_s3_class(st$mse_result, "mse_result")
})

# ========================================================================
# run_self_test: multi-area (spatial)
# ========================================================================

test_that("self-test completes for multi-area OM", {
  om <- make_spatial_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.05),
    n_sims = 3L, n_proj_years = 6L, seed = 1,
    min_assess_years = 99L
  )
  expect_s3_class(st, "self_test_result")
})

test_that("self-test biomass conservation for multi-area", {
  om <- make_spatial_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.05),
    n_sims = 5L, n_proj_years = 6L, seed = 1,
    min_assess_years = 99L
  )
  # B_initial = c(3000, 2000), total = 5000
  expect_true(st$biomass_conserved)
})

test_that("OM conserves total biomass under movement", {
  # With no catch (or minimal), total biomass should be conserved
  # by the movement kernel
  om <- make_spatial_om()
  B_initial <- om$true_params$B_initial
  dm <- om$movement_cost_matrix
  kernel <- build_movement_kernel(dm)

  # Movement redistributes but conserves total
  B_new <- project_biomass(B_initial, c(0, 0), om, kernel,
    process_noise = FALSE
  )
  expect_equal(sum(B_new), sum(B_initial), tolerance = 1e-6)
})

# ========================================================================
# Self-test depletion sanity
# ========================================================================

test_that("self-test with low F maintains healthy stock", {
  om <- make_simple_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.02),
    n_sims = 10L, n_proj_years = 10L, seed = 42,
    min_assess_years = 99L
  )
  s <- st$summary
  depletion <- s$value[s$metric == "mean_B_K" & s$scope == "aggregate"]
  # Low F should keep stock well above 50% K

  expect_true(depletion > 0.5)
})

test_that("self-test with high F depletes stock", {
  om <- make_simple_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.20),
    n_sims = 10L, n_proj_years = 15L, seed = 42,
    min_assess_years = 5L
  )
  s <- st$summary
  depletion <- s$value[s$metric == "final_B_K" & s$scope == "aggregate"]
  # High F should deplete stock below K once assessments begin applying the HCR
  expect_true(depletion < 0.9)
})

# ========================================================================
# Print method
# ========================================================================

test_that("print.self_test_result works", {
  om <- make_simple_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.05),
    n_sims = 3L, n_proj_years = 5L, seed = 1,
    min_assess_years = 99L
  )
  out <- capture.output(print(st))
  expect_true(any(grepl("Self-Test", out)))
  expect_true(any(grepl("conserved", out, ignore.case = TRUE)))
})

# ========================================================================
# compare_scenarios integrated with self-test
# ========================================================================

test_that("compare_scenarios works with self-test result", {
  om <- make_simple_om()
  st <- run_self_test(initial_tac = 0, om, hcr_constant_f(0.05),
    n_sims = 3L, n_proj_years = 5L, seed = 1,
    min_assess_years = 99L
  )
  # self_test_result$mse_result is an mse_result
  p <- compare_scenarios(st$mse_result, "mean_catch", "AAV")
  expect_s3_class(p, "ggplot")
})

# ========================================================================
# .extract_scenario_metrics
# ========================================================================

test_that(".extract_scenario_metrics returns correct columns", {
  om <- make_simple_om()
  sc1 <- create_scenario("s1", hcr_constant_f(0.03))
  sc2 <- create_scenario("s2", hcr_constant_f(0.10))
  res <- run_quick_mse(om, list(sc1, sc2))

  df <- SurplusProductionModelMSE:::.extract_scenario_metrics(res, "aggregate")
  expect_true(all(c("metric", "scope", "value", "scenario") %in% names(df)))
  expect_true(all(c("s1", "s2") %in% df$scenario))
})

test_that(".extract_scenario_metrics pools across mse_result list", {
  om <- make_simple_om()
  sc1 <- create_scenario("s1", hcr_constant_f(0.03))
  sc2 <- create_scenario("s2", hcr_constant_f(0.10))
  res1 <- run_quick_mse(om, sc1, seed = 1)
  res2 <- run_quick_mse(om, sc2, seed = 2)

  df <- SurplusProductionModelMSE:::.extract_scenario_metrics(list(res1, res2), "aggregate")
  expect_true(all(c("s1", "s2") %in% df$scenario))
})

# ========================================================================
# OM/EM configuration compatibility
# ========================================================================

test_that("OM/EM configuration: self-test (NULL em) works", {
  om <- make_simple_om()
  # estimation_model = NULL should be accepted
  expect_no_error(
    mse_simulation(initial_tac = 0, 
      om,
      estimation_model = NULL,
      scenarios = create_scenario("t", hcr_constant_f(0.05)),
      n_sims = 2L, n_proj_years = 3L, min_assess_years = 99L, seed = 1
    )
  )
})

test_that("OM/EM configuration: explicit em_config works", {
  om <- make_simple_om()
  em <- em_config(n_areas = 1L)
  expect_no_error(
    mse_simulation(initial_tac = 0, 
      om,
      estimation_model = em,
      scenarios = create_scenario("t", hcr_constant_f(0.05)),
      n_sims = 2L, n_proj_years = 3L, min_assess_years = 99L, seed = 1
    )
  )
})

# ========================================================================
# Relative performance: higher F → more catch, more risk
# ========================================================================

test_that("higher catch produces higher mean catch metric", {
  om <- make_simple_om()
  # Use constant_catch with different initial_tac to guarantee different catch
  sc_lo <- create_scenario("low", hcr_constant_catch(100))
  sc_hi <- create_scenario("high", hcr_constant_catch(300))
  # Provide initial_tac matching each scenario — run separately
  res_lo <- mse_simulation(om,
    scenarios = sc_lo, n_sims = 5L,
    n_proj_years = 10L, min_assess_years = 99L,
    initial_tac = 100, seed = 1
  )
  res_hi <- mse_simulation(om,
    scenarios = sc_hi, n_sims = 5L,
    n_proj_years = 10L, min_assess_years = 99L,
    initial_tac = 300, seed = 1
  )

  df <- SurplusProductionModelMSE:::.extract_scenario_metrics(list(res_lo, res_hi), "aggregate")
  catch_lo <- df$value[df$scenario == "low" & df$metric == "mean_catch"]
  catch_hi <- df$value[df$scenario == "high" & df$metric == "mean_catch"]
  expect_true(catch_hi > catch_lo)
})

test_that("higher catch produces lower final depletion", {
  om <- make_simple_om()
  sc_lo <- create_scenario("low", hcr_constant_catch(100))
  sc_hi <- create_scenario("high", hcr_constant_catch(300))
  res_lo <- mse_simulation(om,
    scenarios = sc_lo, n_sims = 5L,
    n_proj_years = 10L, min_assess_years = 99L,
    initial_tac = 100, seed = 1
  )
  res_hi <- mse_simulation(om,
    scenarios = sc_hi, n_sims = 5L,
    n_proj_years = 10L, min_assess_years = 99L,
    initial_tac = 300, seed = 1
  )

  df <- SurplusProductionModelMSE:::.extract_scenario_metrics(list(res_lo, res_hi), "aggregate")
  dep_lo <- df$value[df$scenario == "low" & df$metric == "final_B_K"]
  dep_hi <- df$value[df$scenario == "high" & df$metric == "final_B_K"]
  expect_true(dep_lo > dep_hi)
})
