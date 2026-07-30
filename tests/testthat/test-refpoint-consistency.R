# Cross-package consistency: the OM truth, the performance-metric reference
# points, and the HCR reference points must all agree with the assessment
# package's canonical Pella-Tomlinson reference points and dynamics.

test_that("performance-metric reference points match the canonical helper", {
  for (m in c(1.5, 2, 3)) {
    r <- 0.3
    K <- 5000
    cfg <- om_config(
      n_areas = 1,
      true_params = list(
        r = r, K = K, m = m, sigma_obs = 0.2, q = 1e-4, B_initial = K
      )
    )
    traj <- list(
      biomass = matrix(0.5 * K, 4, 4),
      catch = matrix(10, 4, 4)
    )
    perf <- calculate_performance_metrics(traj, cfg)
    rp <- SurplusProductionModel::pella_tomlinson_reference_points(r, K, m)

    expect_equal(perf$reference$BMSY, rp$bmsy, tolerance = 1e-10)
    expect_equal(perf$reference$FMSY, rp$fmsy, tolerance = 1e-10)
    expect_equal(perf$reference$MSY, rp$msy, tolerance = 1e-10)
  }
})

test_that("HCR reference points require a successful estimation model", {
  tp <- list(r = 0.3, K = 5000, m = 2)
  expect_error(
    SurplusProductionModelMSE:::.build_hcr_ref_points(NULL, tp),
    "successful estimation-model result"
  )
})

test_that("OM production at BMSY equals MSY (dynamics match reference points)", {
  r <- 0.3
  K <- 5000
  m <- 2
  cfg <- om_config(
    n_areas = 1,
    true_params = list(
      r = r, K = K, m = m, sigma_obs = 0.2, q = 1e-4, B_initial = K
    )
  )
  rp <- SurplusProductionModel::pella_tomlinson_reference_points(r, K, m)
  b_next <- project_biomass(rp$bmsy, 0, cfg, process_noise = FALSE)
  expect_equal(b_next - rp$bmsy, rp$msy, tolerance = 1e-6)
})

test_that("OM movement conserves biomass with asymmetric attractiveness", {
  dm <- matrix(c(0, 100, 200, 100, 0, 100, 200, 100, 0), nrow = 3)
  cfg_move <- om_config(
    n_areas = 3,
    movement_rate = 0.4,
    movement_cost_matrix = dm,
    attractiveness = c(1, 2, 0.5),
    decay = 0.02,
    true_params = list(
      r = 0.3, K = 9000, m = 2, sigma_obs = 0.2,
      q = rep(1e-4, 3), B_initial = c(3000, 3000, 3000)
    )
  )
  cfg_nomove <- cfg_move
  cfg_nomove$movement_rate <- 0

  B <- c(4000, 3000, 2000)
  b_move <- project_biomass(B, c(0, 0, 0), cfg_move, process_noise = FALSE)
  b_nomove <- project_biomass(B, c(0, 0, 0), cfg_nomove, process_noise = FALSE)

  # Movement only redistributes: totals must match the no-movement total.
  expect_equal(sum(b_move), sum(b_nomove), tolerance = 1e-8)
})
