# Tests for Harvest Control Rules
# Tests the four HCR factory functions:
#   hcr_constant_f, hcr_constant_catch, hcr_hockey_stick, hcr_ccamlr_krill

# ── hcr_constant_f ─────────────────────────────────────────────────────────────

test_that("hcr_constant_f returns a function", {
  hcr <- hcr_constant_f(0.1)
  expect_true(is.function(hcr))
  expect_equal(length(formals(hcr)), 2)
})

test_that("hcr_constant_f computes correct TAC", {
  hcr <- hcr_constant_f(0.1)
  expect_equal(hcr(1000, list(BMSY = 500)), 100)
  expect_equal(hcr(5000, list()), 500)
  expect_equal(hcr(0, list()), 0)
})

test_that("hcr_constant_f TAC scales linearly with biomass", {
  hcr <- hcr_constant_f(0.2)
  tac_low <- hcr(100, list())
  tac_high <- hcr(1000, list())
  expect_equal(tac_high / tac_low, 10)
})

test_that("hcr_constant_f rejects invalid f_target", {
  expect_error(hcr_constant_f(0))
  expect_error(hcr_constant_f(-0.1))
  expect_error(hcr_constant_f("abc"))
})

test_that("hcr_constant_f works with create_scenario", {
  hcr <- hcr_constant_f(0.1)
  sc <- create_scenario("ConstF", hcr)
  expect_s3_class(sc, "mse_scenario")
})


# ── hcr_constant_catch ─────────────────────────────────────────────────────────

test_that("hcr_constant_catch returns a function", {
  hcr <- hcr_constant_catch(200)
  expect_true(is.function(hcr))
})

test_that("hcr_constant_catch returns fixed catch when biomass is sufficient", {
  hcr <- hcr_constant_catch(200)
  expect_equal(hcr(1000, list()), 200)
  expect_equal(hcr(500, list()), 200)
  expect_equal(hcr(200, list()), 200)
})

test_that("hcr_constant_catch clamps to biomass when biomass is low", {
  hcr <- hcr_constant_catch(200)
  expect_equal(hcr(50, list()), 50)
  expect_equal(hcr(0, list()), 0)
})

test_that("hcr_constant_catch rejects invalid catch_target", {
  expect_error(hcr_constant_catch(0))
  expect_error(hcr_constant_catch(-100))
})


# ── hcr_hockey_stick ───────────────────────────────────────────────────────────

test_that("hcr_hockey_stick returns a function", {
  hcr <- hcr_hockey_stick(0.1, 0.2, 0.4)
  expect_true(is.function(hcr))
})

test_that("hcr_hockey_stick above target gives full F", {
  hcr <- hcr_hockey_stick(0.1, 0.2, 0.4)
  ref <- list(K = 5000)

  # B/K = 3000/5000 = 0.6 > 0.4
  expect_equal(hcr(3000, ref), 0.1 * 3000)

  # Exactly at target: B/K = 2000/5000 = 0.4
  expect_equal(hcr(2000, ref), 0.1 * 2000)
})

test_that("hcr_hockey_stick in ramp zone gives reduced F", {
  hcr <- hcr_hockey_stick(0.1, 0.2, 0.4)
  ref <- list(K = 5000)

  # B/K = 1500/5000 = 0.3; multiplier = (0.3 - 0.2)/(0.4 - 0.2) = 0.5
  expected <- 0.1 * 1500 * 0.5
  expect_equal(hcr(1500, ref), expected)
})

test_that("hcr_hockey_stick ramp is linear", {
  hcr <- hcr_hockey_stick(0.1, 0.2, 0.4)
  ref <- list(K = 10000)

  # Test two points in the ramp
  # B = 2500, B/K = 0.25, mult = (0.25-0.2)/(0.4-0.2) = 0.25
  tac1 <- hcr(2500, ref)
  expected1 <- 0.1 * 2500 * 0.25
  expect_equal(tac1, expected1)

  # B = 3500, B/K = 0.35, mult = (0.35-0.2)/(0.4-0.2) = 0.75
  tac2 <- hcr(3500, ref)
  expected2 <- 0.1 * 3500 * 0.75
  expect_equal(tac2, expected2)
})

test_that("hcr_hockey_stick below limit gives zero catch", {
  hcr <- hcr_hockey_stick(0.1, 0.2, 0.4)
  ref <- list(K = 5000)

  # B/K = 500/5000 = 0.1 < 0.2
  expect_equal(hcr(500, ref), 0)

  # Exactly at limit: B/K = 1000/5000 = 0.2
  expect_equal(hcr(1000, ref), 0)

  # Zero biomass
  expect_equal(hcr(0, ref), 0)
})

test_that("hcr_hockey_stick requires K in reference_points", {
  hcr <- hcr_hockey_stick(0.1, 0.2, 0.4)
  expect_error(hcr(1000, list()), "K")
  expect_error(hcr(1000, list(K = NULL)), "K")
})

test_that("hcr_hockey_stick rejects b_limit >= b_target", {
  expect_error(hcr_hockey_stick(0.1, 0.4, 0.4), "b_limit")
  expect_error(hcr_hockey_stick(0.1, 0.5, 0.4), "b_limit")
})

test_that("hcr_hockey_stick rejects invalid parameters", {
  expect_error(hcr_hockey_stick(0, 0.2, 0.4))
  expect_error(hcr_hockey_stick(0.1, 0, 0.4))
  expect_error(hcr_hockey_stick(0.1, 0.2, 0))
  expect_error(hcr_hockey_stick(0.1, 1.5, 0.4))
})


# ── hcr_ccamlr_krill ──────────────────────────────────────────────────────────

test_that("hcr_ccamlr_krill returns a function", {
  hcr <- hcr_ccamlr_krill()
  expect_true(is.function(hcr))
})

test_that("hcr_ccamlr_krill above target gives gamma * MSY", {
  hcr <- hcr_ccamlr_krill(gamma = 0.5)
  ref <- list(K = 5000, MSY = 400)

  # B/K = 3000/5000 = 0.6 > 0.5
  expect_equal(hcr(3000, ref), 0.5 * 400)

  # Exactly at target: B/K = 2500/5000 = 0.5
  expect_equal(hcr(2500, ref), 0.5 * 400)
})

test_that("hcr_ccamlr_krill in ramp zone gives reduced TAC", {
  hcr <- hcr_ccamlr_krill(
    gamma = 0.5, target_depletion = 0.5,
    limit_depletion = 0.2
  )
  ref <- list(K = 5000, MSY = 400)

  # B/K = 1750/5000 = 0.35; mult = (0.35 - 0.2)/(0.5 - 0.2) = 0.5
  expected <- 0.5 * 400 * 0.5
  expect_equal(hcr(1750, ref), expected)
})

test_that("hcr_ccamlr_krill below limit gives zero", {
  hcr <- hcr_ccamlr_krill()
  ref <- list(K = 5000, MSY = 400)

  # B/K = 500/5000 = 0.1 < 0.2
  expect_equal(hcr(500, ref), 0)

  # Exactly at limit: B/K = 1000/5000 = 0.2
  expect_equal(hcr(1000, ref), 0)
})

test_that("hcr_ccamlr_krill requires K and MSY", {
  hcr <- hcr_ccamlr_krill()
  expect_error(hcr(1000, list(MSY = 400)), "K")
  expect_error(hcr(1000, list(K = 5000)), "MSY")
})

test_that("hcr_ccamlr_krill rejects limit >= target", {
  expect_error(hcr_ccamlr_krill(target_depletion = 0.3, limit_depletion = 0.3))
  expect_error(hcr_ccamlr_krill(target_depletion = 0.3, limit_depletion = 0.5))
})

test_that("hcr_ccamlr_krill respects custom gamma", {
  hcr <- hcr_ccamlr_krill(gamma = 0.75)
  ref <- list(K = 5000, MSY = 400)
  # Above target
  expect_equal(hcr(3000, ref), 0.75 * 400)
})

test_that("hcr_ccamlr_krill rejects invalid gamma", {
  expect_error(hcr_ccamlr_krill(gamma = 0))
  expect_error(hcr_ccamlr_krill(gamma = -1))
})


# ── cross-cutting ──────────────────────────────────────────────────────────────

test_that("all HCRs return non-negative TAC", {
  hcr_f <- hcr_constant_f(0.1)
  hcr_c <- hcr_constant_catch(200)
  hcr_h <- hcr_hockey_stick(0.1, 0.2, 0.4)
  hcr_k <- hcr_ccamlr_krill()
  ref <- list(K = 5000, MSY = 400, BMSY = 2500)

  for (b in c(0, 100, 500, 2500, 5000)) {
    expect_gte(hcr_f(b, ref), 0)
    expect_gte(hcr_c(b, ref), 0)
    expect_gte(hcr_h(b, ref), 0)
    expect_gte(hcr_k(b, ref), 0)
  }
})

test_that("all HCRs compatible with create_scenario", {
  ref <- list(K = 5000, MSY = 400)
  hcrs <- list(
    hcr_constant_f(0.1),
    hcr_constant_catch(200),
    hcr_hockey_stick(0.1, 0.2, 0.4),
    hcr_ccamlr_krill()
  )
  for (hcr in hcrs) {
    sc <- create_scenario("test", hcr)
    expect_s3_class(sc, "mse_scenario")
  }
})
