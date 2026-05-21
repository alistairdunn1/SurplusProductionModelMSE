# MSE Simulation Loop (Module 2.5)
#
# Core closed-loop simulation engine connecting the operating model,
# observation error, estimation model, harvest control rules, and
# implementation error into a full management strategy evaluation.


#' Run MSE Simulation
#'
#' Core closed-loop simulation engine that iterates over projection years:
#' operating model projection, observation error generation, estimation
#' model fitting (at assessment frequency), harvest control rule
#' application, and implementation error.
#'
#' @param operating_model An \code{\link{om_config}} object with
#'   \code{true_params} set.
#' @param estimation_model An \code{\link{em_config}} object, or
#'   \code{NULL} (default) for self-test mode where the EM inherits
#'   the OM structure.
#' @param scenarios A single \code{\link{create_scenario}} object or a
#'   list of them.
#' @param n_sims Integer. Number of simulation replicates (default 100).
#' @param n_proj_years Integer. Number of projection years (default 20).
#' @param obs_error_params An \code{obs_error_params} object from
#'   \code{\link{calibrate_observation_error}}, or \code{NULL} to
#'   construct from \code{operating_model$true_params$sigma_obs}
#'   with zero autocorrelation.
#' @param initial_tac Numeric scalar. TAC used before the first
#'   assessment. If \code{NULL} (default), the true MSY is used.
#' @param min_assess_years Integer. Minimum number of accumulated
#'   data years before the first EM fit (default 5).
#' @param parallel Logical. Use \pkg{future.apply} for parallel
#'   replicates (default \code{FALSE}). Requires the \code{future} and
#'   \code{future.apply} packages (listed in Suggests).
#' @param seed Integer random seed for reproducibility, or \code{NULL}.
#' @param verbose Logical. Print progress messages (default \code{FALSE}).
#'
#' @return An S3 object of class \code{mse_result} containing:
#'   \describe{
#'     \item{results}{Named list (one entry per scenario), each
#'       containing \code{trajectories} and \code{performance}.}
#'     \item{scenarios}{List of \code{mse_scenario} objects.}
#'     \item{om_config}{Operating model configuration.}
#'     \item{em_config}{Estimation model configuration (or \code{NULL}).}
#'     \item{n_sims}{Number of replicates.}
#'     \item{n_proj_years}{Number of projection years.}
#'     \item{call}{The matched call.}
#'   }
#'
#' Trajectory arrays have dimensions
#' \code{[n_sims x n_proj_years x n_areas]} for per-area metrics
#' (biomass, catch, harvest_rate) and
#' \code{[n_sims x n_proj_years]} for aggregate metrics
#' (tac, estimated_biomass).
#'
#' @details
#' The simulation loop for each replicate and year proceeds as:
#' \enumerate{
#'   \item Record true biomass.
#'   \item If an assessment year and enough data have accumulated,
#'     fit the estimation model to accumulated CPUE and catch,
#'     then apply the harvest control rule to set the TAC.
#'     Convergence failures are handled gracefully by retaining
#'     the previous assessment result.
#'   \item Apply implementation error to the TAC to obtain
#'     realized catch.
#'   \item Generate a CPUE observation and record catch for the
#'     next assessment.
#'   \item Project the operating model biomass forward one year.
#' }
#'
#' Between assessments the most recent TAC is carried forward.
#' Before the first assessment, \code{initial_tac} is used.
#'
#' @examples
#' \dontrun{
#' om <- om_config(
#'   n_areas = 1,
#'   true_params = list(
#'     r = 0.3, K = 5000, m = 2,
#'     sigma_obs = 0.2, q = 1e-4, B0 = 5000
#'   )
#' )
#' sc <- create_scenario("constant_f", hcr_constant_f(0.05))
#' res <- mse_simulation(om,
#'   scenarios = sc, n_sims = 10,
#'   n_proj_years = 15, seed = 42
#' )
#' res
#' }
#'
#' @export
mse_simulation <- function(operating_model,
                           estimation_model = NULL,
                           scenarios,
                           n_sims = 100L,
                           n_proj_years = 20L,
                           obs_error_params = NULL,
                           initial_tac = NULL,
                           min_assess_years = 5L,
                           parallel = FALSE,
                           seed = NULL,
                           verbose = FALSE) {
  cl <- match.call()

  # --- Input validation ---
  if (!inherits(operating_model, "om_config")) {
    stop("operating_model must be an 'om_config' object", call. = FALSE)
  }
  tp <- operating_model$true_params
  if (is.null(tp)) {
    stop("operating_model$true_params must be set", call. = FALSE)
  }
  if (!is.null(estimation_model) && !inherits(estimation_model, "em_config")) {
    stop("estimation_model must be an 'em_config' object or NULL",
      call. = FALSE
    )
  }

  # Accept single scenario or list
  if (inherits(scenarios, "mse_scenario")) {
    scenarios <- list(scenarios)
  }
  assert_list(scenarios, min.len = 1, .var.name = "scenarios")
  for (i in seq_along(scenarios)) {
    if (!inherits(scenarios[[i]], "mse_scenario")) {
      stop("scenarios[[", i, "]] must be an 'mse_scenario' object",
        call. = FALSE
      )
    }
  }

  assert_count(n_sims, positive = TRUE, .var.name = "n_sims")
  assert_count(n_proj_years, positive = TRUE, .var.name = "n_proj_years")
  assert_count(min_assess_years, positive = TRUE, .var.name = "min_assess_years")
  assert_flag(parallel, .var.name = "parallel")
  assert_flag(verbose, .var.name = "verbose")

  n_sims <- as.integer(n_sims)
  n_proj_years <- as.integer(n_proj_years)
  min_assess_years <- as.integer(min_assess_years)
  n_areas <- operating_model$n_areas

  # --- Construct obs_error_params if needed ---
  if (is.null(obs_error_params)) {
    obs_error_params <- .make_obs_params(tp)
  }

  # --- Default initial TAC from true MSY ---
  if (is.null(initial_tac)) {
    initial_tac <- .compute_true_msy(tp)
  }
  assert_number(initial_tac, lower = 0, .var.name = "initial_tac")

  # --- Pre-compute movement kernel ---
  movement_kernel <- NULL
  if (n_areas > 1 && operating_model$movement_rate > 0) {
    dm <- operating_model$distance_matrix
    if (is.null(dm)) {
      dm <- matrix(1, n_areas, n_areas)
      diag(dm) <- 0
    }
    movement_kernel <- build_movement_kernel(
      dm, operating_model$attractiveness, operating_model$decay
    )
  }

  # --- Generate per-sim seeds for reproducibility ---
  if (!is.null(seed)) {
    set.seed(seed)
    sim_seeds <- sample.int(
      .Machine$integer.max,
      n_sims * length(scenarios)
    )
  } else {
    sim_seeds <- NULL
  }

  # --- Run scenarios ---
  results <- list()

  for (s in seq_along(scenarios)) {
    scenario <- scenarios[[s]]
    if (verbose) message("Running scenario: ", scenario$name)

    seed_offset <- (s - 1L) * n_sims

    # Seeds for this scenario
    s_seeds <- if (!is.null(sim_seeds)) {
      sim_seeds[seed_offset + seq_len(n_sims)]
    } else {
      NULL
    }

    sim_fn <- function(i) {
      .run_single_sim(
        sim_seed         = if (!is.null(s_seeds)) s_seeds[i] else NULL,
        scenario         = scenario,
        om_config        = operating_model,
        em_config        = estimation_model,
        n_proj_years     = n_proj_years,
        n_areas          = n_areas,
        obs_error_params = obs_error_params,
        movement_kernel  = movement_kernel,
        initial_tac      = initial_tac,
        min_assess_years = min_assess_years,
        tp               = tp
      )
    }

    if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
      sim_results <- future.apply::future_lapply(
        seq_len(n_sims), sim_fn,
        future.seed = TRUE
      )
    } else {
      if (parallel) {
        message("future.apply not available; running sequentially")
      }
      sim_results <- lapply(seq_len(n_sims), sim_fn)
    }

    # --- Combine into trajectory arrays ---
    trajectories <- .combine_trajectories(
      sim_results, n_sims, n_proj_years, n_areas
    )

    # --- Performance metrics ---
    perf <- calculate_performance_metrics(trajectories, operating_model)

    results[[scenario$name]] <- list(
      trajectories = trajectories,
      performance  = perf
    )
  }

  structure(
    list(
      results      = results,
      scenarios    = scenarios,
      om_config    = operating_model,
      em_config    = estimation_model,
      n_sims       = n_sims,
      n_proj_years = n_proj_years,
      call         = cl
    ),
    class = "mse_result"
  )
}


# ========================================================================
# Internal helpers
# ========================================================================

#' Construct obs_error_params from true parameters
#' @noRd
.make_obs_params <- function(tp) {
  structure(
    list(
      sigma  = tp$sigma_obs,
      rho    = 0,
      labels = NULL,
      n_obs  = NA_integer_
    ),
    class = "obs_error_params"
  )
}


#' Compute true MSY from Pella-Tomlinson parameters
#' @noRd
.compute_true_msy <- function(tp) {
  r <- tp$r
  K <- tp$K
  m <- tp$m
  if (abs(m - 1) < 1e-10) {
    fmsy <- r / exp(1)
    bmsy <- K / exp(1)
  } else {
    bmsy <- K * (1 / m)^(1 / (m - 1))
    fmsy <- r * (1 - 1 / m) / m
  }
  fmsy * bmsy
}


#' Build reference_points list for HCR calls
#'
#' Translates EM output (or true parameter fallback) into the named list
#' expected by HCR closures: K, MSY, BMSY, FMSY.
#' @noRd
.build_hcr_ref_points <- function(em_result, tp) {
  if (!is.null(em_result)) {
    list(
      K    = em_result$K,
      MSY  = em_result$msy,
      BMSY = em_result$bmsy,
      FMSY = em_result$fmsy
    )
  } else {
    # Fallback to true parameters
    r <- tp$r
    K <- tp$K
    m <- tp$m
    if (abs(m - 1) < 1e-10) {
      fmsy <- r / exp(1)
      bmsy <- K / exp(1)
    } else {
      bmsy <- K * (1 / m)^(1 / (m - 1))
      fmsy <- r * (1 - 1 / m) / m
    }
    list(K = K, MSY = fmsy * bmsy, BMSY = bmsy, FMSY = fmsy)
  }
}


#' Safely fit the estimation model
#'
#' Wraps \code{fit_pella_tomlinson_model} in tryCatch so convergence
#' failures return \code{NULL} instead of propagating errors.
#' @noRd
.fit_em_safely <- function(cpue_history, catch_history, em_config,
                           om_config) {
  data <- list(
    cpue_data  = cpue_history,
    catch_data = catch_history
  )

  tryCatch(
    withCallingHandlers(
      {
        fit <- fit_pella_tomlinson_model(
          data,
          options = list(silent = TRUE, validate_data = FALSE)
        )

        if (!fit$fitted) {
          return(NULL)
        }

        ref <- calculate_reference_points(fit)
        bio <- estimate_biomass(fit)
        current_b <- tail(bio$biomass, 1)

        if (length(current_b) != 1 || !is.finite(current_b)) {
          return(NULL)
        }

        list(
          est_biomass = current_b,
          K           = fit$parameters[["K"]],
          msy         = ref$msy,
          bmsy        = ref$bmsy,
          fmsy        = ref$fmsy,
          converged   = TRUE
        )
      },
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
}


#' Run a single MSE simulation replicate
#' @noRd
.run_single_sim <- function(sim_seed, scenario, om_config, em_config,
                            n_proj_years, n_areas, obs_error_params,
                            movement_kernel, initial_tac, min_assess_years,
                            tp) {
  if (!is.null(sim_seed)) set.seed(sim_seed)

  B0 <- rep_len(tp$B0, n_areas)

  # State
  biomass <- B0
  impl_eps <- rep(0, n_areas)

  # Storage
  biomass_store <- matrix(NA_real_, n_proj_years, n_areas)
  catch_store <- matrix(NA_real_, n_proj_years, n_areas)
  tac_store <- rep(NA_real_, n_proj_years)
  est_bio_store <- rep(NA_real_, n_proj_years)

  # Accumulated histories for EM
  cpue_history <- data.frame(year = integer(0), cpue = numeric(0))
  catch_history <- data.frame(year = integer(0), catch = numeric(0))

  # Management state
  current_tac <- initial_tac
  last_em_result <- NULL
  assess_freq <- scenario$assessment_frequency

  # q for aggregate CPUE simulation
  q_agg <- if (length(tp$q) > 1) mean(tp$q) else tp$q

  for (yr in seq_len(n_proj_years)) {
    # Record true biomass
    biomass_store[yr, ] <- biomass

    # 1. Assessment decision (uses data accumulated from years 1..yr-1)
    enough_data <- nrow(cpue_history) >= min_assess_years
    on_schedule <- (yr %% assess_freq) == 0
    if (on_schedule && enough_data) {
      em_result <- .fit_em_safely(
        cpue_history, catch_history, em_config, om_config
      )

      if (!is.null(em_result)) {
        last_em_result <- em_result
      }
      # else: keep previous result (convergence failure handling)

      # Estimated biomass for HCR
      est_b <- if (!is.null(last_em_result)) {
        last_em_result$est_biomass
      } else {
        sum(biomass) # fallback: perfect knowledge
      }

      ref_pts <- .build_hcr_ref_points(last_em_result, tp)
      current_tac <- scenario$harvest_control_rule(est_b, ref_pts)
    }

    tac_store[yr] <- current_tac
    eb <- if (!is.null(last_em_result)) last_em_result$est_biomass else NULL
    est_bio_store[yr] <- if (length(eb) == 1) eb else NA_real_

    # 2. Implementation error → realized catch
    if (n_areas > 1) {
      tac_area <- current_tac * (biomass / pmax(sum(biomass), 1e-8))
    } else {
      tac_area <- current_tac
    }

    impl_result <- apply_implementation_error(
      tac_area, scenario$implementation_error, impl_eps
    )
    catch <- impl_result$catch
    impl_eps <- impl_result$eps
    catch_store[yr, ] <- catch

    # 3. Generate CPUE observation from total biomass (pre-fishing)
    cpue_yr <- simulate_cpue(
      true_biomass = sum(biomass),
      obs_error_params = obs_error_params,
      q = q_agg,
      years = as.integer(yr)
    )
    cpue_history <- rbind(cpue_history, cpue_yr[, c("year", "cpue")])
    catch_history <- rbind(
      catch_history,
      data.frame(year = as.integer(yr), catch = sum(catch))
    )

    # 4. Project biomass forward
    biomass <- project_biomass(
      biomass, catch, om_config, movement_kernel
    )
  }

  # Derived: harvest rate
  harvest_rate_store <- catch_store / pmax(biomass_store, 1e-8)

  list(
    biomass           = biomass_store,
    catch             = catch_store,
    tac               = tac_store,
    estimated_biomass = est_bio_store,
    harvest_rate      = harvest_rate_store
  )
}


#' Combine per-replicate results into trajectory arrays
#' @noRd
.combine_trajectories <- function(sim_results, n_sims, n_proj_years,
                                  n_areas) {
  biomass <- array(NA_real_, dim = c(n_sims, n_proj_years, n_areas))
  catch_arr <- array(NA_real_, dim = c(n_sims, n_proj_years, n_areas))
  hr_arr <- array(NA_real_, dim = c(n_sims, n_proj_years, n_areas))
  tac_mat <- matrix(NA_real_, n_sims, n_proj_years)
  est_bio <- matrix(NA_real_, n_sims, n_proj_years)

  for (i in seq_len(n_sims)) {
    sim <- sim_results[[i]]
    biomass[i, , ] <- sim$biomass
    catch_arr[i, , ] <- sim$catch
    hr_arr[i, , ] <- sim$harvest_rate
    tac_mat[i, ] <- sim$tac
    est_bio[i, ] <- sim$estimated_biomass
  }

  list(
    biomass           = biomass,
    catch             = catch_arr,
    harvest_rate      = hr_arr,
    tac               = tac_mat,
    estimated_biomass = est_bio
  )
}


# ========================================================================
# S3 methods
# ========================================================================

#' @export
print.mse_result <- function(x, ...) {
  cat("MSE Simulation Result\n")
  cat("=====================\n")
  cat("Scenarios:       ", length(x$scenarios), "\n")
  cat("Simulations:     ", x$n_sims, "\n")
  cat("Projection years:", x$n_proj_years, "\n")
  cat("Areas:           ", x$om_config$n_areas, "\n\n")

  cat("Scenarios:\n")
  for (nm in names(x$results)) {
    cat("  -", nm, "\n")
  }
  invisible(x)
}
