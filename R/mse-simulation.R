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
#'   assessment. This must be supplied explicitly by the caller.
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
#'     sigma_obs = 0.2, q = 1e-4, B_initial = 5000
#'   )
#' )
#' sc <- create_scenario("constant_f", hcr_constant_f(0.05))
#' res <- mse_simulation(om,
#'   scenarios = sc, n_sims = 10,
#'   n_proj_years = 15, initial_tac = 0, seed = 42
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
                           initial_tac,
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

  if (missing(initial_tac) || is.null(initial_tac)) {
    stop(
      "initial_tac must be supplied explicitly; no default is applied",
      call. = FALSE
    )
  }
  assert_number(initial_tac, lower = 0, .var.name = "initial_tac")

  # --- Pre-compute movement kernel ---
  movement_kernel <- NULL
  if (n_areas > 1 && operating_model$movement_rate > 0) {
    dm <- operating_model$movement_cost_matrix
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
    if (verbose) {
      cat(
        "Running scenario ", s, " of ", length(scenarios),
        ": ", scenario$name, "\n",
        sep = ""
      )
      flush.console()
    }

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
      if (!requireNamespace("future", quietly = TRUE)) {
        stop(
          "parallel=TRUE requires package 'future', but it is not installed."
        )
      }
      active_workers <- future::nbrOfWorkers()
      if (!is.finite(active_workers) || active_workers <= 1L) {
        stop(
          "parallel=TRUE requires an active future plan with >1 worker. ",
          "Current plan reports ", active_workers, " worker(s)."
        )
      }
      sim_results <- future.apply::future_lapply(
        seq_len(n_sims), sim_fn,
        future.seed = TRUE
      )
    } else {
      if (parallel) {
        stop(
          "parallel=TRUE requires package 'future.apply', but it is not installed."
        )
      }
      if (verbose) {
        sim_results <- vector("list", n_sims)
        progress_every <- max(1L, floor(n_sims / 20L))
        for (i in seq_len(n_sims)) {
          if (i == 1L || i == n_sims || (i %% progress_every) == 0L) {
            pct_done <- round(100 * i / n_sims, 1)
            cat(
              "  Simulation ", i, " of ", n_sims,
              " (", pct_done, "%)",
              " [", scenario$name, "]\n",
              sep = ""
            )
            flush.console()
          }
          sim_results[[i]] <- sim_fn(i)
        }
      } else {
        sim_results <- lapply(seq_len(n_sims), sim_fn)
      }
    }

    # --- Combine into trajectory arrays ---
    trajectories <- .combine_trajectories(
      sim_results, n_sims, n_proj_years, n_areas
    )

    first_assessment_year <- as.integer(min_assess_years + 1L)
    while ((first_assessment_year %% scenario$assessment_frequency) != 0L) {
      first_assessment_year <- first_assessment_year + 1L
    }
    # When the first assessment falls beyond the projection horizon (short
    # runs, or scenarios with no in-horizon assessment) clamp the metric
    # window so performance is still summarised rather than erroring.
    metric_start_year <- max(1L, min(first_assessment_year, n_proj_years))

    # --- Performance metrics ---
    perf <- calculate_performance_metrics(
      trajectories,
      operating_model,
      start_year = metric_start_year
    )

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
    # Fallback to true parameters, using the canonical Pella-Tomlinson
    # reference points shared with the assessment package.
    rp <- SurplusProductionModel::pella_tomlinson_reference_points(
      r = tp$r, K = tp$K, m = tp$m
    )
    list(K = tp$K, MSY = rp$msy, BMSY = rp$bmsy, FMSY = rp$fmsy)
  }
}


#' Normalise EM fixed parameters for SurplusProductionModel fitting
#'
#' Converts natural-scale parameter names (for example, m = 2) to the
#' log-scale names expected by fit_pella_tomlinson_model options.
#' @noRd
.normalise_em_fixed_params <- function(fixed_params) {
  if (is.null(fixed_params) || length(fixed_params) == 0) {
    return(NULL)
  }

  out <- list()
  for (nm in names(fixed_params)) {
    if (!nzchar(nm)) next

    val <- fixed_params[[nm]]
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val)) {
      stop("All fixed_params values must be finite numeric scalars", call. = FALSE)
    }

    if (grepl("^log_", nm)) {
      out[[nm]] <- as.numeric(val)
      next
    }

    # Core SPM parameters are estimated on the log scale.
    if (val <= 0) {
      stop("Natural-scale fixed parameter values must be positive", call. = FALSE)
    }
    out[[paste0("log_", nm)]] <- log(as.numeric(val))
  }

  if (length(out) == 0) {
    return(NULL)
  }
  out
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

  fit_options <- list(
    silent = TRUE,
    validate_data = FALSE,
    show_starting_values_message = FALSE
  )

  if (!is.null(em_config)) {
    fit_options$process_noise <- isTRUE(em_config$process_noise)
    fit_options$process_error_structure <-
      tolower(as.character(em_config$process_error_structure %||% "iid"))

    fixed_params <- .normalise_em_fixed_params(em_config$fixed_params)
    if (!is.null(fixed_params)) {
      fit_options$fixed_params <- fixed_params
    }

    # Add control options if specified in EM config
    if (!is.null(em_config$control)) {
      fit_options$control <- em_config$control
    }

    # Add n_starts if specified
    if (!is.null(em_config$n_starts)) {
      fit_options$n_starts <- em_config$n_starts
    }

    # Add calculate_se if specified
    if (!is.null(em_config$calculate_se)) {
      fit_options$calculate_se <- em_config$calculate_se
    }

    # Add priors if specified
    if (!is.null(em_config$priors)) {
      fit_options$priors <- em_config$priors
    }
  }

  tryCatch(
    withCallingHandlers(
      {
        fit <- fit_pella_tomlinson_model(
          data,
          options = fit_options
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

        # Extract K. Accept both aggregate (K) and area-specific (K.A1, K.A2, ...)
        # parameterisations and sum to total carrying capacity.
        k_names <- grep("^K($|\\.)", names(fit$parameters), value = TRUE)
        if (length(k_names) == 0) {
          return(NULL)
        }
        total_K <- sum(fit$parameters[k_names])

        list(
          est_biomass = current_b,
          K           = total_K,
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

  B_initial <- rep_len(tp$B_initial, n_areas)

  # State
  biomass <- B_initial
  impl_eps <- rep(0, n_areas)
  process_state <- rep(0, n_areas)
  obs_eps <- NULL # AR(1) observation-error state; NULL triggers stationary initialisation

  # Storage
  biomass_store <- matrix(NA_real_, n_proj_years, n_areas)
  catch_store <- matrix(NA_real_, n_proj_years, n_areas)
  tac_store <- rep(NA_real_, n_proj_years)
  est_bio_store <- rep(NA_real_, n_proj_years)
  hcr_bio_store <- rep(NA_real_, n_proj_years)
  em_fit_success_store <- rep(NA_real_, n_proj_years)

  # Accumulated histories for EM
  cpue_history <- data.frame(year = integer(0), cpue = numeric(0))
  catch_history <- data.frame(year = integer(0), catch = numeric(0))

  # Management state
  current_tac <- initial_tac
  last_em_result <- NULL
  assess_freq <- scenario$assessment_frequency

  # Optional fixed catch-allocation weights by area.
  # If omitted, catch is allocated by current biomass share.
  alloc_weights <- NULL
  if (!is.null(scenario$catch_allocation)) {
    alloc_raw <- scenario$catch_allocation

    if (length(alloc_raw) != n_areas) {
      stop(
        "scenario$catch_allocation length (", length(alloc_raw),
        ") must equal n_areas (", n_areas, ")",
        call. = FALSE
      )
    }

    area_names <- NULL
    if (!is.null(names(tp$B_initial)) && length(tp$B_initial) == n_areas) {
      area_names <- names(tp$B_initial)
    } else if (!is.null(om_config$movement_cost_matrix) &&
      !is.null(rownames(om_config$movement_cost_matrix)) &&
      nrow(om_config$movement_cost_matrix) == n_areas) {
      area_names <- rownames(om_config$movement_cost_matrix)
    } else if (!is.null(names(tp$q)) && length(tp$q) == n_areas) {
      area_names <- names(tp$q)
    }

    if (n_areas > 1) {
      nm <- names(alloc_raw)
      if (is.null(nm) || any(!nzchar(nm))) {
        stop(
          "scenario$catch_allocation must be a named vector for multi-area OMs",
          call. = FALSE
        )
      }
      if (is.null(area_names) || any(!nzchar(area_names))) {
        stop(
          "Could not determine OM area names to align catch_allocation",
          call. = FALSE
        )
      }

      missing_names <- setdiff(area_names, nm)
      extra_names <- setdiff(nm, area_names)
      if (length(missing_names) > 0 || length(extra_names) > 0) {
        stop(
          "scenario$catch_allocation names must match OM areas exactly. Missing: ",
          paste(missing_names, collapse = ", "),
          "; Extra: ",
          paste(extra_names, collapse = ", "),
          call. = FALSE
        )
      }

      alloc_weights <- as.numeric(alloc_raw[area_names])
    } else {
      alloc_weights <- as.numeric(alloc_raw)
    }

    if (sum(alloc_weights) <= 0 || any(!is.finite(alloc_weights))) {
      stop("scenario$catch_allocation must contain finite non-negative values with positive sum",
        call. = FALSE
      )
    }
    alloc_weights <- alloc_weights / sum(alloc_weights)
  }

  q_vec <- rep_len(as.numeric(tp$q), n_areas)

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

      em_fit_success_store[yr] <- if (!is.null(em_result)) 1 else 0

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
    } else {
      # Carry forward latest assessment outcome between assessment years.
      if (yr > 1) {
        em_fit_success_store[yr] <- em_fit_success_store[yr - 1]
      }
    }

    tac_store[yr] <- current_tac
    em_est_vec <- if (!is.null(last_em_result)) suppressWarnings(as.numeric(last_em_result$est_biomass)) else numeric(0)
    est_bio_store[yr] <- if (length(em_est_vec) > 0 && all(is.finite(em_est_vec))) {
      sum(em_est_vec)
    } else {
      NA_real_
    }

    hcr_est_vec <- if (is.finite(est_bio_store[yr])) {
      est_bio_store[yr]
    } else {
      sum(biomass)
    }
    hcr_bio_store[yr] <- if (length(hcr_est_vec) > 0 && all(is.finite(hcr_est_vec))) {
      sum(hcr_est_vec)
    } else {
      NA_real_
    }

    if (yr == 1 && is.na(em_fit_success_store[yr])) {
      em_fit_success_store[yr] <- 0
    }

    if (!is.finite(em_fit_success_store[yr])) {
      em_fit_success_store[yr] <- if (yr > 1) em_fit_success_store[yr - 1] else 0
    }

    est_vec <- suppressWarnings(as.numeric(est_bio_store[yr]))
    est_bio_store[yr] <- if (length(est_vec) > 0 && all(is.finite(est_vec))) {
      sum(est_vec)
    } else {
      NA_real_
    }

    # 2. Implementation error → realized catch
    if (n_areas > 1) {
      if (!is.null(alloc_weights)) {
        tac_area <- current_tac * alloc_weights
      } else {
        tac_area <- current_tac * (biomass / pmax(sum(biomass), 1e-8))
      }
    } else {
      tac_area <- current_tac
    }

    impl_result <- apply_implementation_error(
      tac_area, scenario$implementation_error, impl_eps
    )
    catch <- impl_result$catch
    impl_eps <- impl_result$eps
    catch_store[yr, ] <- catch

    # 3. Generate CPUE observation from total biomass (pre-fishing).
    # Catchability is catch-weighted across areas so that the observed
    # index is dominated by the areas where fishing effort is concentrated.
    # If total catch is zero, fall back to biomass weighting.
    b_total <- sum(biomass)
    q_agg <- if (n_areas > 1) {
      catch_total <- sum(catch)
      if (catch_total > 0) {
        sum(q_vec * catch) / catch_total
      } else {
        sum(q_vec * biomass) / pmax(b_total, 1e-8)
      }
    } else {
      q_vec[[1L]]
    }
    cpue_step <- .simulate_cpue_step(
      true_biomass_total = b_total,
      obs_error_params   = obs_error_params,
      q                  = q_agg,
      year               = yr,
      previous_eps       = obs_eps
    )
    obs_eps <- cpue_step$eps
    cpue_history <- rbind(cpue_history, cpue_step$cpue_row)
    catch_history <- rbind(
      catch_history,
      data.frame(year = as.integer(yr), catch = sum(catch))
    )

    # 4. Project biomass forward
    step <- project_biomass(
      biomass,
      catch,
      om_config,
      movement_kernel,
      process_state = process_state,
      return_process_state = TRUE
    )
    biomass <- step$biomass
    process_state <- step$process_state
  }

  # Derived: harvest rate
  harvest_rate_store <- catch_store / pmax(biomass_store, 1e-8)

  list(
    biomass           = biomass_store,
    catch             = catch_store,
    tac               = tac_store,
    estimated_biomass = est_bio_store,
    hcr_biomass_used  = hcr_bio_store,
    em_fit_success    = em_fit_success_store,
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
  hcr_bio <- matrix(NA_real_, n_sims, n_proj_years)
  em_fit_success <- matrix(NA_real_, n_sims, n_proj_years)

  for (i in seq_len(n_sims)) {
    sim <- sim_results[[i]]
    biomass[i, , ] <- sim$biomass
    catch_arr[i, , ] <- sim$catch
    hr_arr[i, , ] <- sim$harvest_rate
    tac_mat[i, ] <- sim$tac
    est_bio[i, ] <- sim$estimated_biomass
    hcr_bio[i, ] <- sim$hcr_biomass_used
    em_fit_success[i, ] <- sim$em_fit_success
  }

  list(
    biomass           = biomass,
    catch             = catch_arr,
    harvest_rate      = hr_arr,
    tac               = tac_mat,
    estimated_biomass = est_bio,
    hcr_biomass_used  = hcr_bio,
    em_fit_success    = em_fit_success
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
