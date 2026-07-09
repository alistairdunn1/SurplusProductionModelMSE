# Performance Metrics (Module 2.4)
#
# Calculate MSE performance indicators from simulation trajectories.
# Includes probability-based risk metrics, catch statistics, and
# biomass reference point ratios.


#' Compute Equilibrium F at a Given Depletion Level
#'
#' For the Pella-Tomlinson model, calculate the fishing mortality that
#' would produce an equilibrium biomass of \code{x * K}.
#'
#' @param x Numeric scalar or vector. Depletion fraction(s) (e.g. 0.5
#'   for 50 percent of \code{K}).
#' @param r Numeric scalar. Intrinsic growth rate.
#' @param K Numeric scalar. Carrying capacity.
#' @param m Numeric scalar. Shape parameter.
#' @param B_unfished Numeric scalar. Unfished biomass (= K for the
#'   Pella-Tomlinson model). Default \code{K}.
#'
#' @return Named numeric vector of equilibrium F values. At equilibrium
#'   \eqn{F \cdot B = P(B)}, so
#'   \eqn{F = P(B) / B = \frac{r}{m-1} (1 - (B/K)^{m-1})}.
#'
#' @details
#' The equilibrium depletion F is:
#' \deqn{F_{x} = \frac{r}{m-1} \left(1 - (x B_\mathrm{unfished} / K)^{m-1}\right)}
#'
#' Returns 0 when \eqn{x B_\mathrm{unfished} \geq K} (no production surplus).
#'
#' @examples
#' # F that produces 50% depletion under Schaefer model
#' equilibrium_f(0.5, r = 0.3, K = 5000, m = 2)
#'
#' @export
equilibrium_f <- function(x, r, K, m, B_unfished = K) {
  assert_numeric(x,
    lower = 0, upper = 1, any.missing = FALSE,
    .var.name = "x"
  )
  assert_number(r, lower = 0)
  assert_number(K, lower = 0)
  assert_number(m, lower = 0)
  assert_number(B_unfished, lower = 0)

  B_target <- x * B_unfished
  f_eq <- ifelse(
    B_target >= K,
    0,
    SurplusProductionModel::pt_equilibrium_f(r = r, K = K, m = m, biomass = B_target)
  )
  names(f_eq) <- paste0("F", x * 100, "%K")
  f_eq
}


#' Calculate Performance Metrics from MSE Trajectories
#'
#' Compute a comprehensive set of performance indicators from simulation
#' output, including risk probabilities, catch statistics, and biomass
#' ratios. All probability metrics are computed per year, for the final
#' year, and as probability of ever breaching the threshold.
#'
#' @param trajectories A list containing simulation output matrices. Each
#'   element should be a matrix \code{[n_sims x n_years]} for single-area
#'   models, or a 3D array \code{[n_sims x n_years x n_areas]} for
#'   multi-area models. Required elements:
#'   \describe{
#'     \item{biomass}{True biomass from the operating model.}
#'     \item{catch}{Realized catch.}
#'   }
#'   Optional elements: \code{harvest_rate} (computed from
#'   \code{catch / biomass} if absent).
#' @param om_config An \code{\link{om_config}} object with
#'   \code{true_params} set (needed for reference points).
#' @param thresholds Numeric vector of depletion thresholds for risk
#'   metrics, as fractions of the unfished biomass \code{B0}. Without movement
#'   \code{B0} equals the sum of per-area carrying capacities; under movement it
#'   is the joint unfished spatial equilibrium (which differs from the sum of
#'   per-area \code{K}). Default \code{c(0.5, 0.2)}.
#' @param start_year Integer scalar. First projection year to include in
#'   performance summaries.
#'
#' @return An object of class \code{mse_performance} containing:
#'   \describe{
#'     \item{biomass_risk}{List of biomass risk metrics per threshold.
#'       Each has \code{per_year}, \code{final_year}, and \code{ever}
#'       components (aggregate and per-area).}
#'     \item{f_risk}{List of F risk metrics per threshold (same
#'       structure as \code{biomass_risk}).}
#'     \item{catch_stats}{List with \code{mean_catch},
#'       \code{median_catch}, \code{sd_catch}, and \code{aav}
#'       (aggregate and per-area).}
#'     \item{biomass_ratios}{List with \code{mean_depletion} (B/K),
#'       \code{mean_b_bmsy} (B/BMSY), \code{final_depletion},
#'       \code{final_b_bmsy} (aggregate and per-area).}
#'     \item{reference}{List of reference values: \code{K}, \code{K_area},
#'       \code{B_initial}, \code{BMSY}, \code{FMSY}, \code{MSY},
#'       \code{thresholds}, \code{f_thresholds}.}
#'     \item{n_sims}{Number of simulations.}
#'     \item{n_years}{Number of projection years.}
#'     \item{n_areas}{Number of areas.}
#'   }
#'
#' @export
calculate_performance_metrics <- function(trajectories,
                                          om_config,
                                          thresholds = c(0.5, 0.2),
                                          start_year = 1L) {
  assert_list(trajectories, .var.name = "trajectories")
  if (is.null(trajectories$biomass)) {
    stop("trajectories must contain a 'biomass' element", call. = FALSE)
  }
  if (is.null(trajectories$catch)) {
    stop("trajectories must contain a 'catch' element", call. = FALSE)
  }
  if (!inherits(om_config, "om_config")) {
    stop("om_config must be an 'om_config' object", call. = FALSE)
  }
  tp <- om_config$true_params
  if (is.null(tp)) {
    stop("om_config$true_params must be set", call. = FALSE)
  }

  assert_numeric(thresholds,
    lower = 0, upper = 1, any.missing = FALSE,
    min.len = 1, .var.name = "thresholds"
  )
  assert_count(start_year, positive = TRUE, .var.name = "start_year")

  # --- Coerce inputs to 3D: [n_sims x n_years x n_areas] ---
  biomass <- .to_3d(trajectories$biomass)
  catch <- .to_3d(trajectories$catch)

  n_sims <- dim(biomass)[1]
  n_years <- dim(biomass)[2]
  n_areas <- dim(biomass)[3]

  start_year <- as.integer(start_year)
  if (start_year > n_years) {
    stop("start_year must be less than or equal to the number of projected years", call. = FALSE)
  }
  eval_index <- seq.int(from = start_year, to = n_years)
  biomass <- biomass[, eval_index, , drop = FALSE]
  catch <- catch[, eval_index, , drop = FALSE]

  # Harvest rate: use provided or compute
  if (!is.null(trajectories$harvest_rate)) {
    harvest_rate <- .to_3d(trajectories$harvest_rate)
    harvest_rate <- harvest_rate[, eval_index, , drop = FALSE]
  } else {
    harvest_rate <- catch / pmax(biomass, 1e-8)
  }

  n_years <- dim(biomass)[2]

  # --- Reference points ---
  r <- tp$r
  K <- tp$K
  m <- tp$m

  # Status baseline is the unfished biomass B0 (the quantity denoted B0 in
  # standard surplus-production notation). Under movement the joint unfished
  # equilibrium differs from the sum of per-area carrying capacities, so
  # depletion and biomass risk are expressed relative to that equilibrium;
  # without movement B0 = sum of per-area K. The aggregate carrying capacity
  # K_cc is retained separately for the (production-based) reference points.
  K_cc <- tp$K
  B_initial_area <- rep_len(tp$B_initial, n_areas)
  K_area <- .om_unfished_biomass(om_config) # per-area unfished biomass (B0_a)
  K_total <- sum(K_area) # total unfished biomass (B0)

  # Standard Pella-Tomlinson reference points, shared with the assessment
  # package so OM truth and EM estimates use one implementation:
  # BMSY = K * m^(-1/(m-1)), FMSY = r/m, MSY = FMSY * BMSY
  # (Fox limit m -> 1: BMSY = K/e, FMSY = r, MSY = rK/e).
  rp <- SurplusProductionModel::pella_tomlinson_reference_points(
    r = r, K = K, m = m
  )
  BMSY <- rp$bmsy
  FMSY <- rp$fmsy
  MSY <- rp$msy

  # F thresholds: equilibrium F at each depletion level, using the
  # status baseline biomass reference.
  f_thresholds <- equilibrium_f(thresholds, r = r, K = K, m = m, B_unfished = K_total)

  # --- Aggregate biomass/catch/F: sum across areas ---
  biomass_agg <- apply(biomass, c(1, 2), sum) # [n_sims x n_years]
  catch_agg <- apply(catch, c(1, 2), sum)
  # Aggregate F: biomass-weighted mean across areas
  harvest_rate_agg <- apply(harvest_rate * biomass, c(1, 2), sum) /
    pmax(biomass_agg, 1e-8)

  # --- Biomass risk metrics ---
  biomass_risk <- list()
  for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]
    lbl <- paste0(thresh * 100, "%K")

    # Aggregate
    B_limit_agg <- thresh * K_total
    agg <- .risk_metrics(biomass_agg, B_limit_agg)

    # Per-area
    area_results <- vector("list", n_areas)
    for (a in seq_len(n_areas)) {
      B_limit_a <- thresh * K_area[a]
      area_results[[a]] <- .risk_metrics(biomass[, , a], B_limit_a)
    }
    names(area_results) <- paste0("A", seq_len(n_areas))

    biomass_risk[[lbl]] <- list(aggregate = agg, by_area = area_results)
  }

  # --- F risk metrics ---
  f_risk <- list()
  for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]
    lbl <- paste0("F", thresh * 100, "%K")
    f_limit <- f_thresholds[i]

    # For F risk, we want Pr(F > F_limit) — probability that fishing
    # mortality EXCEEDS the limit F. This is the "risk" interpretation:
    # the F threshold is a limit, and exceeding it is dangerous.
    agg <- .risk_metrics_above(harvest_rate_agg, f_limit)

    area_results <- vector("list", n_areas)
    for (a in seq_len(n_areas)) {
      area_results[[a]] <- .risk_metrics_above(harvest_rate[, , a], f_limit)
    }
    names(area_results) <- paste0("A", seq_len(n_areas))

    f_risk[[lbl]] <- list(aggregate = agg, by_area = area_results)
  }

  # --- Catch statistics ---
  catch_stats <- list(
    aggregate = .catch_stats(catch_agg),
    by_area = lapply(seq_len(n_areas), function(a) {
      .catch_stats(catch[, , a])
    })
  )
  names(catch_stats$by_area) <- paste0("A", seq_len(n_areas))

  # --- Biomass ratios ---
  depletion_agg <- biomass_agg / K_total
  b_bmsy_agg <- biomass_agg / BMSY

  biomass_ratios <- list(
    aggregate = list(
      mean_depletion  = mean(depletion_agg),
      mean_b_bmsy     = mean(b_bmsy_agg),
      final_depletion = mean(depletion_agg[, n_years, drop = FALSE]),
      final_b_bmsy    = mean(b_bmsy_agg[, n_years, drop = FALSE])
    ),
    by_area = lapply(seq_len(n_areas), function(a) {
      # Extract area, ensuring we get a [n_sims x n_years] matrix
      biomass_a <- biomass[, , a, drop = TRUE]
      if (!is.matrix(biomass_a) && length(biomass_a) > 0) {
        biomass_a <- matrix(biomass_a, nrow = 1)
      }

      dep_a <- biomass_a / K_area[a]
      bb_a <- biomass_a / (BMSY * K_area[a] / K_total)
      list(
        mean_depletion  = mean(dep_a, na.rm = TRUE),
        mean_b_bmsy     = mean(bb_a, na.rm = TRUE),
        final_depletion = mean(dep_a[, n_years, drop = FALSE], na.rm = TRUE),
        final_b_bmsy    = mean(bb_a[, n_years, drop = FALSE], na.rm = TRUE)
      )
    })
  )
  names(biomass_ratios$by_area) <- paste0("A", seq_len(n_areas))

  structure(
    list(
      biomass_risk = biomass_risk,
      f_risk = f_risk,
      catch_stats = catch_stats,
      biomass_ratios = biomass_ratios,
      reference = list(
        K                 = K_total, # status baseline B0 (unfished biomass)
        B0                = K_total, # explicit unfished-biomass alias
        K_area            = K_area, # per-area unfished biomass (B0_a)
        carrying_capacity = K_cc, # aggregate carrying capacity (production K)
        B_initial         = sum(B_initial_area), # OM initial biomass (start of series)
        BMSY              = BMSY,
        FMSY              = FMSY,
        MSY               = MSY,
        thresholds        = thresholds,
        f_thresholds      = f_thresholds
      ),
      n_sims = n_sims,
      n_years = n_years,
      n_areas = n_areas
    ),
    class = "mse_performance"
  )
}


# --- Internal helpers ---

#' Operating-model unfished spatial equilibrium (B0 per area)
#'
#' Deterministic, zero-catch biomass per area at equilibrium. Under movement
#' this is the joint redistributed equilibrium (production balances net
#' movement) and differs from the per-area carrying capacities; without movement
#' each area rests at its own K. Used as the unfished baseline B0 for depletion
#' and status so that these are measured against the true unfished level.
#' @noRd
.om_unfished_biomass <- function(om_config, n_spinup = 100L) {
  n_areas <- om_config$n_areas
  tp <- om_config$true_params

  B_initial_area <- rep_len(tp$B_initial, n_areas)
  if (!is.null(tp$K_area)) {
    K_area <- rep_len(as.numeric(tp$K_area), n_areas)
  } else if (n_areas > 1) {
    K_area <- tp$K * (B_initial_area / sum(B_initial_area))
  } else {
    K_area <- rep_len(as.numeric(tp$K), n_areas)
  }

  # No movement: each area rests at its own carrying capacity (B0 = K).
  if (n_areas == 1L || is.null(om_config$movement_rate) ||
    om_config$movement_rate <= 0) {
    return(K_area)
  }

  # With movement, iterate the deterministic zero-catch dynamics to the joint
  # unfished equilibrium using the operating model's own transition function.
  traj <- project_trajectory(
    B_initial     = K_area,
    catch_series  = matrix(0, nrow = n_spinup, ncol = n_areas),
    om_config     = om_config,
    n_years       = n_spinup,
    process_noise = FALSE
  )
  as.numeric(traj[nrow(traj), , drop = TRUE])
}


#' Coerce to 3D array \[n_sims x n_years x n_areas\]
#' @noRd
.to_3d <- function(x) {
  if (is.array(x) && length(dim(x)) == 3) {
    return(x)
  }
  if (is.matrix(x)) {
    # [n_sims x n_years] -> [n_sims x n_years x 1]
    return(array(x, dim = c(nrow(x), ncol(x), 1)))
  }
  stop("trajectory element must be a matrix [n_sims x n_years] or ",
    "3D array [n_sims x n_years x n_areas]",
    call. = FALSE
  )
}


#' Compute risk metrics: Pr(value < threshold) per year, final year, ever
#' @param mat Matrix \[n_sims x n_years\]
#' @param threshold Numeric scalar
#' @return List with per_year, final_year, ever
#' @noRd
.risk_metrics <- function(mat, threshold) {
  # Ensure mat is a matrix (handle dimension dropping edge cases)
  if (!is.matrix(mat)) {
    if (is.numeric(mat) && length(mat) > 0) {
      mat <- matrix(mat, nrow = 1)
    } else {
      stop("Invalid input to .risk_metrics: mat must be a numeric matrix")
    }
  }

  # Handle empty matrix
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(list(per_year = numeric(0), final_year = NA_real_, ever = NA_real_))
  }

  below <- mat < threshold
  per_year <- colMeans(below, na.rm = TRUE) # probability per year
  final_year <- per_year[ncol(mat)] # final year only
  ever <- mean(apply(below, 1, any, na.rm = TRUE), na.rm = TRUE) # prob of ever breaching

  list(per_year = per_year, final_year = final_year, ever = ever)
}


#' Compute risk metrics: Pr(value > threshold) per year, final year, ever
#' @param mat Matrix \[n_sims x n_years\]
#' @param threshold Numeric scalar
#' @return List with per_year, final_year, ever
#' @noRd
.risk_metrics_above <- function(mat, threshold) {
  # Ensure mat is a matrix (handle dimension dropping edge cases)
  if (!is.matrix(mat)) {
    if (is.numeric(mat) && length(mat) > 0) {
      mat <- matrix(mat, nrow = 1)
    } else {
      stop("Invalid input to .risk_metrics_above: mat must be a numeric matrix")
    }
  }

  # Handle empty matrix
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(list(per_year = numeric(0), final_year = NA_real_, ever = NA_real_))
  }

  above <- mat > threshold
  per_year <- colMeans(above, na.rm = TRUE)
  final_year <- per_year[ncol(mat)]
  ever <- mean(apply(above, 1, any, na.rm = TRUE), na.rm = TRUE)

  list(per_year = per_year, final_year = final_year, ever = ever)
}


#' Compute catch statistics from a \[n_sims x n_years\] matrix
#' @noRd
.catch_stats <- function(mat) {
  # Ensure mat is a matrix (handle dimension dropping edge cases)
  if (!is.matrix(mat)) {
    if (is.numeric(mat) && length(mat) > 0) {
      mat <- matrix(mat, nrow = 1)
    } else {
      stop("Invalid input to .catch_stats: mat must be a numeric matrix")
    }
  }

  # Handle empty matrix
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(list(
      mean_catch = NA_real_, median_catch = NA_real_,
      sd_catch = NA_real_, aav = NA_real_
    ))
  }

  # Average Annual Variation in catch. Computed per simulation as the sum of
  # absolute year-on-year changes divided by the total catch, then averaged
  # across simulations: AAV = mean_s( sum_t |C_{t+1} - C_t| / sum_t C_t ). The
  # sum-based form is robust to individual near-zero catch years, which make the
  # per-year ratio form (dividing each change by that year's catch) explode.
  n_years <- ncol(mat)
  if (n_years > 1) {
    abs_change <- rowSums(
      abs(mat[, -1, drop = FALSE] - mat[, -n_years, drop = FALSE]),
      na.rm = TRUE
    )
    total_catch <- rowSums(mat, na.rm = TRUE)
    aav <- mean(abs_change / pmax(total_catch, 1e-8), na.rm = TRUE)
  } else {
    aav <- 0
  }

  list(
    mean_catch   = mean(mat, na.rm = TRUE),
    median_catch = median(mat, na.rm = TRUE),
    sd_catch     = sd(as.vector(mat), na.rm = TRUE),
    aav          = aav
  )
}


#' @export
print.mse_performance <- function(x, ...) {
  cat("MSE Performance Metrics\n")
  cat("=======================\n")
  cat(
    "Simulations:", x$n_sims, " | Years:", x$n_years,
    " | Areas:", x$n_areas, "\n\n"
  )

  # Biomass risk
  cat("Biomass Risk (Pr(B < threshold)):\n")
  for (nm in names(x$biomass_risk)) {
    agg <- x$biomass_risk[[nm]]$aggregate
    cat(sprintf(
      "  %-10s  final_year: %.3f  ever: %.3f\n",
      nm, agg$final_year, agg$ever
    ))
  }

  # F risk
  cat("\nF Risk (Pr(F > F_limit)):\n")
  for (nm in names(x$f_risk)) {
    agg <- x$f_risk[[nm]]$aggregate
    f_val <- x$reference$f_thresholds[match(nm, names(x$reference$f_thresholds))]
    cat(sprintf(
      "  %-12s (F=%.4f)  final_year: %.3f  ever: %.3f\n",
      nm, f_val, agg$final_year, agg$ever
    ))
  }

  # Catch
  cs <- x$catch_stats$aggregate
  cat(sprintf(
    "\nCatch: mean=%.1f  median=%.1f  sd=%.1f  AAV=%.3f\n",
    cs$mean_catch, cs$median_catch, cs$sd_catch, cs$aav
  ))

  # Biomass ratios
  br <- x$biomass_ratios$aggregate
  cat(sprintf(
    "\nBiomass Ratios (aggregate): mean B/K=%.3f  ",
    br$mean_depletion
  ))
  cat(sprintf("mean B/BMSY=%.3f\n", br$mean_b_bmsy))
  cat(sprintf(
    "  Final year: B/K=%.3f  B/BMSY=%.3f\n",
    br$final_depletion, br$final_b_bmsy
  ))

  # Reference points
  ref <- x$reference
  cat(sprintf(
    "\nReference: K=%.0f  BMSY=%.0f  FMSY=%.4f  MSY=%.0f\n",
    ref$K, ref$BMSY, ref$FMSY, ref$MSY
  ))

  invisible(x)
}


#' Summarise Performance Metrics as a Data Frame
#'
#' Convert an \code{mse_performance} object into a tidy data frame
#' suitable for comparison across scenarios.
#'
#' @param object An \code{mse_performance} object.
#' @param ... Additional arguments (not used).
#'
#' @return A data frame with one row per metric, including columns
#'   \code{metric}, \code{scope} (aggregate or area name),
#'   \code{value}.
#'
#' @export
summary.mse_performance <- function(object, ...) {
  rows <- list()

  # Biomass risk
  for (nm in names(object$biomass_risk)) {
    agg <- object$biomass_risk[[nm]]$aggregate
    rows[[length(rows) + 1]] <- data.frame(
      metric = paste0("Pr(B<", nm, ")_final"),
      scope = "aggregate", value = agg$final_year,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      metric = paste0("Pr(B<", nm, ")_ever"),
      scope = "aggregate", value = agg$ever,
      stringsAsFactors = FALSE
    )
    for (anm in names(object$biomass_risk[[nm]]$by_area)) {
      ar <- object$biomass_risk[[nm]]$by_area[[anm]]
      rows[[length(rows) + 1]] <- data.frame(
        metric = paste0("Pr(B<", nm, ")_final"),
        scope = anm, value = ar$final_year,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- data.frame(
        metric = paste0("Pr(B<", nm, ")_ever"),
        scope = anm, value = ar$ever,
        stringsAsFactors = FALSE
      )
    }
  }

  # F risk
  for (nm in names(object$f_risk)) {
    agg <- object$f_risk[[nm]]$aggregate
    rows[[length(rows) + 1]] <- data.frame(
      metric = paste0("Pr(F>", nm, ")_final"),
      scope = "aggregate", value = agg$final_year,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      metric = paste0("Pr(F>", nm, ")_ever"),
      scope = "aggregate", value = agg$ever,
      stringsAsFactors = FALSE
    )
    for (anm in names(object$f_risk[[nm]]$by_area)) {
      ar <- object$f_risk[[nm]]$by_area[[anm]]
      rows[[length(rows) + 1]] <- data.frame(
        metric = paste0("Pr(F>", nm, ")_final"),
        scope = anm, value = ar$final_year,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- data.frame(
        metric = paste0("Pr(F>", nm, ")_ever"),
        scope = anm, value = ar$ever,
        stringsAsFactors = FALSE
      )
    }
  }

  # Catch stats
  cs <- object$catch_stats$aggregate
  rows[[length(rows) + 1]] <- data.frame(
    metric = "mean_catch", scope = "aggregate",
    value = cs$mean_catch, stringsAsFactors = FALSE
  )
  rows[[length(rows) + 1]] <- data.frame(
    metric = "AAV", scope = "aggregate",
    value = cs$aav, stringsAsFactors = FALSE
  )

  # Biomass ratios
  br <- object$biomass_ratios$aggregate
  rows[[length(rows) + 1]] <- data.frame(
    metric = "mean_B_K", scope = "aggregate",
    value = br$mean_depletion, stringsAsFactors = FALSE
  )
  rows[[length(rows) + 1]] <- data.frame(
    metric = "mean_B_BMSY", scope = "aggregate",
    value = br$mean_b_bmsy, stringsAsFactors = FALSE
  )
  rows[[length(rows) + 1]] <- data.frame(
    metric = "final_B_K", scope = "aggregate",
    value = br$final_depletion, stringsAsFactors = FALSE
  )
  rows[[length(rows) + 1]] <- data.frame(
    metric = "final_B_BMSY", scope = "aggregate",
    value = br$final_b_bmsy, stringsAsFactors = FALSE
  )

  do.call(rbind, rows)
}
