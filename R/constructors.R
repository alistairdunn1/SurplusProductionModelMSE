# Constructor and print methods for MSE S3 configuration classes


# ── om_config ──────────────────────────────────────────────────────────────────

#' Configure the Operating Model
#'
#' Create an operating model configuration specifying the "true" population
#' dynamics for MSE simulations.
#'
#' @param n_areas Integer >= 1. Number of spatial areas (default 1 = non-spatial).
#' @param movement_rate Numeric in \[0, 1\]. Fraction of biomass redistributed per
#'   time step (default 0 = no movement).
#' @param movement_cost_matrix Square matrix of inter-area movement costs,
#'   or \code{NULL} for single area. Diagonal must be zero.
#' @param attractiveness Numeric vector of length \code{n_areas} giving habitat
#'   quality weights, or \code{NULL} (all areas equal).
#' @param decay Numeric >= 0. Distance decay parameter for the gravity movement
#'   kernel (default 0).
#' @param true_params Named list of true parameter values containing at minimum
#'   \code{r}, \code{K}, \code{m}, \code{sigma_obs}, \code{q}, \code{B_initial}.
#'   \code{q} and \code{B_initial} may be scalar or per-area vectors. If \code{NULL},
#'   parameters must be supplied separately.
#' @param max_harvest_rate Optional maximum realised annual exploitation rate.
#'   \code{NULL} (default) applies no additional catch cap. Supply a value in
#'   \code{(0, 1]} only to represent a documented implementation constraint.
#'
#' @return An S3 object of class \code{om_config}.
#'
#' @examples
#' # Single-area operating model
#' om <- om_config(
#'   true_params = list(
#'     r = 0.3, K = 5000, m = 2,
#'     sigma_obs = 0.2, q = 1e-4, B_initial = 4000
#'   )
#' )
#'
#' # Multi-area with movement
#' om3 <- om_config(
#'   n_areas = 3,
#'   movement_rate = 0.1,
#'   movement_cost_matrix = matrix(c(0, 100, 200, 100, 0, 100, 200, 100, 0), 3, 3),
#'   attractiveness = c(1, 1.2, 0.8),
#'   decay = 0.01,
#'   true_params = list(
#'     r = 0.3, K = 5000, m = 2,
#'     sigma_obs = 0.2, q = rep(1e-4, 3),
#'     B_initial = c(2000, 2500, 1500)
#'   )
#' )
#'
#' @export
om_config <- function(n_areas = 1L,
                      movement_rate = 0,
                      movement_cost_matrix = NULL,
                      attractiveness = NULL,
                      decay = 0,
                      true_params = NULL,
                      max_harvest_rate = NULL) {
  assert_count(n_areas, positive = TRUE, .var.name = "n_areas")

  obj <- structure(
    list(
      n_areas = as.integer(n_areas),
      movement_rate = movement_rate,
      movement_cost_matrix = movement_cost_matrix,
      attractiveness = attractiveness,
      decay = decay,
      true_params = true_params,
      max_harvest_rate = max_harvest_rate
    ),
    class = "om_config"
  )
  validate_om_config(obj)
  obj
}


#' @export
print.om_config <- function(x, ...) {
  cat("Operating Model Configuration\n")
  cat("-----------------------------\n")
  cat("  Areas:         ", x$n_areas, "\n")
  cat("  Movement rate: ", x$movement_rate, "\n")
  if (!is.null(x$movement_cost_matrix)) {
    cat("  Movement cost matrix: ", x$n_areas, "x", x$n_areas, " supplied\n")
  }
  if (!is.null(x$attractiveness)) {
    cat(
      "  Attractiveness:", paste(round(x$attractiveness, 3), collapse = ", "),
      "\n"
    )
  }
  cat("  Decay:         ", x$decay, "\n")
  if (!is.null(x$max_harvest_rate)) {
    cat("  Maximum realised U:", x$max_harvest_rate, "\n")
  }
  if (!is.null(x$true_params)) {
    tp <- x$true_params
    cat("  True parameters:\n")
    cat("    r =", tp$r, " K =", tp$K, " m =", tp$m, "\n")
    cat("    sigma_obs =", tp$sigma_obs, "\n")
    cat("    q =", paste(signif(tp$q, 4), collapse = ", "), "\n")
    cat("    B_initial =", paste(round(tp$B_initial, 1), collapse = ", "), "\n")
  } else {
    cat("  True parameters: not specified\n")
  }
  invisible(x)
}


# ── em_config ──────────────────────────────────────────────────────────────────

#' Configure the Estimation Model
#'
#' Create an estimation model configuration for MSE simulations. The estimation
#' model structure can differ from the operating model to test robustness to
#' model misspecification.
#'
#' @param n_areas Integer >= 1. Number of areas in the estimation model. May
#'   differ from the operating model for misspecification testing.
#' @param estimate_movement Logical. Whether the EM estimates movement parameters
#'   (default \code{FALSE}).
#' @param aggregate_areas Logical. If \code{TRUE} and the OM has multiple areas,
#'   aggregate observed data to a single-area EM input (default \code{FALSE}).
#' @param process_noise Logical. Whether to include process noise in EM fitting
#'   (default \code{FALSE}).
#' @param process_error_structure Character, one of \code{"iid"} or
#'   \code{"ar1"}. Used when \code{process_noise = TRUE}.
#' @param fixed_params Named list of parameters to fix (not estimate) in the EM.
#'   For example, \code{list(m = 2)} to fix the Schaefer shape. \code{NULL}
#'   means all parameters are estimated.
#' @param control Optional optimiser control list passed through to
#'   \code{\link[SurplusProductionModel]{fit_pella_tomlinson_model}}.
#' @param n_starts Integer >= 1. Number of random-restart optimiser runs
#'   passed through to
#'   \code{\link[SurplusProductionModel]{fit_pella_tomlinson_model}} (default
#'   \code{1}).
#' @param calculate_se Logical. Whether to compute parameter standard errors
#'   during EM fitting (default \code{FALSE}).
#' @param priors Optional named list of priors passed through to
#'   \code{\link[SurplusProductionModel]{fit_pella_tomlinson_model}}.
#' @param initial_depletion Optional starting biomass depletion for the
#'   estimation model. Supply an assessment-derived value when the simulated
#'   series begins from an already exploited stock. If \code{NULL}, initial
#'   depletion is estimated or controlled through \code{fixed_params}/priors.
#'
#' @return An S3 object of class \code{em_config}.
#'
#' @examples
#' # Single-area EM that fixes shape to Schaefer
#' em <- em_config(fixed_params = list(m = 2))
#'
#' # Aggregated EM (ignores spatial structure)
#' em_agg <- em_config(n_areas = 1, aggregate_areas = TRUE)
#'
#' @export
em_config <- function(n_areas = 1L,
                      estimate_movement = FALSE,
                      aggregate_areas = FALSE,
                      process_noise = FALSE,
                      process_error_structure = "iid",
                      fixed_params = NULL,
                      control = NULL,
                      n_starts = 1L,
                      calculate_se = FALSE,
                      priors = NULL,
                      initial_depletion = NULL) {
  assert_count(n_areas, positive = TRUE, .var.name = "n_areas")
  assert_flag(process_noise, .var.name = "process_noise")
  process_error_structure <- match.arg(
    tolower(as.character(process_error_structure)),
    c("iid", "ar1")
  )

  obj <- structure(
    list(
      n_areas = as.integer(n_areas),
      estimate_movement = estimate_movement,
      aggregate_areas = aggregate_areas,
      process_noise = process_noise,
      process_error_structure = process_error_structure,
      fixed_params = fixed_params,
      control = control,
      n_starts = as.integer(n_starts),
      calculate_se = calculate_se,
      priors = priors,
      initial_depletion = initial_depletion
    ),
    class = "em_config"
  )
  validate_em_config(obj)
  obj
}


#' @export
print.em_config <- function(x, ...) {
  cat("Estimation Model Configuration\n")
  cat("------------------------------\n")
  cat("  Areas:             ", x$n_areas, "\n")
  cat("  Estimate movement: ", x$estimate_movement, "\n")
  cat("  Aggregate areas:   ", x$aggregate_areas, "\n")
  cat("  Process noise:     ", x$process_noise, "\n")
  cat("  Process structure: ", x$process_error_structure, "\n")
  if (!is.null(x$fixed_params)) {
    cat(
      "  Fixed parameters:  ",
      paste(names(x$fixed_params), "=",
        unlist(x$fixed_params),
        collapse = ", "
      ), "\n"
    )
  } else {
    cat("  Fixed parameters:   none\n")
  }
  invisible(x)
}


# ── impl_error ─────────────────────────────────────────────────────────────────

#' Configure Implementation Error
#'
#' Specify the stochastic discrepancy between the TAC (total allowable catch)
#' recommended by the harvest control rule and the realized catch.
#'
#' @param cv Numeric > 0. Coefficient of variation for the lognormal error
#'   (default 0.1).
#' @param bias Numeric >= 0. Multiplicative bias: 1 = unbiased, > 1 = systematic
#'   over-catch, < 1 = under-catch (default 1).
#' @param autocorr Numeric in \[0, 1\]. Temporal autocorrelation (AR1 coefficient)
#'   in errors (default 0 = independent).
#' @param max_overage Numeric >= 1. Maximum allowed catch / TAC ratio (default
#'   1.1, i.e. 10 percent maximum overage).
#'
#' @return An S3 object of class \code{impl_error}.
#'
#' @examples
#' # Default: 10 pct CV, unbiased, no autocorrelation, max 10 pct overage
#' ie <- impl_error()
#'
#' # High autocorrelation with bias
#' ie2 <- impl_error(cv = 0.15, bias = 1.05, autocorr = 0.5)
#'
#' @export
impl_error <- function(cv = 0.1,
                       bias = 1,
                       autocorr = 0,
                       max_overage = 1.1) {
  obj <- structure(
    list(
      cv = cv,
      bias = bias,
      autocorr = autocorr,
      max_overage = max_overage
    ),
    class = "impl_error"
  )
  validate_impl_error(obj)
  obj
}


#' @export
print.impl_error <- function(x, ...) {
  cat("Implementation Error Specification\n")
  cat("----------------------------------\n")
  cat("  CV:            ", x$cv, "\n")
  cat(
    "  Bias:          ", x$bias,
    if (x$bias == 1) "(unbiased)" else if (x$bias > 1) "(over-catch)" else "(under-catch)", "\n"
  )
  cat("  Autocorrelation:", x$autocorr, "\n")
  cat("  Max overage:   ", x$max_overage, "\n")
  invisible(x)
}


# ── mse_scenario ───────────────────────────────────────────────────────────────

#' Create an MSE Scenario
#'
#' Define an MSE scenario consisting of a harvest control rule, optional
#' implementation error, and assessment frequency.
#'
#' @param name Character. Scenario identifier.
#' @param harvest_control_rule A function with signature
#'   \code{function(biomass, reference_points)} that returns a TAC value.
#' @param implementation_error An \code{\link{impl_error}} object, or
#'   \code{NULL} for perfect implementation (catch == TAC).
#' @param assessment_frequency Integer >= 1. Years between assessments
#'   (default 1 = annual).
#' @param catch_allocation Optional numeric vector of non-negative area weights
#'   used to split TAC across operating-model areas. If \code{NULL} (default),
#'   TAC is allocated in proportion to current biomass in each area. For
#'   multi-area operating models, provide this as a named vector so weights
#'   can be matched to area names.
#'
#' @return An S3 object of class \code{mse_scenario}.
#'
#' @examples
#' # Simple constant-F scenario with annual assessment
#' my_hcr <- function(biomass, reference_points) biomass * 0.1
#' sc <- create_scenario("ConstantF_0.1", my_hcr)
#'
#' # With implementation error and biennial assessment
#' sc2 <- create_scenario(
#'   "HockeyStick",
#'   my_hcr,
#'   implementation_error = impl_error(cv = 0.15),
#'   assessment_frequency = 2
#' )
#'
#' @export
create_scenario <- function(name,
                            harvest_control_rule,
                            implementation_error = NULL,
                            assessment_frequency = 1L,
                            catch_allocation = NULL) {
  obj <- structure(
    list(
      name                 = name,
      harvest_control_rule = harvest_control_rule,
      implementation_error = implementation_error,
      assessment_frequency = as.integer(assessment_frequency),
      catch_allocation     = catch_allocation
    ),
    class = "mse_scenario"
  )
  validate_mse_scenario(obj)
  obj
}


#' @export
print.mse_scenario <- function(x, ...) {
  cat("MSE Scenario:", x$name, "\n")
  cat("-----------------------------------\n")
  cat("  HCR:                  <function>\n")
  if (!is.null(x$implementation_error)) {
    cat(
      "  Implementation error:  CV =", x$implementation_error$cv,
      " bias =", x$implementation_error$bias, "\n"
    )
  } else {
    cat("  Implementation error:  none (perfect)\n")
  }
  cat(
    "  Assessment frequency: ", x$assessment_frequency,
    if (x$assessment_frequency == 1) "year" else "years", "\n"
  )
  if (!is.null(x$catch_allocation)) {
    cat(
      "  Catch allocation:     ",
      paste(round(x$catch_allocation, 4), collapse = ", "),
      "\n"
    )
  } else {
    cat("  Catch allocation:      biomass-proportional\n")
  }
  invisible(x)
}


#' Create baseline scenarios for MSE runs
#'
#' Construct the standard baseline scenario set used in workflow scripts:
#' no-catch, constant annual exploitation rate, and two hockey-stick variants.
#'
#' @param f_target Positive numeric scalar for target annual exploitation rate.
#'   The argument name is retained for backwards compatibility.
#' @param assessment_frequency Integer >= 1 giving assessment interval in years.
#' @param catch_allocation Optional numeric vector of non-negative area weights.
#'
#' @return A named list of `mse_scenario` objects.
#' @export
create_baseline_scenarios <- function(f_target,
                                      assessment_frequency = 1L,
                                      catch_allocation = NULL) {
  assert_number(f_target, lower = .Machine$double.eps, .var.name = "f_target")
  assert_count(assessment_frequency, positive = TRUE, .var.name = "assessment_frequency")

  if (!is.null(catch_allocation)) {
    assert_numeric(
      catch_allocation,
      any.missing = FALSE,
      lower = 0,
      finite = TRUE,
      min.len = 1,
      .var.name = "catch_allocation"
    )
    if (all(catch_allocation == 0)) {
      stop("catch_allocation must contain at least one positive value", call. = FALSE)
    }
  }

  list(
    create_scenario(
      name = "no_catch",
      harvest_control_rule = function(biomass, reference_points) {
        0
      },
      assessment_frequency = assessment_frequency,
      catch_allocation = catch_allocation
    ),
    create_scenario(
      name = sprintf("constant_f_%0.2f", f_target),
      harvest_control_rule = hcr_constant_f(f_target),
      assessment_frequency = assessment_frequency,
      catch_allocation = catch_allocation
    ),
    create_scenario(
      name = "hockey_20_50",
      harvest_control_rule = hcr_hockey_stick(
        f_target = f_target,
        b_limit = 0.20,
        b_target = 0.50
      ),
      assessment_frequency = assessment_frequency,
      catch_allocation = catch_allocation
    ),
    create_scenario(
      name = "hockey_10_40",
      harvest_control_rule = hcr_hockey_stick(
        f_target = f_target,
        b_limit = 0.10,
        b_target = 0.40
      ),
      assessment_frequency = assessment_frequency,
      catch_allocation = catch_allocation
    )
  )
}


#' Create an HCR grid scenario set for tuning
#'
#' Construct a flattened list of scenarios over an `f_grid` for selected
#' harvest-control-rule families.
#'
#' @param f_grid Numeric vector of positive annual exploitation-rate targets.
#'   The argument and family names are retained for backwards compatibility.
#' @param assessment_frequency Integer >= 1 giving assessment interval in years.
#' @param catch_allocation Optional numeric vector of non-negative area weights.
#' @param include Character vector of HCR families to include. Allowed values are
#'   `"constant_f"`, `"hockey_20_50"`, and `"hockey_10_40"`.
#'
#' @return A named list of `mse_scenario` objects.
#' @export
create_hcr_grid_scenarios <- function(f_grid,
                                      assessment_frequency = 1L,
                                      catch_allocation = NULL,
                                      include = c("constant_f", "hockey_20_50", "hockey_10_40")) {
  assert_numeric(
    f_grid,
    any.missing = FALSE,
    lower = .Machine$double.eps,
    finite = TRUE,
    min.len = 1,
    .var.name = "f_grid"
  )
  assert_count(assessment_frequency, positive = TRUE, .var.name = "assessment_frequency")
  assert_character(include, any.missing = FALSE, min.len = 1, .var.name = "include")

  if (!is.null(catch_allocation)) {
    assert_numeric(
      catch_allocation,
      any.missing = FALSE,
      lower = 0,
      finite = TRUE,
      min.len = 1,
      .var.name = "catch_allocation"
    )
    if (all(catch_allocation == 0)) {
      stop("catch_allocation must contain at least one positive value", call. = FALSE)
    }
  }

  include <- unique(include)
  allowed <- c("constant_f", "hockey_20_50", "hockey_10_40")
  bad <- setdiff(include, allowed)
  if (length(bad) > 0) {
    stop(
      "Unknown include option(s): ",
      paste(bad, collapse = ", "),
      ". Allowed values are: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }

  out <- lapply(f_grid, function(f_val) {
    scenarios <- list()

    if ("constant_f" %in% include) {
      scenarios[[length(scenarios) + 1L]] <- create_scenario(
        name = sprintf("cf_%0.3f", f_val),
        harvest_control_rule = hcr_constant_f(f_val),
        assessment_frequency = assessment_frequency,
        catch_allocation = catch_allocation
      )
    }

    if ("hockey_20_50" %in% include) {
      scenarios[[length(scenarios) + 1L]] <- create_scenario(
        name = sprintf("hs2050_%0.3f", f_val),
        harvest_control_rule = hcr_hockey_stick(
          f_target = f_val,
          b_limit = 0.20,
          b_target = 0.50
        ),
        assessment_frequency = assessment_frequency,
        catch_allocation = catch_allocation
      )
    }

    if ("hockey_10_40" %in% include) {
      scenarios[[length(scenarios) + 1L]] <- create_scenario(
        name = sprintf("hs1040_%0.3f", f_val),
        harvest_control_rule = hcr_hockey_stick(
          f_target = f_val,
          b_limit = 0.10,
          b_target = 0.40
        ),
        assessment_frequency = assessment_frequency,
        catch_allocation = catch_allocation
      )
    }

    scenarios
  })

  unlist(out, recursive = FALSE)
}
