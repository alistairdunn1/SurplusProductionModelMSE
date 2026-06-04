# Spatial Robustness Testing (Module 2.4)
#
# Compare MSE scenario performance and validate the simulation framework
# via self-tests and spatial misspecification diagnostics.


#' Compare Scenario Performance
#'
#' Extract performance metrics from one or more \code{mse_result} objects
#' and produce a trade-off plot (Pareto frontier) between two selected
#' metrics.
#'
#' @param mse_results An \code{mse_result} object (from
#'   \code{\link{mse_simulation}}), or a named list of such objects.
#'   When a list is supplied, scenarios from all elements are pooled.
#' @param metric_x Character. Performance metric for the x-axis (see
#'   Details).
#' @param metric_y Character. Performance metric for the y-axis.
#' @param scope Character. Spatial scope to extract: \code{"aggregate"}
#'   (default) or an area name (e.g. \code{"A1"}).
#'
#' @return A \code{ggplot} object showing one point per scenario,
#'   labelled by name, with Pareto-optimal scenarios highlighted.
#'
#' @details
#' Metric names correspond to the \code{metric} column from
#' \code{\link{summary.mse_performance}}.  Common choices include:
#' \describe{
#'   \item{Catch metrics}{\code{mean_catch}, \code{AAV}}
#'   \item{Biomass risk}{\code{Pr(B<50\%B0)_final},
#'     \code{Pr(B<50\%B0)_ever}, \code{Pr(B<20\%B0)_final},
#'     \code{Pr(B<20\%B0)_ever}}
#'   \item{F risk}{\code{Pr(F>F50\%B0)_final},
#'     \code{Pr(F>F50\%B0)_ever}}
#'   \item{Biomass ratios}{\code{mean_B_B0}, \code{mean_B_BMSY},
#'     \code{final_B_B0}, \code{final_B_BMSY}}
#' }
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
#' s1 <- create_scenario("low_f", hcr_constant_f(0.03))
#' s2 <- create_scenario("high_f", hcr_constant_f(0.10))
#' res <- mse_simulation(om,
#'   scenarios = list(s1, s2),
#'   n_sims = 20, n_proj_years = 15, seed = 1
#' )
#' compare_scenarios(res, "mean_catch", "Pr(B<20%B0)_final")
#' }
#'
#' @export
compare_scenarios <- function(mse_results,
                              metric_x,
                              metric_y,
                              scope = "aggregate") {
  assert_character(metric_x, len = 1, .var.name = "metric_x")
  assert_character(metric_y, len = 1, .var.name = "metric_y")
  assert_character(scope, len = 1, .var.name = "scope")

  # Collect per-scenario summary data frames
  df <- .extract_scenario_metrics(mse_results, scope)

  # Pivot to wide-ish form with one row per scenario
  vals_x <- df[df$metric == metric_x, , drop = FALSE]
  vals_y <- df[df$metric == metric_y, , drop = FALSE]

  if (nrow(vals_x) == 0) {
    stop("metric_x '", metric_x, "' not found. Available: ",
      paste(unique(df$metric), collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(vals_y) == 0) {
    stop("metric_y '", metric_y, "' not found. Available: ",
      paste(unique(df$metric), collapse = ", "),
      call. = FALSE
    )
  }

  plot_df <- merge(
    vals_x[, c("scenario", "value")],
    vals_y[, c("scenario", "value")],
    by = "scenario",
    suffixes = c("_x", "_y")
  )

  # Identify Pareto-optimal points (non-dominated in BOTH dimensions)
  plot_df$pareto <- .is_pareto_optimal(
    plot_df$value_x, plot_df$value_y,
    metric_x, metric_y
  )

  ggplot(plot_df, aes(
    x = .data$value_x, y = .data$value_y,
    label = .data$scenario
  )) +
    geom_point(aes(colour = .data$pareto), size = 3) +
    geom_text(
      nudge_y = diff(range(plot_df$value_y, na.rm = TRUE)) * 0.04,
      size = 3
    ) +
    scale_colour_manual(
      values = c("FALSE" = "grey50", "TRUE" = "firebrick"),
      name = "Pareto optimal"
    ) +
    labs(x = metric_x, y = metric_y) +
    theme_bw()
}


#' Run Self-Test Validation
#'
#' Execute a self-test where the estimation model matches the operating
#' model structure, confirming that the MSE framework recovers true
#' dynamics. This is a diagnostic: under perfect structural match the
#' HCR should perform near-optimally.
#'
#' @param operating_model An \code{\link{om_config}} object with
#'   \code{true_params} set.
#' @param harvest_control_rule A harvest control rule function (from
#'   \code{\link{hcr_constant_f}} or similar).
#' @param n_sims Integer. Number of replicates (default 20).
#' @param n_proj_years Integer. Projection years (default 20).
#' @param seed Integer random seed, or \code{NULL}.
#' @param ... Additional arguments passed to
#'   \code{\link{mse_simulation}}.
#'
#' @return An S3 object of class \code{self_test_result} containing:
#'   \describe{
#'     \item{mse_result}{The underlying \code{mse_result} object.}
#'     \item{biomass_conserved}{Logical. TRUE if total biomass at
#'       year 1 matches B0 across all replicates.}
#'     \item{positive_biomass}{Logical. TRUE if all biomass values
#'       are positive.}
#'     \item{summary}{A data frame of aggregate performance metrics.}
#'   }
#'
#' @details
#' The self-test passes \code{estimation_model = NULL} to
#' \code{mse_simulation()}, meaning the EM inherits the OM structure.
#' The returned diagnostics help verify that:
#' \enumerate{
#'   \item Total biomass is conserved at initialisation.
#'   \item Biomass remains positive throughout the projection.
#'   \item Performance metrics are sensible (depletion stays
#'     reasonable under a moderate HCR).
#' }
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
#' st <- run_self_test(om, hcr_constant_f(0.05), n_sims = 10, seed = 1)
#' st
#' }
#'
#' @export
run_self_test <- function(operating_model,
                          harvest_control_rule,
                          n_sims = 20L,
                          n_proj_years = 20L,
                          seed = NULL,
                          ...) {
  if (!inherits(operating_model, "om_config")) {
    stop("operating_model must be an 'om_config' object", call. = FALSE)
  }
  assert_function(harvest_control_rule, .var.name = "harvest_control_rule")

  sc <- create_scenario(
    name = "self_test",
    harvest_control_rule = harvest_control_rule
  )

  res <- mse_simulation(
    operating_model  = operating_model,
    estimation_model = NULL,
    scenarios        = sc,
    n_sims           = n_sims,
    n_proj_years     = n_proj_years,
    seed             = seed,
    ...
  )

  traj <- res$results$self_test$trajectories

  # Check biomass conservation at initialisation: year-1 biomass must equal the
  # initial conditions (sum of tp$B_initial, which may differ from K if the stock is
  # not unfished at the start of the projection).
  b_initial_total <- sum(rep_len(
    operating_model$true_params$B_initial,
    operating_model$n_areas
  ))
  total_b_yr1 <- apply(traj$biomass[, 1, , drop = FALSE], 1, sum)
  biomass_conserved <- all(abs(total_b_yr1 - b_initial_total) < 1e-6)

  # Check positive biomass
  positive_biomass <- all(traj$biomass > 0)

  # Performance summary
  perf_summary <- summary(res$results$self_test$performance)

  structure(
    list(
      mse_result        = res,
      biomass_conserved = biomass_conserved,
      positive_biomass  = positive_biomass,
      summary           = perf_summary
    ),
    class = "self_test_result"
  )
}


#' @export
print.self_test_result <- function(x, ...) {
  cat("MSE Self-Test Result\n")
  cat("====================\n")
  cat("Biomass conserved at init:", x$biomass_conserved, "\n")
  cat("All biomass positive:     ", x$positive_biomass, "\n\n")

  cat("Aggregate Performance Summary:\n")
  s <- x$summary[x$summary$scope == "aggregate", ]
  for (i in seq_len(nrow(s))) {
    cat("  ", s$metric[i], "=", round(s$value[i], 4), "\n")
  }
  invisible(x)
}


# ========================================================================
# Internal helpers
# ========================================================================

#' Extract per-scenario metrics from mse_result(s)
#' @noRd
.extract_scenario_metrics <- function(mse_results, scope) {
  if (inherits(mse_results, "mse_result")) {
    mse_list <- list(mse_results)
  } else if (is.list(mse_results)) {
    mse_list <- mse_results
  } else {
    stop("mse_results must be an 'mse_result' object or a list of them",
      call. = FALSE
    )
  }

  all_rows <- list()

  for (res in mse_list) {
    if (!inherits(res, "mse_result")) {
      stop("Each element of mse_results must be an 'mse_result' object",
        call. = FALSE
      )
    }
    for (sc_name in names(res$results)) {
      perf <- res$results[[sc_name]]$performance
      s <- summary(perf)
      s_scope <- s[s$scope == scope, , drop = FALSE]
      s_scope$scenario <- sc_name
      all_rows[[length(all_rows) + 1]] <- s_scope
    }
  }

  do.call(rbind, all_rows)
}


#' Determine Pareto optimality
#'
#' A point is Pareto-optimal if no other point is at least as good in
#' both metrics and strictly better in at least one. The "direction"
#' depends on the metric: risk and AAV are better when lower; catch
#' and biomass ratios are better when higher.
#' @noRd
.is_pareto_optimal <- function(x, y, metric_x, metric_y) {
  n <- length(x)
  if (n <= 1) {
    return(rep(TRUE, n))
  }

  # Determine preferred direction (lower is better for risk/AAV)
  dir_x <- .metric_direction(metric_x) # +1 = higher is better, -1 = lower

  dir_y <- .metric_direction(metric_y)

  # Transform so that higher is always better
  xx <- dir_x * x
  yy <- dir_y * y

  pareto <- logical(n)
  for (i in seq_len(n)) {
    dominated <- FALSE
    for (j in seq_len(n)) {
      if (j == i) next
      if (xx[j] >= xx[i] && yy[j] >= yy[i] &&
        (xx[j] > xx[i] || yy[j] > yy[i])) {
        dominated <- TRUE
        break
      }
    }
    pareto[i] <- !dominated
  }
  pareto
}


#' Determine whether higher or lower is better for a metric
#'
#' Returns +1 if higher is better, -1 if lower is better.
#' @noRd
.metric_direction <- function(metric_name) {
  # Risk metrics and AAV: lower is better

  lower_patterns <- c("^Pr\\(", "^AAV$")
  for (pat in lower_patterns) {
    if (grepl(pat, metric_name)) {
      return(-1)
    }
  }
  # Everything else (catch, biomass ratios): higher is better
  1
}
