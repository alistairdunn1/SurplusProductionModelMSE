# Observation Error Model (Module 2.1)
#
# Calibrate observation error from fitted model residuals and simulate
# future CPUE with realistic error structure.


#' Calibrate Observation Error from a Fitted Model
#'
#' Extract residual statistics from a fitted \code{ProductionModel} to
#' parameterise observation error for MSE simulations.
#'
#' @param model_fit A fitted \code{ProductionModel} object (from
#'   \code{\link[SurplusProductionModel]{fit_pella_tomlinson_model}}).
#'
#' @return A list of class \code{obs_error_params} containing:
#'   \describe{
#'     \item{sigma}{Numeric (scalar or per-index vector). Residual standard
#'       deviation on the log scale.}
#'     \item{rho}{Numeric (scalar or per-index vector). Lag-1 autocorrelation
#'       of residuals.}
#'     \item{labels}{Character vector of CPUE index labels, or \code{NULL}
#'       for single-index models.}
#'     \item{n_obs}{Integer (scalar or per-index vector). Number of
#'       non-missing observations used.}
#'     \item{by_series}{Data frame of the within-area (and, where relevant,
#'       within-index) estimates used to construct the aggregate parameters.}
#'   }
#'
#' @details
#' Temporal autocorrelation is estimated separately for each area/index time
#' series, retaining only consecutive finite residual pairs. Area-specific
#' residual variances are pooled using their residual degrees of freedom.
#' Lag-one correlations are pooled on the Fisher z scale, weighted by the
#' number of consecutive residual pairs. This avoids treating transitions
#' between areas as temporal observations.
#'
#' @examples
#' \dontrun{
#' fit <- fit_pella_tomlinson_model(cpue, catch)
#' obs_params <- calibrate_observation_error(fit)
#' obs_params$sigma
#' }
#'
#' @export
calibrate_observation_error <- function(model_fit) {
  if (!inherits(model_fit, "ProductionModel")) {
    stop("model_fit must be a 'ProductionModel' object", call. = FALSE)
  }
  if (!model_fit$fitted) {
    stop("model_fit must be a fitted model (run fit_pella_tomlinson_model first)",
      call. = FALSE
    )
  }

  resid <- model_fit$results$residuals

  # Estimate each temporal series independently so area transitions are not
  # interpreted as successive residuals.
  extract_series_stats <- function(r) {
    finite <- is.finite(r)
    n_obs <- sum(finite)
    sigma <- if (n_obs >= 2) stats::sd(r[finite]) else NA_real_

    paired <- if (length(r) >= 2) {
      finite[-length(r)] & finite[-1]
    } else {
      logical()
    }
    n_pairs <- sum(paired)
    rho <- if (n_pairs >= 2) {
      stats::cor(r[-length(r)][paired], r[-1][paired])
    } else {
      0
    }
    if (!is.finite(rho)) rho <- 0

    list(sigma = sigma, rho = rho, n_obs = n_obs, n_pairs = n_pairs)
  }

  combine_series_stats <- function(series_stats) {
    n_obs <- vapply(series_stats, `[[`, integer(1), "n_obs")
    n_pairs <- vapply(series_stats, `[[`, integer(1), "n_pairs")
    sigma <- vapply(series_stats, `[[`, numeric(1), "sigma")
    rho <- vapply(series_stats, `[[`, numeric(1), "rho")

    variance_weights <- pmax(n_obs - 1, 0)
    valid_sigma <- is.finite(sigma) & variance_weights > 0
    pooled_sigma <- if (any(valid_sigma)) {
      sqrt(sum(variance_weights[valid_sigma] * sigma[valid_sigma]^2) /
        sum(variance_weights[valid_sigma]))
    } else {
      NA_real_
    }

    valid_rho <- is.finite(rho) & n_pairs > 0
    pooled_rho <- if (any(valid_rho)) {
      rho_bounded <- pmin(pmax(rho[valid_rho], -0.999999), 0.999999)
      tanh(stats::weighted.mean(atanh(rho_bounded), n_pairs[valid_rho]))
    } else {
      0
    }

    list(
      sigma = pooled_sigma,
      rho = pooled_rho,
      n_obs = sum(n_obs),
      n_pairs = sum(n_pairs)
    )
  }

  make_series_table <- function(series_stats, areas, label = NULL) {
    label_column <- if (is.null(label)) {
      rep(NA_character_, length(areas))
    } else {
      rep(label, length(areas))
    }
    data.frame(
      area = areas,
      label = label_column,
      sigma = vapply(series_stats, `[[`, numeric(1), "sigma"),
      rho = vapply(series_stats, `[[`, numeric(1), "rho"),
      n_obs = vapply(series_stats, `[[`, integer(1), "n_obs"),
      n_pairs = vapply(series_stats, `[[`, integer(1), "n_pairs"),
      row.names = NULL
    )
  }

  if (is.array(resid) && length(dim(resid)) == 3) {
    # Multi-index: dim = [year, area, label]
    labels <- dimnames(resid)[[3]]
    areas <- dimnames(resid)[[2]] %||% paste0("A", seq_len(dim(resid)[2]))
    # Pool area-specific estimates for each index label.
    sigmas <- numeric(length(labels))
    rhos <- numeric(length(labels))
    n_obs <- integer(length(labels))
    series_tables <- vector("list", length(labels))
    for (i in seq_along(labels)) {
      series_stats <- lapply(seq_along(areas), function(area) {
        extract_series_stats(resid[, area, i])
      })
      stats_i <- combine_series_stats(series_stats)
      sigmas[i] <- stats_i$sigma
      rhos[i] <- stats_i$rho
      n_obs[i] <- stats_i$n_obs
      series_tables[[i]] <- make_series_table(series_stats, areas, labels[i])
    }
    names(sigmas) <- labels
    names(rhos) <- labels
    names(n_obs) <- labels
    by_series <- do.call(rbind, series_tables)
  } else {
    # Single-index: vector or matrix [year x area].
    r_matrix <- as.matrix(resid)
    areas <- colnames(r_matrix) %||% paste0("A", seq_len(ncol(r_matrix)))
    series_stats <- lapply(seq_along(areas), function(area) {
      extract_series_stats(r_matrix[, area])
    })
    stats_all <- combine_series_stats(series_stats)
    sigmas <- stats_all$sigma
    rhos <- stats_all$rho
    n_obs <- stats_all$n_obs
    labels <- NULL
    by_series <- make_series_table(series_stats, areas)
  }

  structure(
    list(
      sigma  = sigmas,
      rho    = rhos,
      labels = labels,
      n_obs  = n_obs,
      by_series = by_series
    ),
    class = "obs_error_params"
  )
}


#' @export
print.obs_error_params <- function(x, ...) {
  cat("Observation Error Parameters\n")
  cat("----------------------------\n")
  if (!is.null(x$labels)) {
    for (i in seq_along(x$labels)) {
      cat(
        "  ", x$labels[i], ": sigma =", round(x$sigma[i], 4),
        " rho =", round(x$rho[i], 3),
        " (n =", x$n_obs[i], ")\n"
      )
    }
  } else {
    cat(
      "  sigma =", round(x$sigma, 4),
      " rho =", round(x$rho, 3),
      " (n =", x$n_obs, ")\n"
    )
  }
  invisible(x)
}


#' Simulate Future CPUE
#'
#' Generate simulated CPUE observations from true biomass using lognormal
#' observation error, optionally with temporal autocorrelation and missing
#' data gaps.
#'
#' @param true_biomass Numeric vector or matrix. True biomass from the
#'   operating model. If a matrix, rows are years and columns are areas.
#' @param obs_error_params An \code{obs_error_params} object from
#'   \code{\link{calibrate_observation_error}}, or a list with at least
#'   \code{sigma} (scalar or per-index).
#' @param q Numeric. Catchability coefficient (scalar or per-area vector).
#' @param years Integer vector of years corresponding to \code{true_biomass}
#'   rows.
#' @param labels Character vector of CPUE index labels to simulate, or
#'   \code{NULL} for a single unnamed index (default).
#' @param missing_prob Numeric in \eqn{[0, 1)}. Probability that any
#'   year/area/label observation is missing (default 0 = no gaps).
#' @param seed Integer. Random seed for reproducibility, or \code{NULL}.
#'
#' @return A data frame with columns \code{year}, \code{cpue}, and optionally
#'   \code{area} and \code{label}, matching the format expected by
#'   \code{\link[SurplusProductionModel]{fit_pella_tomlinson_model}}.
#'
#' @details
#' CPUE is generated as:
#' \deqn{I_{t} = q \cdot B_{t} \cdot \exp(\varepsilon_t - \sigma^2 / 2)}
#' where \eqn{\varepsilon_t = \rho \cdot \varepsilon_{t-1} + \sqrt{1 - \rho^2} \cdot \eta_t}
#' and \eqn{\eta_t \sim N(0, \sigma^2)}.
#'
#' The bias-correction term \eqn{-\sigma^2/2} ensures \eqn{E[I_t] = q \cdot B_t}.
#'
#' @examples
#' \dontrun{
#' # Simulate from known biomass
#' B <- c(4000, 3800, 3600, 3500, 3400)
#' cpue_sim <- simulate_cpue(B, list(sigma = 0.2, rho = 0),
#'   q = 1e-4, years = 2020:2024
#' )
#' }
#'
#' @export
simulate_cpue <- function(true_biomass,
                          obs_error_params,
                          q,
                          years,
                          labels = NULL,
                          missing_prob = 0,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Coerce biomass to matrix [n_years x n_areas]
  if (is.vector(true_biomass) || (is.matrix(true_biomass) && ncol(true_biomass) == 1)) {
    B <- matrix(as.numeric(true_biomass), ncol = 1)
    area_names <- "A1"
  } else {
    B <- as.matrix(true_biomass)
    area_names <- if (!is.null(colnames(B))) {
      colnames(B)
    } else {
      paste0("A", seq_len(ncol(B)))
    }
  }
  n_years <- nrow(B)
  n_areas <- ncol(B)

  assert_integerish(years,
    len = n_years, any.missing = FALSE,
    .var.name = "years"
  )
  assert_numeric(q,
    lower = 0, any.missing = FALSE, min.len = 1,
    max.len = n_areas, .var.name = "q"
  )
  assert_number(missing_prob, lower = 0, upper = 1, .var.name = "missing_prob")

  # Recycle q to n_areas
  q_vec <- rep_len(q, n_areas)

  # Determine sigma/rho per label
  if (inherits(obs_error_params, "obs_error_params")) {
    sigma_all <- obs_error_params$sigma
    rho_all <- obs_error_params$rho
    if (is.null(labels) && !is.null(obs_error_params$labels)) {
      labels <- obs_error_params$labels
    }
  } else {
    assert_list(obs_error_params, .var.name = "obs_error_params")
    sigma_all <- obs_error_params$sigma
    rho_all <- if (!is.null(obs_error_params$rho)) obs_error_params$rho else 0
  }

  if (is.null(labels)) {
    labels_vec <- NULL
    n_labels <- 1L
    sigma_vec <- rep_len(sigma_all, 1)
    rho_vec <- rep_len(rho_all, 1)
  } else {
    labels_vec <- labels
    n_labels <- length(labels)
    sigma_vec <- rep_len(sigma_all, n_labels)
    rho_vec <- rep_len(rho_all, n_labels)
  }

  # Build output data frame
  out_list <- vector("list", n_labels * n_areas)
  out_idx  <- 0L

  for (il in seq_len(n_labels)) {
    sigma <- sigma_vec[il]
    rho <- rho_vec[il]

    for (ia in seq_len(n_areas)) {
      # Generate AR(1) residuals.
      # Year 1 is drawn from the stationary distribution N(0, sigma^2);
      # subsequent years follow the AR(1) recurrence.
      eps <- numeric(n_years)
      innovation_sd <- sigma * sqrt(max(0, 1 - rho^2))
      eps[1] <- rnorm(1L, mean = 0, sd = sigma)
      if (n_years > 1) {
        for (t in 2:n_years) {
          eps[t] <- rho * eps[t - 1] + rnorm(1L, mean = 0, sd = innovation_sd)
        }
      }

      # Bias-corrected lognormal CPUE
      cpue_vals <- q_vec[ia] * B[, ia] * exp(eps - sigma^2 / 2)

      # Apply missing data
      if (missing_prob > 0) {
        is_missing <- runif(n_years) < missing_prob
        cpue_vals[is_missing] <- NA
      }

      df_chunk <- data.frame(
        year = as.integer(years),
        cpue = cpue_vals,
        stringsAsFactors = FALSE
      )

      if (n_areas > 1) {
        df_chunk$area <- area_names[ia]
      }
      if (!is.null(labels_vec)) {
        df_chunk$label <- labels_vec[il]
      }

      out_idx <- out_idx + 1L
      out_list[[out_idx]] <- df_chunk
    }
  }

  result <- do.call(rbind, out_list)
  # Remove NA rows (missing observations)
  result <- result[!is.na(result$cpue), , drop = FALSE]
  rownames(result) <- NULL
  result
}


#' Generate a single-year CPUE observation with AR(1) state
#'
#' Called by the MSE closed-loop to generate one observation per year while
#' propagating the observation-error AR(1) state across years.
#'
#' @param true_biomass_total Numeric scalar. Aggregate (sum across areas) true
#'   biomass before fishing in the current year.
#' @param obs_error_params An \code{obs_error_params} object.
#' @param q Numeric scalar. Effective catchability coefficient.
#' @param expected_index Optional positive scalar giving the deterministic
#'   expected index before observation error. When supplied, it replaces
#'   \code{q * true_biomass_total}; this supports spatially weighted aggregate
#'   indices in the MSE loop.
#' @param year Integer. Simulation year (stored in the returned data frame).
#' @param previous_eps Numeric scalar. AR(1) state from the previous year, or
#'   \code{NULL} for the first call (draws from the stationary distribution).
#'
#' @return A list with:
#'   \describe{
#'     \item{cpue_row}{Single-row data frame with columns \code{year} and
#'       \code{cpue}.}
#'     \item{eps}{Updated AR(1) state to pass as \code{previous_eps} next year.}
#'   }
#'
#' @noRd
.simulate_cpue_step <- function(true_biomass_total, obs_error_params, q,
                                year, previous_eps = NULL,
                                expected_index = NULL) {
  sigma <- as.numeric(obs_error_params$sigma[[1]])
  rho   <- as.numeric(obs_error_params$rho[[1]])

  if (is.null(previous_eps)) {
    # Initialise from the stationary marginal distribution N(0, sigma^2).
    eps_new <- rnorm(1L, mean = 0, sd = sigma)
  } else {
    innovation_sd <- sigma * sqrt(max(0, 1 - rho^2))
    eps_new <- rho * as.numeric(previous_eps) +
      rnorm(1L, mean = 0, sd = innovation_sd)
  }

  index_mean <- if (is.null(expected_index)) {
    q * true_biomass_total
  } else {
    as.numeric(expected_index)
  }
  if (length(index_mean) != 1 || !is.finite(index_mean) || index_mean <= 0) {
    stop("expected index must be a positive finite scalar", call. = FALSE)
  }
  cpue_val <- index_mean * exp(eps_new - sigma^2 / 2)

  list(
    cpue_row = data.frame(
      year = as.integer(year),
      cpue = cpue_val,
      stringsAsFactors = FALSE
    ),
    eps = eps_new
  )
}
