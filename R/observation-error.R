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
#'   }
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

  # Helper: compute sigma and rho from a residual vector (with NAs)
  extract_stats <- function(r) {
    r <- r[is.finite(r)]
    n <- length(r)
    if (n < 3) {
      return(list(sigma = NA_real_, rho = 0, n_obs = n))
    }
    s <- sd(r)
    # Lag-1 autocorrelation (sample)
    rho_val <- if (n > 3) {
      stats::cor(r[-n], r[-1])
    } else {
      0
    }
    if (!is.finite(rho_val)) rho_val <- 0
    list(sigma = s, rho = rho_val, n_obs = n)
  }

  if (is.array(resid) && length(dim(resid)) == 3) {
    # Multi-index: dim = [year, area, label]
    labels <- dimnames(resid)[[3]]
    # Pool across areas for each label
    sigmas <- numeric(length(labels))
    rhos <- numeric(length(labels))
    n_obs <- integer(length(labels))
    for (i in seq_along(labels)) {
      r_vec <- as.vector(resid[, , i])
      stats_i <- extract_stats(r_vec)
      sigmas[i] <- stats_i$sigma
      rhos[i] <- stats_i$rho
      n_obs[i] <- stats_i$n_obs
    }
    names(sigmas) <- labels
    names(rhos) <- labels
    names(n_obs) <- labels
  } else {
    # Single-index: vector or matrix [year x area]
    r_vec <- as.vector(resid)
    stats_all <- extract_stats(r_vec)
    sigmas <- stats_all$sigma
    rhos <- stats_all$rho
    n_obs <- stats_all$n_obs
    labels <- NULL
  }

  structure(
    list(
      sigma  = sigmas,
      rho    = rhos,
      labels = labels,
      n_obs  = n_obs
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
  out_list <- vector("list", n_labels)

  for (il in seq_len(n_labels)) {
    sigma <- sigma_vec[il]
    rho <- rho_vec[il]

    for (ia in seq_len(n_areas)) {
      # Generate AR(1) residuals
      eps <- numeric(n_years)
      innovation_sd <- sigma * sqrt(max(0, 1 - rho^2))
      eta <- rnorm(n_years, mean = 0, sd = innovation_sd)
      eps[1] <- eta[1]
      if (n_years > 1) {
        for (t in 2:n_years) {
          eps[t] <- rho * eps[t - 1] + eta[t]
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

      out_list[[length(out_list) + 1]] <- df_chunk
    }
  }

  result <- do.call(rbind, out_list)
  # Remove NA rows (missing observations)
  result <- result[!is.na(result$cpue), , drop = FALSE]
  rownames(result) <- NULL
  result
}
