# Operating Model Projection (Module 2.3)
#
# Forward-project true biomass using the Pella-Tomlinson production model,
# with optional process noise and spatial movement.


#' Build Gravity Movement Kernel
#'
#' Compute a row-normalized transition matrix for spatial biomass
#' redistribution using a gravity model.
#'
#' @param distance_matrix Numeric matrix \eqn{[n \times n]}. Pairwise
#'   distances between areas. Must be symmetric with zero diagonal.
#' @param attractiveness Numeric vector of length \eqn{n}. Habitat quality
#'   weights per area. Default \code{NULL} uses equal weights.
#' @param decay Numeric scalar \eqn{>= 0}. Distance decay parameter.
#'
#' @return A row-normalized \eqn{n \times n} transition matrix where entry
#'   \eqn{(a, b)} is the fraction of biomass in area \eqn{a} that is
#'   attracted to area \eqn{b}.
#'
#' @details
#' Movement weights are:
#' \deqn{W_{ab} = A_b \cdot \exp(-\delta \cdot D_{ab})}
#'
#' Each row is then normalized to sum to 1:
#' \deqn{K_{ab} = W_{ab} / \sum_b W_{ab}}
#'
#' @export
build_movement_kernel <- function(distance_matrix,
                                  attractiveness = NULL,
                                  decay = 0) {
  n <- nrow(distance_matrix)
  if (is.null(attractiveness)) {
    attractiveness <- rep(1, n)
  }

  # Weight matrix: W[a,b] = A[b] * exp(-decay * D[a,b])
  W <- matrix(0, n, n)
  for (a in seq_len(n)) {
    for (b in seq_len(n)) {
      W[a, b] <- attractiveness[b] * exp(-decay * distance_matrix[a, b])
    }
  }

  # Row-normalize
  row_sums <- rowSums(W)
  K <- W / row_sums

  K
}


#' Project Biomass One Year Forward
#'
#' Advance the operating model biomass state by one time step using the
#' Pella-Tomlinson surplus production model, with optional process noise
#' and spatial movement.
#'
#' @param biomass Numeric scalar or vector. Current biomass per area.
#' @param catch Numeric scalar or vector. Realized catch per area (same
#'   length as \code{biomass}).
#' @param om_config An \code{\link{om_config}} object containing true
#'   population parameters.
#' @param movement_kernel Optional pre-computed movement kernel matrix from
#'   \code{\link{build_movement_kernel}}. If \code{NULL} and
#'   \code{om_config$n_areas > 1}, it will be computed from the config.
#' @param process_noise Logical. Whether to add stochastic process noise.
#'   Default \code{TRUE} if \code{sigma_process} is present in
#'   \code{true_params}.
#' @param seed Integer random seed, or \code{NULL}.
#'
#' @return A numeric vector of biomass per area at \eqn{t + 1}.
#'
#' @details
#' The deterministic Pella-Tomlinson update is:
#' \deqn{B_{t+1,a} = B_{t,a} + P(B_{t,a}) - C_{t,a}}
#' where
#' \deqn{P(B) = \frac{r \cdot B \cdot \bigl(1 - (B/K_a)^{m-1}\bigr)}{m}}
#'
#' For multi-area models, \eqn{K} is distributed across areas proportional
#' to initial biomass \eqn{B_0}.
#'
#' If \code{sigma_process > 0}, log-scale process noise is added:
#' \deqn{B_{t+1,a} = (B_{t,a} + P(B_{t,a}) - C_{t,a}) \cdot
#'   \exp(\varepsilon_a - \sigma_p^2/2)}
#' where \eqn{\varepsilon_a \sim N(0, \sigma_p^2)}.
#'
#' Spatial redistribution (when \code{movement_rate > 0}):
#' \deqn{B'_a = (1 - \rho) \cdot B_a + \rho \cdot \sum_b K_{ab} \cdot B_b}
#' where \eqn{K} is the gravity movement kernel.
#'
#' Biomass is floored at a small positive value (0.01) to prevent
#' numerical issues.
#'
#' @examples
#' \dontrun{
#' cfg <- om_config(
#'   n_areas = 1,
#'   true_params = list(
#'     r = 0.3, K = 5000, m = 2,
#'     sigma_obs = 0.2, q = 1e-4, B0 = 5000
#'   )
#' )
#' B1 <- project_biomass(5000, 200, cfg, seed = 1)
#' }
#'
#' @export
project_biomass <- function(biomass,
                            catch,
                            om_config,
                            movement_kernel = NULL,
                            process_noise = NULL,
                            seed = NULL) {
  if (!inherits(om_config, "om_config")) {
    stop("om_config must be an 'om_config' object", call. = FALSE)
  }

  tp <- om_config$true_params
  if (is.null(tp)) {
    stop("om_config$true_params must be set for projection", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  n_areas <- om_config$n_areas
  biomass <- rep_len(as.numeric(biomass), n_areas)
  catch <- rep_len(as.numeric(catch), n_areas)

  r <- tp$r
  K <- tp$K
  m <- tp$m

  # Distribute K across areas proportional to B0
  B0 <- rep_len(tp$B0, n_areas)
  if (n_areas > 1) {
    K_area <- K * (B0 / sum(B0))
  } else {
    K_area <- K
  }

  # Pella-Tomlinson production per area
  production <- numeric(n_areas)
  for (a in seq_len(n_areas)) {
    if (biomass[a] > 0 && K_area[a] > 0) {
      production[a] <- r * biomass[a] *
        (1 - (biomass[a] / K_area[a])^(m - 1)) / m
    }
  }

  # State update: B + P - C
  B_new <- biomass + production - catch

  # Process noise
  sigma_p <- if (!is.null(tp$sigma_process)) tp$sigma_process else 0
  if (is.null(process_noise)) {
    process_noise <- sigma_p > 0
  }

  if (process_noise && sigma_p > 0) {
    eps <- rnorm(n_areas, mean = 0, sd = sigma_p)
    # Bias-corrected lognormal noise on positive biomass
    B_positive <- pmax(B_new, 1e-8)
    B_new <- B_positive * exp(eps - sigma_p^2 / 2)
  }

  # Floor at small positive value
  B_new <- pmax(B_new, 0.01)

  # Spatial movement
  if (n_areas > 1 && om_config$movement_rate > 0) {
    if (is.null(movement_kernel)) {
      dm <- om_config$distance_matrix
      if (is.null(dm)) {
        dm <- matrix(1, n_areas, n_areas)
        diag(dm) <- 0
      }
      movement_kernel <- build_movement_kernel(
        dm, om_config$attractiveness, om_config$decay
      )
    }
    rho <- om_config$movement_rate
    # Transpose: K[a,b] = fraction FROM a TO b (row-stochastic)
    # so t(K) %*% B gives biomass ARRIVING in each area => conserves total
    B_moved <- (1 - rho) * B_new + rho * as.vector(t(movement_kernel) %*% B_new)
    B_new <- B_moved
  }

  # Final floor
  B_new <- pmax(B_new, 0.01)

  B_new
}


#' Project Biomass Trajectory
#'
#' Run a multi-year forward projection of the operating model.
#'
#' @param B0 Numeric scalar or vector. Initial biomass per area.
#' @param catch_series Numeric vector or matrix. Catch time series. If a
#'   vector, used for all areas (single-area) or recycled. If a matrix,
#'   rows are years and columns are areas.
#' @param om_config An \code{\link{om_config}} object.
#' @param n_years Integer. Number of years to project (derived from
#'   \code{catch_series} if not specified).
#' @param process_noise Logical. Whether to add process noise. Default
#'   \code{TRUE} if \code{sigma_process} is in \code{true_params}.
#' @param seed Integer random seed, or \code{NULL}.
#'
#' @return A matrix of biomass with \code{n_years + 1} rows (including
#'   the initial year) and \code{n_areas} columns. Row 1 is the initial
#'   biomass.
#'
#' @examples
#' \dontrun{
#' cfg <- om_config(
#'   n_areas = 1,
#'   true_params = list(
#'     r = 0.3, K = 5000, m = 2,
#'     sigma_obs = 0.2, q = 1e-4, B0 = 5000
#'   )
#' )
#' catches <- rep(200, 10)
#' traj <- project_trajectory(5000, catches, cfg, seed = 42)
#' }
#'
#' @export
project_trajectory <- function(B0,
                               catch_series,
                               om_config,
                               n_years = NULL,
                               process_noise = NULL,
                               seed = NULL) {
  if (!inherits(om_config, "om_config")) {
    stop("om_config must be an 'om_config' object", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  n_areas <- om_config$n_areas
  B0 <- rep_len(as.numeric(B0), n_areas)


  # Coerce catch to matrix [n_years x n_areas]
  if (is.vector(catch_series) || (is.matrix(catch_series) && ncol(catch_series) == 1)) {
    catch_mat <- matrix(as.numeric(catch_series), ncol = 1)
    if (n_areas > 1) {
      catch_mat <- matrix(rep(catch_mat[, 1], n_areas),
        ncol = n_areas
      )
    }
  } else {
    catch_mat <- as.matrix(catch_series)
  }

  if (is.null(n_years)) {
    n_years <- nrow(catch_mat)
  }
  if (nrow(catch_mat) < n_years) {
    stop("catch_series has fewer rows (", nrow(catch_mat),
      ") than n_years (", n_years, ")",
      call. = FALSE
    )
  }

  # Pre-compute movement kernel
  movement_kernel <- NULL
  if (n_areas > 1 && om_config$movement_rate > 0) {
    dm <- om_config$distance_matrix
    if (is.null(dm)) {
      dm <- matrix(1, n_areas, n_areas)
      diag(dm) <- 0
    }
    movement_kernel <- build_movement_kernel(
      dm, om_config$attractiveness, om_config$decay
    )
  }

  # Trajectory matrix: (n_years + 1) rows x n_areas cols
  traj <- matrix(NA_real_, nrow = n_years + 1, ncol = n_areas)
  traj[1, ] <- B0

  for (t in seq_len(n_years)) {
    traj[t + 1, ] <- project_biomass(
      biomass = traj[t, ],
      catch = catch_mat[t, ],
      om_config = om_config,
      movement_kernel = movement_kernel,
      process_noise = process_noise
    )
  }

  if (n_areas == 1) {
    colnames(traj) <- "biomass"
  } else {
    colnames(traj) <- paste0("A", seq_len(n_areas))
  }

  traj
}
