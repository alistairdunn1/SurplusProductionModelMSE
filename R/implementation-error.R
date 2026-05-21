# Implementation Error (Module 2.2)
#
# Apply implementation error to TAC to get realized catch.
# Catch != TAC due to lognormal error, bias, and overage constraints.


#' Apply Implementation Error to TAC
#'
#' Simulate realized catch from an intended TAC by applying lognormal
#' implementation error with optional bias, temporal autocorrelation,
#' and overage constraints.
#'
#' @param tac Numeric scalar or vector (multi-area). The intended TAC.
#'   Must be non-negative.
#' @param impl_config An \code{\link{impl_error}} configuration object,
#'   or \code{NULL} for perfect implementation (catch = TAC).
#' @param previous_eps Numeric scalar or vector matching \code{tac}.
#'   The error state from the previous year, used for AR(1) continuity.
#'   Default 0 (no history).
#' @param seed Integer random seed for reproducibility, or \code{NULL}.
#'
#' @return A list with:
#'   \describe{
#'     \item{catch}{Numeric vector of realized catch (same length as
#'       \code{tac}).}
#'     \item{eps}{Numeric vector of error states to pass as
#'       \code{previous_eps} in the next year.}
#'   }
#'
#' @details
#' The realised catch is generated as:
#' \deqn{C_t = b \cdot \mathrm{TAC}_t \cdot \exp(\varepsilon_t - \sigma^2 / 2)}
#' where \eqn{b} is the multiplicative bias,
#' \eqn{\sigma = \sqrt{\log(1 + \mathrm{cv}^2)}} converts the coefficient
#' of variation to a log-scale standard deviation, and
#' \eqn{\varepsilon_t = \rho \cdot \varepsilon_{t-1} +
#' \sqrt{1 - \rho^2} \cdot \eta_t} with
#' \eqn{\eta_t \sim N(0, \sigma^2)}.
#'
#' The bias-correction term \eqn{-\sigma^2/2} ensures that for
#' \eqn{b = 1} (unbiased), the expected catch equals the TAC.
#'
#' After applying the random error, catch is capped at
#' \code{max_overage * tac} and floored at zero. When TAC is zero
#' (fishery closure), catch is zero and the error state resets.
#'
#' @examples
#' \dontrun{
#' cfg <- impl_error(cv = 0.1, bias = 1, autocorr = 0, max_overage = 1.1)
#' result <- apply_implementation_error(1000, cfg, seed = 42)
#' result$catch
#'
#' # Chain across years for AR(1)
#' r1 <- apply_implementation_error(1000, cfg, previous_eps = 0)
#' r2 <- apply_implementation_error(1200, cfg, previous_eps = r1$eps)
#' }
#'
#' @export
apply_implementation_error <- function(tac,
                                       impl_config,
                                       previous_eps = 0,
                                       seed = NULL) {
  # Perfect implementation
  if (is.null(impl_config)) {
    return(list(catch = tac, eps = rep(0, length(tac))))
  }

  if (!inherits(impl_config, "impl_error")) {
    stop("impl_config must be an 'impl_error' object or NULL", call. = FALSE)
  }

  assert_numeric(tac,
    lower = 0, any.missing = FALSE, min.len = 1,
    .var.name = "tac"
  )

  if (!is.null(seed)) set.seed(seed)

  n <- length(tac)
  previous_eps <- rep_len(previous_eps, n)

  # Identify closures (zero TAC)
  is_zero <- tac <= 0

  # Extract configuration
  cv <- impl_config$cv
  sigma <- sqrt(log(1 + cv^2))
  rho <- impl_config$autocorr
  bias <- impl_config$bias
  max_overage <- impl_config$max_overage

  # AR(1) lognormal error
  innovation_sd <- sigma * sqrt(max(0, 1 - rho^2))
  eta <- rnorm(n, mean = 0, sd = innovation_sd)
  eps <- rho * previous_eps + eta

  # Bias-corrected lognormal realisation
  catch_val <- bias * tac * exp(eps - sigma^2 / 2)

  # Enforce overage cap
  catch_val <- pmin(catch_val, max_overage * tac)

  # Floor at zero
  catch_val <- pmax(catch_val, 0)

  # Zero TAC => zero catch, reset error state
  catch_val[is_zero] <- 0
  eps[is_zero] <- 0

  list(catch = catch_val, eps = eps)
}
