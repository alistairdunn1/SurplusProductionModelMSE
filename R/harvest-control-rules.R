# Harvest Control Rule factory functions
#
# Each function returns a closure with signature function(biomass, reference_points)
# that returns a TAC (Total Allowable Catch) value.


#' Constant Exploitation-Rate HCR
#'
#' Returns a harvest control rule that applies a constant annual exploitation
#' rate to the current biomass estimate.
#'
#' @param u_target Numeric > 0. Target exploitation rate (fraction of biomass
#'   removed per year).
#'
#' @return A function with signature \code{function(biomass, reference_points)}
#'   that returns TAC = \code{u_target * biomass}.
#'
#' @examples
#' hcr <- hcr_constant_u(0.1)
#' hcr(1000, list(BMSY = 500))
#' # Returns 100
#'
#' @export
hcr_constant_u <- function(u_target) {
  assert_number(u_target, lower = .Machine$double.eps, .var.name = "u_target")

  function(biomass, reference_points) {
    tac <- u_target * biomass
    max(0, tac)
  }
}


#' Legacy Alias for Constant Exploitation-Rate HCR
#'
#' `hcr_constant_f()` is retained for backwards compatibility. Its argument is
#' an annual exploitation rate, not instantaneous fishing mortality. Use
#' \code{hcr_constant_u()} in new analyses.
#'
#' @param f_target Numeric > 0. Target annual exploitation rate; the name is
#'   retained for backwards compatibility only and is not an instantaneous
#'   fishing mortality.
#'
#' @export
hcr_constant_f <- function(f_target) {
  hcr_constant_u(f_target)
}


#' Constant Catch HCR
#'
#' Returns a harvest control rule that recommends a fixed catch each year,
#' clamped to the available biomass.
#'
#' @param catch_target Numeric > 0. Fixed annual catch target.
#'
#' @return A function with signature \code{function(biomass, reference_points)}
#'   that returns TAC = \code{min(catch_target, biomass)}.
#'
#' @examples
#' hcr <- hcr_constant_catch(200)
#' hcr(1000, list(BMSY = 500))
#' # Returns 200
#'
#' hcr(50, list(BMSY = 500))
#' # Returns 50 (clamped to biomass)
#'
#' @export
hcr_constant_catch <- function(catch_target) {
  assert_number(catch_target,
    lower = .Machine$double.eps,
    .var.name = "catch_target"
  )

  function(biomass, reference_points) {
    min(catch_target, max(0, biomass))
  }
}


#' Hockey-Stick Exploitation-Rate HCR
#'
#' Returns a harvest control rule that applies a target harvest rate when
#' biomass is above \code{b_target}, linearly ramps down between
#' \code{b_limit} and \code{b_target}, and sets catch to zero below
#' \code{b_limit}. Biomass thresholds are expressed as fractions of the
#' observable estimation-model proxy \code{reference_points$B0}. The legacy
#' \code{K} element is accepted for backwards compatibility.
#'
#' @param u_target Numeric > 0. Target exploitation rate at and above \code{b_target}.
#' @param b_limit Numeric in (0, 1). Biomass depletion level (B/B0) below which
#'   catch is zero. Must be less than \code{b_target}.
#' @param b_target Numeric in (0, 1). Biomass depletion level (B/B0) above which
#'   the full \code{u_target} is applied.
#'
#' @return A function with signature \code{function(biomass, reference_points)}
#'   that returns a TAC value. \code{reference_points} must contain \code{B0},
#'   the EM-derived observable biomass reference.
#'
#' @details
#' The TAC is computed as:
#' \itemize{
#'   \item If \code{B/B0 >= b_target}: TAC = \code{u_target * B}
#'   \item If \code{b_limit < B/B0 < b_target}: TAC = \code{u_target * B *
#'     (B/B0 - b_limit) / (b_target - b_limit)}
#'   \item If \code{B/B0 <= b_limit}: TAC = 0
#' }
#'
#' @examples
#' hcr <- hcr_hockey_stick(f_target = 0.1, b_limit = 0.2, b_target = 0.4)
#'
#' # Above target: full F
#' hcr(3000, list(K = 5000))
#' # B/K = 0.6 > 0.4, TAC = 0.1 * 3000 = 300
#'
#' # In ramp zone
#' hcr(1500, list(K = 5000))
#' # B/K = 0.3, TAC = 0.1 * 1500 * (0.3 - 0.2)/(0.4 - 0.2) = 75
#'
#' # Below limit: zero catch
#' hcr(500, list(K = 5000))
#' # B/K = 0.1 < 0.2, TAC = 0
#'
#' @export
hcr_hockey_stick_u <- function(u_target, b_limit, b_target) {
  assert_number(u_target, lower = .Machine$double.eps, .var.name = "u_target")
  assert_number(b_limit,
    lower = .Machine$double.eps, upper = 1,
    .var.name = "b_limit"
  )
  assert_number(b_target,
    lower = .Machine$double.eps, upper = 1,
    .var.name = "b_target"
  )

  if (b_limit >= b_target) {
    stop("b_limit (", b_limit, ") must be less than b_target (", b_target, ")",
      call. = FALSE
    )
  }

  function(biomass, reference_points) {
    B0 <- reference_points$B0 %||% reference_points$K
    if (is.null(B0) || !is.finite(B0) || B0 <= 0) {
      stop("reference_points must contain positive finite 'B0' (legacy alias: 'K')", call. = FALSE)
    }

    depletion <- biomass / B0

    if (depletion >= b_target) {
      tac <- u_target * biomass
    } else if (depletion > b_limit) {
      multiplier <- (depletion - b_limit) / (b_target - b_limit)
      tac <- u_target * biomass * multiplier
    } else {
      tac <- 0
    }

    max(0, tac)
  }
}


#' Legacy Alias for Hockey-Stick Exploitation-Rate HCR
#'
#' `hcr_hockey_stick()` is retained for backwards compatibility. Its first
#' argument is an annual exploitation rate. Use \code{hcr_hockey_stick_u()} in new
#' analyses.
#'
#' @param f_target Numeric > 0. Target annual exploitation rate at and above
#'   \code{b_target}; the name is retained for backwards compatibility only.
#' @param b_limit Numeric in (0, 1). Biomass depletion level (B/B0) below
#'   which catch is zero. Must be less than \code{b_target}.
#' @param b_target Numeric in (0, 1). Biomass depletion level (B/B0) above
#'   which the full \code{f_target} exploitation rate is applied.
#'
#' @export
hcr_hockey_stick <- function(f_target, b_limit, b_target) {
  hcr_hockey_stick_u(f_target, b_limit, b_target)
}


#' Convert an annual exploitation rate to instantaneous fishing mortality
#'
#' This conversion is required when a selected annual exploitation-rate HCR is
#' implemented in Casal2 using a Baranov fishing-mortality parametrisation.
#'
#' @param u Annual exploitation rate in [0, 1).
#' @return Instantaneous fishing mortality, `-log(1 - u)`.
#' @export
u_to_f <- function(u) {
  assert_numeric(u, lower = 0, upper = 1, .var.name = "u")
  if (any(u >= 1)) stop("u must be less than one for conversion to F", call. = FALSE)
  -log1p(-u)
}


#' Convert instantaneous fishing mortality to an annual exploitation rate
#'
#' @param f Instantaneous fishing mortality, non-negative.
#' @return Annual exploitation rate, `1 - exp(-f)`.
#' @export
f_to_u <- function(f) {
  assert_numeric(f, lower = 0, .var.name = "f")
  -expm1(-f)
}


#' CCAMLR-Style Precautionary HCR
#'
#' Returns a harvest control rule inspired by CCAMLR's precautionary approach.
#' The rule targets a long-term depletion level of \code{target_depletion}
#' (default 0.5, i.e. 50\% of \eqn{K}) with a probability cap on
#' depleting below \code{limit_depletion} (default 0.2).
#'
#' The TAC is set to maintain biomass near \code{target_depletion * K},
#' with a linear ramp-down when biomass falls below the target and a
#' complete closure below the limit.
#'
#' @param gamma Numeric > 0. Fraction of the estimated surplus production to
#'   harvest (default 0.5, i.e. harvest half of MSY-equivalent production).
#' @param target_depletion Numeric in (0, 1). Target B/K depletion level
#'   (default 0.5).
#' @param limit_depletion Numeric in (0, 1). Limit B/K depletion level below
#'   which catch is zero (default 0.2). Must be less than
#'   \code{target_depletion}.
#'
#' @return A function with signature \code{function(biomass, reference_points)}
#'   that returns a TAC value. \code{reference_points} must contain \code{K}
#'   and \code{MSY}.
#'
#' @examples
#' hcr <- hcr_ccamlr_krill()
#' hcr(3000, list(K = 5000, MSY = 400))
#'
#' @export
hcr_ccamlr_krill <- function(gamma = 0.5,
                             target_depletion = 0.5,
                             limit_depletion = 0.2) {
  assert_number(gamma, lower = .Machine$double.eps, .var.name = "gamma")
  assert_number(target_depletion,
    lower = .Machine$double.eps, upper = 1,
    .var.name = "target_depletion"
  )
  assert_number(limit_depletion,
    lower = .Machine$double.eps, upper = 1,
    .var.name = "limit_depletion"
  )

  if (limit_depletion >= target_depletion) {
    stop("limit_depletion (", limit_depletion,
      ") must be less than target_depletion (", target_depletion, ")",
      call. = FALSE
    )
  }

  function(biomass, reference_points) {
    K <- reference_points$K
    MSY <- reference_points$MSY
    if (is.null(K) || !is.finite(K) || K <= 0) {
      stop("reference_points must contain positive finite 'K'", call. = FALSE)
    }
    if (is.null(MSY) || !is.finite(MSY) || MSY <= 0) {
      stop("reference_points must contain positive finite 'MSY'", call. = FALSE)
    }

    depletion <- biomass / K
    base_tac <- gamma * MSY

    if (depletion >= target_depletion) {
      tac <- base_tac
    } else if (depletion > limit_depletion) {
      multiplier <- (depletion - limit_depletion) /
        (target_depletion - limit_depletion)
      tac <- base_tac * multiplier
    } else {
      tac <- 0
    }

    max(0, tac)
  }
}
