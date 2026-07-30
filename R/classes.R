# Validation functions for MSE S3 classes
#
# Each validate_* function takes an object and checks structural integrity.
# Returns TRUE invisibly on success; throws an error on failure.

#' Validate an om_config object
#'
#' @param x An object of class \code{om_config}.
#' @return \code{TRUE} invisibly if valid; otherwise an error is raised.
#' @export
validate_om_config <- function(x) {
  if (!inherits(x, "om_config")) {
    stop("Object is not of class 'om_config'", call. = FALSE)
  }

  required <- c(
    "n_areas", "movement_rate", "movement_cost_matrix",
    "attractiveness", "decay", "true_params"
  )
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0) {
    stop("om_config missing fields: ", paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  assert_count(x$n_areas, positive = TRUE, .var.name = "n_areas")
  assert_number(x$movement_rate,
    lower = 0, upper = 1,
    .var.name = "movement_rate"
  )
  assert_number(x$decay, lower = 0, .var.name = "decay")

  na <- x$n_areas

  # movement_cost_matrix

  if (!is.null(x$movement_cost_matrix)) {
    if (!is.matrix(x$movement_cost_matrix) ||
      nrow(x$movement_cost_matrix) != na ||
      ncol(x$movement_cost_matrix) != na) {
      stop("movement_cost_matrix must be a ", na, " x ", na, " matrix",
        call. = FALSE
      )
    }
    if (any(diag(x$movement_cost_matrix) != 0)) {
      stop("movement_cost_matrix diagonal must be zero", call. = FALSE)
    }
    if (any(x$movement_cost_matrix < 0)) {
      stop("movement_cost_matrix must contain non-negative values", call. = FALSE)
    }
  }

  # attractiveness
  if (!is.null(x$attractiveness)) {
    assert_numeric(x$attractiveness,
      len = na, lower = 0, any.missing = FALSE,
      .var.name = "attractiveness"
    )
  }

  # true_params
  if (!is.null(x$true_params)) {
    assert_list(x$true_params, .var.name = "true_params")
    required_params <- c("r", "K", "m", "sigma_obs", "q", "B_initial")
    missing_params <- setdiff(required_params, names(x$true_params))
    if (length(missing_params) > 0) {
      stop("true_params missing: ", paste(missing_params, collapse = ", "),
        call. = FALSE
      )
    }
    tp <- x$true_params
    assert_number(tp$r,
      lower = .Machine$double.eps, upper = 2,
      .var.name = "true_params$r"
    )
    assert_number(tp$K,
      lower = .Machine$double.eps,
      .var.name = "true_params$K"
    )
    if (!is.null(tp$K_area)) {
      assert_numeric(tp$K_area,
        lower = .Machine$double.eps, any.missing = FALSE,
        min.len = 1, max.len = na, .var.name = "true_params$K_area"
      )
    }
    assert_number(tp$m,
      lower = .Machine$double.eps,
      .var.name = "true_params$m"
    )
    assert_number(tp$sigma_obs,
      lower = .Machine$double.eps,
      .var.name = "true_params$sigma_obs"
    )
    # q and B_initial can be scalar or per-area vector
    assert_numeric(tp$q,
      lower = .Machine$double.eps, any.missing = FALSE,
      min.len = 1, max.len = na, .var.name = "true_params$q"
    )
    assert_numeric(tp$B_initial,
      lower = .Machine$double.eps, any.missing = FALSE,
      min.len = 1, max.len = na, .var.name = "true_params$B_initial"
    )
  }

  if (!is.null(x$max_harvest_rate)) {
    assert_number(x$max_harvest_rate,
      lower = .Machine$double.eps, upper = 1,
      .var.name = "max_harvest_rate"
    )
  }

  invisible(TRUE)
}


#' Validate an em_config object
#'
#' @param x An object of class \code{em_config}.
#' @return \code{TRUE} invisibly if valid; otherwise an error is raised.
#' @export
validate_em_config <- function(x) {
  if (!inherits(x, "em_config")) {
    stop("Object is not of class 'em_config'", call. = FALSE)
  }

  required <- c(
    "n_areas", "estimate_movement", "aggregate_areas",
    "process_noise", "process_error_structure", "fixed_params",
    "initial_depletion"
  )
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0) {
    stop("em_config missing fields: ", paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  assert_count(x$n_areas, positive = TRUE, .var.name = "n_areas")
  assert_flag(x$estimate_movement, .var.name = "estimate_movement")
  assert_flag(x$aggregate_areas, .var.name = "aggregate_areas")
  assert_flag(x$process_noise, .var.name = "process_noise")
  assert_character(
    x$process_error_structure,
    len = 1,
    any.missing = FALSE,
    .var.name = "process_error_structure"
  )
  if (!x$process_error_structure %in% c("iid", "ar1")) {
    stop("process_error_structure must be one of: iid, ar1", call. = FALSE)
  }

  if (!is.null(x$fixed_params)) {
    assert_list(x$fixed_params, names = "named", .var.name = "fixed_params")
  }
  if (!is.null(x$initial_depletion)) {
    assert_number(x$initial_depletion,
      lower = .Machine$double.eps, upper = 1,
      .var.name = "initial_depletion"
    )
  }

  invisible(TRUE)
}


#' Validate an impl_error object
#'
#' @param x An object of class \code{impl_error}.
#' @return \code{TRUE} invisibly if valid; otherwise an error is raised.
#' @export
validate_impl_error <- function(x) {
  if (!inherits(x, "impl_error")) {
    stop("Object is not of class 'impl_error'", call. = FALSE)
  }

  required <- c("cv", "bias", "autocorr", "max_overage")
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0) {
    stop("impl_error missing fields: ", paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  assert_number(x$cv, lower = .Machine$double.eps, .var.name = "cv")
  assert_number(x$bias, lower = 0, .var.name = "bias")
  assert_number(x$autocorr, lower = 0, upper = 1, .var.name = "autocorr")
  assert_number(x$max_overage, lower = 1, .var.name = "max_overage")

  invisible(TRUE)
}


#' Validate an mse_scenario object
#'
#' @param x An object of class \code{mse_scenario}.
#' @return \code{TRUE} invisibly if valid; otherwise an error is raised.
#' @export
validate_mse_scenario <- function(x) {
  if (!inherits(x, "mse_scenario")) {
    stop("Object is not of class 'mse_scenario'", call. = FALSE)
  }

  required <- c(
    "name", "harvest_control_rule", "implementation_error",
    "assessment_frequency", "catch_allocation"
  )
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0) {
    stop("mse_scenario missing fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  assert_character(x$name, len = 1, any.missing = FALSE, .var.name = "name")
  assert_function(x$harvest_control_rule, .var.name = "harvest_control_rule")

  # HCR must accept biomass and reference_points arguments
  hcr_args <- names(formals(x$harvest_control_rule))
  if (length(hcr_args) < 2) {
    stop("harvest_control_rule must accept at least 2 arguments ",
      "(biomass, reference_points)",
      call. = FALSE
    )
  }

  if (!is.null(x$implementation_error)) {
    if (!inherits(x$implementation_error, "impl_error")) {
      stop("implementation_error must be an 'impl_error' object or NULL",
        call. = FALSE
      )
    }
  }

  assert_count(x$assessment_frequency,
    positive = TRUE,
    .var.name = "assessment_frequency"
  )

  if (!is.null(x$catch_allocation)) {
    assert_numeric(
      x$catch_allocation,
      lower = 0,
      any.missing = FALSE,
      min.len = 1,
      .var.name = "catch_allocation"
    )
    if (sum(x$catch_allocation) <= 0) {
      stop("catch_allocation must contain at least one positive value",
        call. = FALSE
      )
    }

    if (length(x$catch_allocation) > 1) {
      nm <- names(x$catch_allocation)
      if (is.null(nm) || any(!nzchar(nm))) {
        stop(
          "For multi-area use, catch_allocation must be a named vector",
          call. = FALSE
        )
      }
      if (anyDuplicated(nm) > 0) {
        stop("catch_allocation names must be unique", call. = FALSE)
      }
    }
  }

  invisible(TRUE)
}
