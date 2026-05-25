#' SurplusProductionModelMSE: Management Strategy Evaluation
#'
#' A closed-loop Management Strategy Evaluation (MSE) framework. Uses the
#' \code{\link[SurplusProductionModel]{SurplusProductionModel}} package as the
#' assessment engine and provides operating model projection, observation and
#' implementation error simulation, harvest control rules, and performance
#' metric evaluation.
#'
#' @section Operating Model:
#'
#' \itemize{
#'   \item Pella-Tomlinson forward projection with spatial movement
#'   \item Configurable process noise
#'   \item Multi-area gravity movement kernel
#' }
#'
#' @section Observation Error:
#'
#' \itemize{
#'   \item Calibrated from fitted model residuals
#'   \item Lognormal CPUE simulation with bias correction
#'   \item Temporal autocorrelation support
#' }
#'
#' @section Implementation Error:
#'
#' \itemize{
#'   \item Lognormal catch vs TAC discrepancy
#'   \item Temporal autocorrelation
#'   \item Maximum overage constraints
#' }
#'
#' @section Harvest Control Rules:
#'
#' \itemize{
#'   \item \code{hcr_constant_f}: Constant harvest rate
#'   \item \code{hcr_constant_catch}: Constant catch
#'   \item \code{hcr_hockey_stick}: Ramp-down when biomass below target
#'   \item \code{hcr_ccamlr_krill}: CCAMLR-style precautionary rule
#' }
#'
#' @section Performance Metrics:
#'
#' \itemize{
#'   \item Average Annual Variation (AAV) in catch
#'   \item Depletion risk: P(B < threshold)
#'   \item Yield metrics: mean, median, CV
#'   \item Biomass metrics: B/BMSY ratio
#' }
#'
#' @section Spatial Robustness:
#'
#' \itemize{
#'   \item OM/EM structural mismatch testing
#'   \item Self-test validation
#'   \item Spatial aggregation scenarios
#' }
#'
#' @docType package
#' @name SurplusProductionModelMSE-package
#' @keywords internal
"_PACKAGE"
