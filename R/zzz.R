.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "SurplusProductionModelMSE v",
    utils::packageVersion("SurplusProductionModelMSE"),
    " loaded."
  )
  packageStartupMessage(
    "Management Strategy Evaluation framework."
  )
}

.onLoad <- function(libname, pkgname) {
  # Check for SurplusProductionModel availability
  if (!requireNamespace("SurplusProductionModel", quietly = TRUE)) {
    warning(
      "SurplusProductionModel package not available. ",
      "Install with: remotes::install_github('alistairdunn1/SurplusProductionModel')"
    )
  }
}

# Clean up on unload
.onUnload <- function(libpath) {
  # Clean up any global variables or connections if needed
}
