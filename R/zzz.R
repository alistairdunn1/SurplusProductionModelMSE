.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "MSE v",
    utils::packageVersion("MSE"),
    " loaded."
  )
  packageStartupMessage(
    "Management Strategy Evaluation framework for Antarctic toothfish."
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
