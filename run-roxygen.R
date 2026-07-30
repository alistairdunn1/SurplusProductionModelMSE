# =====================================
# Generate Documentation with roxygen2
# =====================================

# Load required packages with better error handling
load_package <- function(pkg_name) {
  if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg_name, "package...\n")
    install.packages(pkg_name, repos = "https://cran.r-project.org/")
    if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
      cat("ERROR: Failed to install", pkg_name, "package.\n")
      quit(status = 1)
    }
  }
}

# Load required packages
load_package("devtools")
load_package("roxygen2")

# Generate documentation for the SurplusProductionModelMSE package
cat("Generating documentation for SurplusProductionModelMSE package...\n")

tryCatch(
  {
    # Document the package
    devtools::document(".")
    cat("Documentation generated successfully!\n")

    # Verify documentation was created
    man_files <- list.files("man", pattern = "\\.Rd$")
    cat("Created", length(man_files), "documentation files.\n")
  },
  error = function(e) {
    cat("ERROR generating documentation:", e$message, "\n")
    quit(status = 1)
  }
)

cat("Documentation generation completed successfully.\n")
