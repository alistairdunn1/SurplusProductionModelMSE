# MSE: Management Strategy Evaluation for Antarctic Toothfish

<!-- badges: start -->
[![R-CMD-check](https://github.com/alistairdunn1/MSE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alistairdunn1/MSE/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/alistairdunn1/MSE/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/alistairdunn1/MSE/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

## Overview

**MSE** provides a closed-loop Management Strategy Evaluation framework for Antarctic toothfish (*Dissostichus mawsoni*). It uses the [SurplusProductionModel](https://github.com/alistairdunn1/SurplusProductionModel) package as the assessment engine and supports:
 
- **Operating model** forward projection with Pella-Tomlinson dynamics and spatial movement
- **Observation error** simulation calibrated from fitted model residuals
- **Implementation error** (catch vs TAC) with lognormal errors and temporal autocorrelation
- **Harvest control rules**: constant F, constant catch, hockey-stick, CCAMLR-style
- **Performance metrics**: AAV, depletion risk, yield, biomass ratio
- **Spatial robustness testing**: OM/EM structural mismatch scenarios
- **Parallel processing** via the `future` package

## Installation

```r
# Install from GitHub (requires SurplusProductionModel)
# install.packages("remotes")
remotes::install_github("alistairdunn1/SurplusProductionModel")
remotes::install_github("alistairdunn1/MSE")
```

## Quick Example

```r
library(MSE)

# Configure operating model (single-area Schaefer)
om <- om_config(
  n_areas = 1L,
  true_params = list(
    r = 0.3, K = 5000, m = 2,
    sigma_obs = 0.2, q = 0.001, B0 = 5000
  )
)

# Define two harvest control rules
# 1. Constant catch at 80% of MSY (MSY = rK/8 = 187.5)
# 2. Hockey-stick that reduces F when biomass is low
scenarios <- list(
  create_scenario("Constant catch", hcr_constant_catch(150)),
  create_scenario("Hockey stick",   hcr_hockey_stick(0.075, 0.20, 0.40))
)

# Run MSE (skip estimation model for speed)
result <- mse_simulation(
  operating_model  = om,
  scenarios        = scenarios,
  n_sims           = 50L,
  n_proj_years     = 20L,
  min_assess_years = 99L,
  seed             = 42L
)
print(result)

# Evaluate performance for each scenario
perf <- lapply(result$results, function(sc) {
  calculate_performance_metrics(sc$trajectories, om)
})

# Attach performance and visualise trade-offs
for (nm in names(perf)) {
  result$results[[nm]]$performance <- perf[[nm]]
}
compare_scenarios(result, "mean_catch", "Pr(B<20%B0)_ever")
```

## Package Structure

| Module | Description |
|--------|-------------|
| Operating Model | Forward projection with Pella-Tomlinson dynamics and spatial movement |
| Observation Error | Future CPUE simulation calibrated from historical residuals |
| Implementation Error | Stochastic catch realisation from TAC advice |
| Harvest Control Rules | Predefined and custom HCR functions |
| Performance Metrics | AAV, depletion risk, yield, biomass-based metrics |
| Spatial Robustness | OM/EM mismatch testing scenarios |

## Dependencies

- [SurplusProductionModel](https://github.com/alistairdunn1/SurplusProductionModel) (>= 0.1.0) — assessment engine
- [checkmate](https://CRAN.R-project.org/package=checkmate) — input validation
- [ggplot2](https://CRAN.R-project.org/package=ggplot2) — visualisation
- [future](https://CRAN.R-project.org/package=future) / [future.apply](https://CRAN.R-project.org/package=future.apply) — parallel processing (optional)

## Related Packages

This package is part of the ATO rTMB project:

1. **SurplusProductionModel** — Pella-Tomlinson surplus production model (operating model / estimation model)
2. **MSE** — Management Strategy Evaluation framework (this package)
3. **IntegratedAgelengthModel** — Integrated age-length model (future)

## License

MIT
