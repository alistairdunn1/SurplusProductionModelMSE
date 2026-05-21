# MSE 0.1.0

## Initial Release

### New Features
* Package scaffolding and CI/CD infrastructure
* Closed-loop MSE simulation framework for Antarctic toothfish
* Operating model forward projection with Pella-Tomlinson dynamics
* Observation error simulation calibrated from fitted model residuals
* Implementation error with lognormal distribution and temporal autocorrelation
* Predefined harvest control rules: constant F, constant catch, hockey-stick, CCAMLR-style
* Performance metrics: AAV, depletion risk, yield, biomass ratio
* Spatial robustness testing with OM/EM structural mismatch scenarios
* Parallel processing support via `future` package
* Trade-off visualisation (Pareto frontier) for scenario comparison

### Package Structure
* Configuration classes: `om_config`, `em_config`, `impl_error`, `mse_scenario`
* Core functions: `mse_simulation()`, `performance_metrics()`, `compare_scenarios()`
* Harvest control rules: `hcr_constant_f()`, `hcr_constant_catch()`, `hcr_hockey_stick()`

### Documentation
* Complete roxygen2 documentation for all exported functions
* Getting-started vignette with end-to-end MSE example
* Comprehensive README with examples

### Testing
* Unit tests for all major functions
* Integration tests for full simulation workflows
* Self-test validation (EM = OM structure)
