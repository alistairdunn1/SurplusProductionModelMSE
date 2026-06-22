# SurplusProductionModelMSE 0.1.1

## Robustness fixes

* `mse_simulation()` now caps realised catch at the available biomass. Removals
  in a year cannot exceed a maximum exploitation fraction of the current
  biomass (default 0.95, overridable via `om_config$max_harvest_rate`). This
  keeps recorded catch physical and prevents an implausible catch limit, for
  example from a divergent estimation-model fit, from driving biomass negative.

* The estimation-model fit is now subject to a plausibility guard. An estimated
  biomass or carrying capacity beyond a generous multiple of the operating-model
  carrying capacity (default 10x, overridable via `em_config$max_biomass_factor`)
  is rejected and treated as a convergence failure, so the previous assessment
  is carried forward rather than driving the harvest control rule to an extreme
  catch limit.

* Average annual variation (AAV) in catch is now computed with the sum-based
  form, `mean_s( sum_t |C_{t+1} - C_t| / sum_t C_t )`, which is robust to
  individual near-zero catch years. The previous per-year ratio form divided
  each change by that year's catch and could return extreme values when catch
  approached zero.
