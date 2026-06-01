# STbayes 1.1.0

## New features

- Package now supports PSIS-LFO-LOOCV as another tool for for model comparison (inspiration from [Bürkner, Gabry & Vehtari 2020](https://cran.r-project.org/web/packages/loo/vignettes/loo2-lfo.html)).

- Package now supports time-varying and trial-varying edge-weights when networks are given as posterior draws. Before, varying networks were only supported when edge-weights were point estimates. 

- Estimating acquisition times in generated quantities block has been optimized, now supports when using float event times and when supplying posterior draws of edge-weights.

## Changes

- If importing network edge-weights as draws from posterior distributions, arrays must have named dimensions for my sanity and yours. Please see [this vignette](../articles/advanced_recipes.html#edge-weight-uncertainty) for more details.
- The example code for plotting PPCs in the [Getting Started vignette](../articles/getting_started.html) is now more generalisable, addresses cases where trials may have different numbers of individuals, and clarifies how to handle demonstrators.

## Bug fixes

- Fixed bug where demonstrators were not excluded from prediction set when using ```generate_STb_model(est_AcqTime=T)```
- Fixed bug when estimating acquisition times using Weibull shaped hazards

## Documentation

- Added [vignette for Weibull-shaped hazards](../articles/weibull_shaped_hazards.html).
- Added [vignette for LFO-CV](../articles/LFO-CV_model_comparison.html)

# STbayes 1.0.0

- Initial release accompanying the published paper.
- Note: this release was tagged as V1.0 on GitHub, but the package version field was mistakenly left as `0.0.0.9000`.

# STbayes 0.0.0.9000

- development version
