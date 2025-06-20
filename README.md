<br><br>
![banner](man/figures/stbayes_banner2.png)
<!-- badges: start -->
<!-- badges: end -->

**STbayes** (ST standing for social transmission) is a package for building and running Bayesian inferential models of social transmission across static or dynamic networks. Users may supply their own data in formats given below, or import nbdaData objects directly from the [NBDA package](https://github.com/whoppitt/NBDA).
   
There is now an associated [preprint](https://doi.org/10.1101/2025.06.07.658152). The package is still under development, so if you run into a bug, please open an issue or contact me directly!

## Installation<a name="Installation"></a>

The functions of this package depend on ```cmdstanr```, ```posterior```, ```bayestestR```, ```data.table``` and ```loo```. You can install ```cmdstanr``` by following these [instructions](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

Vignettes use packages NBDA (install with ```devtools::install_github("whoppitt/NBDA"```), igraph, dplyr, ggplot2, and ggpubr.


To install the development version of STbayes:

``` r
# install devtools if not already
if (!require("devtools")) install.packages("devtools")
devtools::install_github("michaelchimento/STbayes")
```

## Getting started

- [Importing data](https://michaelchimento.github.io/STbayes/articles/getting_started.html#step-1-importing-data)
- [Generate a model](https://michaelchimento.github.io/STbayes/articles/getting_started.html#step-2-generate-a-model)
- [Fit and save the model](https://michaelchimento.github.io/STbayes/articles/getting_started.html#step-3-fit-and-save-the-model)
- [Viewing model output](https://michaelchimento.github.io/STbayes/articles/getting_started.html#step-4-viewing-model-output)
- [Posterior predictive checks](https://michaelchimento.github.io/STbayes/articles/getting_started.html#step-5-posterior-predictive-checks)
- [Model comparison](https://michaelchimento.github.io/STbayes/articles/getting_started.html#step-6-model-comparison)

## Advanced recipes

- [Multi-network models](https://michaelchimento.github.io/STbayes/articles/advanced_recipes.html#multi-network-models)
- [Individual-level variables](https://michaelchimento.github.io/STbayes/articles/advanced_recipes.html#individual-level-variables-ilvs)
- [Transmission weights](https://michaelchimento.github.io/STbayes/articles/advanced_recipes.html#transmission-weights)
- [Varying effects](https://michaelchimento.github.io/STbayes/articles/getting_started.html#varying-effects-by-individual)
- [Other model types: OADA and dTADA](https://michaelchimento.github.io/STbayes/articles/getting_started.html#other-model-types-oada-and-dtada)
- [Edge weight uncertainty](https://michaelchimento.github.io/STbayes/articles/getting_started.html#edge-weight-uncertainty)
- [Complex Transmission](https://michaelchimento.github.io/STbayes/articles/getting_started.html#complex-transmission)
- [Setting priors](https://michaelchimento.github.io/STbayes/articles/advanced_recipes.html#setting-priors)
- [High resolution data mode](https://michaelchimento.github.io/STbayes/articles/advanced_recipes.html#high-resolution-data-mode)
- [Import NBDA Objects](https://michaelchimento.github.io/STbayes/articles/advanced_recipes.html#import-nbda-data-objects)

## Features:

 - TADA (acquisition time known) and OADA (only acquisition order known) model types.
 - static and dynamic networks.
 - multi-network comparison (with static or dynamic networks).
 - multiple trials with the same set, subsets, or different sets of individuals.
 - constant and time-varying ILVs for additive and multiplicative transmission models.
 - varying effects by individual for strength of social transmission, baseline hazard rates, and other user defined ILVs.
 - easy workflow for ELPD (loo-psis, waic) model comparison.
 - propagation of uncertainty from network measures to transmission model
 - modeling of complex transmission and contagion