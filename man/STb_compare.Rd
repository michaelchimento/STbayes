% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STb_compare.R
\name{STb_compare}
\alias{STb_compare}
\title{STb_compare() automates the workflow of loo-psis or waic elpd comparisons. Relies on loo to do this.}
\usage{
STb_compare(..., model_names = NULL, method = "loo-psis")
}
\arguments{
\item{...}{CmdStanMCMC model fits for comparison}

\item{model_names}{an optional list of model names, otherwise taken from object names}

\item{method}{a string either "loo-psis" or "waic" to indicate the method used for elpd. Defaults to "loo-psis". Pareto diagnostics are not calculated for WAIC.}
}
\value{
list containing loo_objects, comparison, and pareto_diagnostics if using loo-psis.
}
\description{
STb_compare() automates the workflow of loo-psis or waic elpd comparisons. Relies on loo to do this.
}
