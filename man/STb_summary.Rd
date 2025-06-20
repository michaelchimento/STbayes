% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STb_summary.R
\name{STb_summary}
\alias{STb_summary}
\title{STb_summary()}
\usage{
STb_summary(
  fit,
  depth = 1,
  prob = 0.95,
  ignore_params = c("lp__", "idx", "log_lik", "log_lik_matrix", "count_ST", "psocn_sum",
    "acquisition_time", "z_ID", "Rho_ID", "v_ID", ".chain", ".iteration", ".draw",
    "s_prime"),
  digits = 3,
  CI_method = c("HPDI", "PI")
)
}
\arguments{
\item{fit}{CmdStanMCMC model fit}

\item{depth}{integer depth of multidimensional parameters to extract}

\item{prob}{double limits for HPD of estimates (default = 0.95)}

\item{ignore_params}{character vector of parameters to ignore}

\item{digits}{integer of digits to round to}

\item{CI_method}{"HPDI" for highest density interval or "PI" for quantiles (equal tails). Defaults to HPDI.}
}
\value{
Summary table
}
\description{
Create a useful summary table for STbayes fit. Contains estimates on log and
linear scale for lambda_0 (intrinsic rate) and s' (social transmission rate).
Includes estimated s (relative strength of social transmission).
}
\examples{
data_list <- import_user_STb(STbayes::event_data, STbayes::edge_list)
model_obj <- generate_STb_model(data_list)
fit <- fit_STb(data_list,
    model_obj,
    parallel_chains = 4,
    chains = 4,
    cores = 4,
    iter = 4000,
    refresh = 2000
)
STb_summary(fit, prob = 0.95, digits = 3, CI_method = "HPDI")
}
