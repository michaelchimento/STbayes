% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_user_STb2.R
\name{import_user_STb2}
\alias{import_user_STb2}
\title{Import user data 2}
\usage{
import_user_STb2(
  event_data,
  networks,
  network_type = c("undirected", "directed"),
  ILV_c = NULL,
  ILV_tv = NULL,
  ILVi = NULL,
  ILVs = NULL,
  ILVm = NULL,
  t_weights = NULL,
  high_res = FALSE
)
}
\arguments{
\item{event_data}{dataframe with columns id, trial, time, t_end}

\item{networks}{Either a dataframe, a bisonr fit, or a list of bisonr fits.
If dataframe: with columns trial, from, to, and one or more columns of edge weights named descriptively.}

\item{network_type}{"undirected" or "directed".}

\item{ILV_c}{optional dataframe with columns id, and any constant individual-level variables that might be of interest}

\item{ILV_tv}{optional dataframe with columns trial, id, time and any time-varying variables. Variable values should summarize the variable for each inter-acquisition period.}

\item{ILVi}{Optional character vector of column names from ILV metadata to be considered when estimating intrinsic rate. If not specified, all ILV are applied to both.}

\item{ILVs}{Optional character vector of column names from ILV metadata to be considered when estimating social transmission rate. If not specified, all ILV are applied to both.}

\item{ILVm}{Optional character vector of column names from ILV metadata to be considered in a multiplicative model.}

\item{t_weights}{Optional dataframe with columns trial, id, time and t_weight. Transmission rates represent rates of production/relevant cues per inter-event period.}

\item{high_res}{Boolean indicating whether or not user is providing networks and transmission weights per period duration=1}
}
\value{
A list object containing properly formatted data to run social transmission models.
}
\description{
First step of analysis pipeline: create STbayes data object from user
supplied data to be used for generating and fitting models. This function is
basically used when \emph{fitting models of complex transmission to high-resolution data.}
Rather than pre-process high res data, this will just create a massive
data-list to be sent to Stan. Models created with this function will take much longer to run.
}
