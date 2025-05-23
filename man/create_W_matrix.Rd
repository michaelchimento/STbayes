% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_W_matrix.R
\name{create_W_matrix}
\alias{create_W_matrix}
\title{Helper function to transform transmission weights into same size as C}
\usage{
create_W_matrix(t_weights, max_time)
}
\arguments{
\item{t_weights}{dataframe with t_weight. all ids across all trials must be present.}
}
\value{
W\link{k,t,n}
}
\description{
Helper function to transform transmission weights into same size as C
}
\examples{
t_weights <- data.frame(
  trial = c(rep(1, each = 9), rep(2, each = 9)),
  id = c(rep(LETTERS[1:3], each = 3), rep(LETTERS[4:6], each = 3)),
  time = c(rep(1:3, times = 3), rep(1:3, times = 3)),
  t_weight = exp(rnorm(18))
)
create_W_matrix(t_weights)
}
