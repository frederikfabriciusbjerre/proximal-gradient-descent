% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulate_independent}
\alias{simulate_independent}
\title{Simulate independent design}
\usage{
simulate_independent(
  n,
  p,
  k,
  error_dist = c("normal", "t", "uniform"),
  beta = NULL,
  family = c("gaussian", "binomial"),
  seed = 96
)
}
\arguments{
\item{n}{Number of observations}

\item{p}{Number of features}

\item{k}{Number of non-zero coefficients}

\item{error_dist}{Error distribution: "normal", "t", "uniform"}

\item{beta}{Optional true coefficient vector}

\item{seed}{Random seed}
}
\value{
List with X, y, beta
}
\description{
Simulate independent design
}
