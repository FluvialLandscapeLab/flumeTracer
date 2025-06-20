% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_observations.R
\name{fitFlume}
\alias{fitFlume}
\title{Fit a flume model to observed data}
\usage{
fitFlume(par, obs_times, obs_C_c, m, lower, upper)
}
\arguments{
\item{par}{A named vector of model parameters fit to the observed data
(passed to \code{\link{optim}}; see that documentation for definition).  If
one of the parameters is named "delta_t" an time offset error from t = 0
will be included in the fit.}

\item{obs_times}{A vector of times where concentration values were observed.}

\item{obs_C_c}{A vetor of observed concentration values}

\item{m}{An environment describing a flume.  See \code{\link{createFlume}}.}

\item{lower, upper}{A vector of lower and upper boundaries for par values.
See "L-BFGS-B" method for \code{\link{optim}}, which is used for parameter
estimates.}
}
\description{
Fits any model parameters to observed data, potentially allowing for fitting
a time offset between modeled and observed data.
}
