% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMC.r
\name{SISstep.fading}
\alias{SISstep.fading}
\title{Sequential Importance Sampling Step for Fading Channels}
\usage{
SISstep.fading(mm, xx, logww, yyy, par, xdim2, ydim)
}
\arguments{
\item{mm}{the Monte Carlo sample size \code{m}.}

\item{xx}{the sample in the last iteration.}

\item{logww}{the log weight in the last iteration.}

\item{yyy}{the observations with \code{T} columns and \code{ydim} rows.}

\item{par}{a list of parameter values.  \code{HH} is the state coefficient model, \code{WW*t(WW)} is the state innovation covariance matrix,
\code{VV*t(VV)} is the covariance of the observation noise, \code{GG} is the observation model.}

\item{xdim2}{the dimension of the state variable \code{x_t}.}

\item{ydim}{the dimension of the observation \code{y_t}.}
}
\value{
The function returns a list with the following components:
\item{xx}{the new sample.}
\item{logww}{the log weights.}
}
\description{
This function implements one step of the sequential importance sampling method for fading channels.
}
\references{
Tsay, R. and Chen, R. (2018). Nonlinear Time Series Analysis. John Wiley & Sons, New Jersey.
}
