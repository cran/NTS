% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TAR.r
\name{PRnd}
\alias{PRnd}
\title{ND Test}
\usage{
PRnd(x, m = 10, p = 0, q = 0)
}
\arguments{
\item{x}{time series.}

\item{m}{the maximum number of lag of correlation to test.}

\item{p}{AR order.}

\item{q}{MA order.}
}
\value{
\code{PRnd} function outputs the ND test statistic and its p-value.
}
\description{
Compute the ND test statistic of Pena and Rodriguez (2006, JSPI).
}
\examples{
y=arima.sim(n=500,list(ar=c(0.8,-0.6,0.7)))
PRnd(y,10,3,0)
}
\references{
Pena, D., and Rodriguez, J. (2006) A powerful Portmanteau test of lack of fit for time series. series. \emph{Journal of American Statistical Association}, 97, 601-610.
}
