\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{bconfint}
\alias{bconfint}
\alias{kurtosis}
\alias{moment}
\alias{skewness}
\alias{bconfint.default}
\alias{kurtosis.default}
\alias{moment.default}
\alias{skewness.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kurtosis, Skewness, Moment and Confidence Bands}
\description{Generic function for compute the kurtosis, skewness, moment and confidence bands of class \code{"sde"}.}
\usage{
\method{bconfint}{default}(x, level = 0.95, \dots)
\method{kurtosis}{default}(x, \dots)
\method{moment}{default}(x, order = 1,center = TRUE, \dots)
\method{skewness}{default}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object inheriting from class \code{"sde"}.}  
  \item{order}{order of moment.}
  \item{center}{central moment.}
  \item{level}{the confidence level required.} 
  \item{\dots}{further arguments for (non-default) methods.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\examples{
## Example 1:
## dX(t) = 2*(3-X(t)) *dt + dW(t)
set.seed(1234)

f <- expression( 2*(3-x) )
g <- expression( 1 )
res <- snssde1d(drift=f,diffusion=g,M=25)
x <- rsde1d(res,at=1)
kurtosis(res)
skewness(res)
bconfint(res,level = 0.95)
moment(res,order=4)
}
