\name{test_ks_dexp}
\alias{test_ks_dexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov-Smirnov Tests (Exponential Distribution)
}
\description{
Performs one sample Kolmogorov-Smirnov tests.
}
\usage{
test_ks_dexp(X, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of data values.
}
  \item{lambda}{
vector of rates.
}
}
\details{
see detail \code{ks.test}.
}
\value{
A list with class "htest" containing the following components: 
\item{statistic}{the value of the test statistic.}
\item{p.value}{	the p-value of the test.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{data.name}{a character string giving the name(s) of the data.}
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\examples{
## Example 1
X <- rexp(1000,c(1,2,3))
test_ks_dexp(X, lambda=1)
test_ks_dexp(X, lambda=2)
test_ks_dexp(X, lambda=3)
## Example 2
X <- rexp(1000,3)
test_ks_dexp(X, lambda=3)
test_ks_dweibull(X, shape=1, scale=(1/3))
test_ks_dgamma(X, shape=1, rate=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
