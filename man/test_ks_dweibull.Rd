\name{test_ks_dweibull}
\alias{test_ks_dweibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov-Smirnov Tests (Weibull Distribution)
}
\description{
Performs one sample Kolmogorov-Smirnov tests.
}
\usage{
test_ks_dweibull(X, shape, scale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of data values.
}
  \item{shape}{
shape and scale parameters, the latter defaulting to 1.
}
  \item{scale}{
shape and scale parameters, the latter defaulting to 1.
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
Boukhetala Kamal, Guidoum Arsalane.
}
\examples{
 X <- rweibull(1000,1,4)
 test_ks_dweibull(X, shape=1, scale=4)
 test_ks_dexp(X, lambda=0.25)
 test_ks_dgamma(X, shape=1, rate=0.25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
