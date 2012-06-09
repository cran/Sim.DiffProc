\name{test_ks_dgamma}
\alias{test_ks_dgamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov-Smirnov Tests (Gamma Distribution)
}
\description{
Performs one sample Kolmogorov-Smirnov tests.
}
\usage{
test_ks_dgamma(X, shape, rate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of data values.
}
  \item{shape}{
shape parameters. Must be positive, scale strictly.
}
  \item{rate}{
an alternative way to specify the scale.
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
 X <- rgamma(1000,1,6)
 test_ks_dgamma(X, shape=1, rate=6)
 test_ks_dexp(X, lambda=6)
 test_ks_dweibull(X, shape=1, scale=(1/6))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
