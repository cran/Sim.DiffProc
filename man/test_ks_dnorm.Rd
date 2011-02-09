\name{test_ks_dnorm}
\alias{test_ks_dnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov-Smirnov Tests (Normal Distribution)
}
\description{
Performs one sample Kolmogorov-Smirnov tests.
}
\usage{
test_ks_dnorm(X, mean, sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of data values.
}
  \item{mean}{
mean of the distribution.
}
  \item{sd}{
standard deviation of the distribution.
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
X <- rnorm(1000,1,1)
test_ks_dnorm(X, mean=1, sd=1)
test_ks_dlognorm(exp(X), meanlog=1, sdlog=1)
## Example 2
X = c(runif(100),rt(200,20),rnorm(200))
X = sample(X)
test_ks_dnorm(X, mean=mean(X), sd=sd(X))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
