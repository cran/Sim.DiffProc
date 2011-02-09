\name{test_ks_dlognorm}
\alias{test_ks_dlognorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov-Smirnov Tests (Log Normal Distribution)
}
\description{
Performs one sample Kolmogorov-Smirnov tests.
}
\usage{
test_ks_dlognorm(X, meanlog, sdlog)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of data values.
}
  \item{meanlog}{
mean of the distribution.
}
  \item{sdlog}{
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
X <- rlnorm(1000,1,1)
test_ks_dlognorm(X, meanlog=1, sdlog=1)
test_ks_dnorm(log(X), mean=1, sd=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
