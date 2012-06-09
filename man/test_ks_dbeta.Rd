\name{test_ks_dbeta}
\alias{test_ks_dbeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov-Smirnov Tests (Beta Distribution)
}
\description{
Performs one sample Kolmogorov-Smirnov tests.
}
\usage{
test_ks_dbeta(X, shape1, shape2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of data values.
}
  \item{shape1}{
positive parameters of the Beta distribution.
}
  \item{shape2}{
positive parameters of the Beta distribution.
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
 X <- rbeta(1000,1,1)
 test_ks_dbeta(X, shape1=1, shape2=1)
 test_ks_dbeta(X, shape1=1, shape2=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
