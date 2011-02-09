\name{test_ks_dt}
\alias{test_ks_dt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov-Smirnov Tests (Student t Distribution)
}
\description{
Performs one sample Kolmogorov-Smirnov tests.
}
\usage{
test_ks_dt(X, df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of data values.
}
  \item{df}{
degrees of freedom (> 0, maybe non-integer). \code{df = Inf} is allowed.
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
X <- rt(1000,15)
test_ks_dt(X, df=15)
test_ks_dt(X, df=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
