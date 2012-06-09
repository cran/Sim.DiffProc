\name{test_ks_df}
\alias{test_ks_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kolmogorov-Smirnov Tests (F Distribution)
}
\description{
Performs one sample Kolmogorov-Smirnov tests.
}
\usage{
test_ks_df(X, df1, df2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric vector of data values.
}
  \item{df1}{
degrees of freedom. Inf is allowed.
}
  \item{df2}{
degrees of freedom. Inf is allowed.
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
 X <- rf(1000,10,20)
 test_ks_df(X, df1=10, df2=20)
 test_ks_df(X, df1=5, df2=15)
 test_ks_df(X, df1=15, df2=25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Parametric Estimation }
\keyword{ Statistical Analysis }
