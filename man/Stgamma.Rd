\name{Stgamma}
\alias{Stgamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Stochastic Process The Gamma Distribution
}
\description{
Simulation stochastic process by a gamma distribution.
}
\usage{
Stgamma(N, t0, T, alpha, beta, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{t0}{
initial time.
}
  \item{T}{
final time.
}
  \item{alpha}{
constant positive.
}
  \item{beta}{
an alternative way to specify the scale.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\value{
data.frame(time,x) and plot of process.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{SRW}} Creating Random Walk, \code{\link{Stst}} Stochastic Process The Student Distribution, \code{\link{WNG}} White Noise Gaussian.
}
\examples{
## Stochastic Process The Gamma Distribution
Stgamma(N=1000,t0=0,T=5,alpha=1,beta=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
