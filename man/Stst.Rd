\name{Stst}
\alias{Stst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Stochastic Process The Student Distribution
}
\description{
Simulation stochastic process by a Student distribution.
}
\usage{
Stst(N, t0, T, n, output = FALSE)
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
  \item{n}{
degrees of freedom (\code{> 0,non-integer}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\value{
data.frame(time,x) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{SRW}} Creating Random Walk, \code{\link{Stgamma}} Stochastic Process The Gamma Distribution, \code{\link{WNG}} White Noise Gaussian.
}
\examples{
\donttest{
## Stochastic Process The Student Distribution
Stst(N=1000,t0=0,T=1,n=2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
