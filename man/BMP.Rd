\name{BMP}
\alias{BMP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Brownian Motion Property (trajectories brownian between function (+/-)2*sqrt(C*t))
}
\description{
trajectories Brownian lies between the two curves \code{(+/-)2*sqrt(C*t)}.
}
\usage{
BMP(N, M, T, C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{M}{
number of trajectories.
}
  \item{T}{
final time.
}
  \item{C}{
constant positive (\code{if C = 1 it is standard brownian motion}).
}
}
\details{
A flow of brownian motion lies between the two curves \code{(+/-)2*sqrt(C*t)}, \code{W(dt) - W(0) ~~> N(0,dt)}, \code{N(0,dt)} normal distribution. 
}
\value{
plot of the flow.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{BMscal}} brownian motion property (invariance by scaling), \code{\link{BMinf}} brownian motion Property (time tends towards the infinite), \code{\link{BMcov}} empirical covariance for brownian motion, \code{\link{BMIrt}} brownian motion property(invariance by reversal of time).
}
\examples{
 BMP(N=1000,M=10,T=1,C=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
