\name{BMcov}
\alias{BMcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Empirical Covariance for Brownian Motion
}
\description{
Calculate empirical covariance of the Brownian Motion. 
}
\usage{
BMcov(N, M, T, C)
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
The brownian motion is a process with increase independent of function the covariance \code{cov(BM) = C * min(t,s)}, If \code{t > s} than \code{cov(BM) = C * s} else \code{cov(BM) = C * t}.
}
\value{
contour of the empirical covariance for brownian motion.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{BMN}} simulation brownian motion by the Normal Distribution , \code{\link{BMRW}} simulation brownian motion by a Random Walk, \code{\link{BMinf}} brownian motion property(Time tends towards the infinite), \code{\link{BMIrt}} brownian motion property(invariance by reversal of time), \code{\link{BMscal}} brownian motion property (invariance by scaling).
}
\examples{
## empirical covariance of 200 trajectories brownian standard
BMcov(N=100,M=250,T=1,C=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
