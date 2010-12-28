\name{BMRWF}
\alias{BMRWF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Flow of Brownian Motion (by a Random Walk)
}
\description{
Simulation flow of the brownian motion model by a Random Walk.
}
\usage{
BMRWF(N, M, t0, T, C, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{M}{
number of trajectories.
}
  \item{t0}{
initial time.
}
  \item{T}{
final time.
}
  \item{C}{
constant positive (\code{if C = 1 it is standard brownian motion}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
One characterization of the Brownian motion says that it can be seen as the limit of a random walk in the following sense.

Given a sequence of independent and identically distributed random variables \code{X1, X2, . . . , Xn}, taking only two values \code{+1} and \code{-1} with equal probability and considering the partial \code{sum}, \code{Sn = X1+ X2+ . . . + Xn}. then, as \code{n --> lnf},\code{P(Sn/sqrt(N) < x) = P( W(t) < x)}.

Where \code{[x]} is the integer part of the real number \code{x}. Please note that this result is a refinement of the central limit theorem that, in our case, asserts that \code{Sn/sqrt(n) ~~> N(0,1)}.
}
\value{
data.frame(time,x) and plot of process. 
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{BMN}} simulation brownian motion by the normal distribution, \code{\link{BMRW}} simulation brownian motion by a random walk, \code{\link{BB}} Simulation of brownian bridge model, \code{\link{GBM}} simulation geometric brownian motion Model.
}
\examples{
##
BMRWF(N=1000,M=5,t0=0,T=1,C=1)
BMRWF(N=1000,M=5,t0=0,T=1,C=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
