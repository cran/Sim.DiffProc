\name{BMNF}
\alias{BMNF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Flow of Brownian Motion (by the Normal Distribution)
}
\description{
Simulation flow of the brownian motion model by the normal distribution.
}
\usage{
BMNF(N, M, t0, T, C, output = FALSE)
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
Given a fixed time increment \code{dt = (T-t0)/N}, one can easily simulate a flow of the Wiener process in the time interval \code{[t0,T]}. Indeed, for \code{W(dt)} it holds true that \code{W(dt) = W(dt) - W(0) ~> N(0,dt) ~> sqrt(dt) * N(0,1)}, \code{N(0,1)} normal distribution.
}
\value{
data.frame(time,x) and plot of process. 
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{BMRW}} simulation brownian motion by a random walk, \code{\link{BMN}} simulation of brownian motion by the normal distribution, \code{\link{BMRWF}} simulation flow of brownian motion by a random walk, \code{\link{BB}} Simulation of brownian bridge model, \code{\link{GBM}} simulation geometric brownian motion Model.
}
\examples{
##
BMNF(N=1000,M=5,t0=0,T=1,C=1)
BMNF(N=1000,M=5,t0=0,T=1,C=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ financial models }
