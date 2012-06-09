\name{SRW}
\alias{SRW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating Random Walk
}
\description{
Simulation random walk.
}
\usage{
SRW(N, t0, T, p, output = FALSE)
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
  \item{p}{
probability of choosing \code{X = -1 or +1}.
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
\code{\link{Stgamma}} Stochastic Process The Gamma Distribution, \code{\link{Stst}} Stochastic Process The Student Distribution, \code{\link{WNG}} White Noise Gaussian.
}
\examples{
## Random Walk
\donttest{
 SRW(N=1000,t0=0,T=1,p=0.5)
 SRW(N=1000,t0=0,T=1,p=0.25)
 SRW(N=1000,t0=0,T=1,p=0.75)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
