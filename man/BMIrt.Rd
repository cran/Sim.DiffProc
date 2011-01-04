\name{BMIrt}
\alias{BMIrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Brownian Motion Property (Invariance by reversal of time)
}
\description{
Brownian motion is invariance by reversal of time.
}
\usage{
BMIrt(N, T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{T}{
final time.
}
}
\details{
Brownian motion is invariance by reversal of time,i.e \code{W(t) = W(T-t) - W(T)}. 
}
\value{
plot of \code{W(T-t) - W(T)}. 
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{BMN}} simulation brownian motion by the Normal Distribution , \code{\link{BMRW}} simulation brownian motion by a Random Walk, \code{\link{BMinf}} Brownian Motion Property (time tends towards the infinite), \code{\link{BMscal}} brownian motion property (invariance by scaling), \code{\link{BMcov}} empirical covariance for brownian motion.
}
\examples{
BMIrt(N=1000,T=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
