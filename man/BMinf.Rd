\name{BMinf}
\alias{BMinf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Brownian Motion Property 
}
\description{
Calculated the limit of standard brownian motion \code{limit(W(t)/t,0,T)}. 
}
\usage{
BMinf(N,T)
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
Calculated the limit of standard brownian motion if the time tends towards the infinite,i.e the \code{limit(W(t)/t,0,T) = 0} .
}
\value{
plot of \code{limit(W(t)/t)}.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{BMN}} simulation brownian motion by the Normal Distribution , \code{\link{BMRW}} simulation brownian motion by a Random Walk, \code{\link{BMIrt}} brownian motion property(invariance by reversal of time), \code{\link{BMscal}} brownian motion property (invariance by scaling), \code{\link{BMcov}} empirical covariance for brownian motion.
}
\examples{
\donttest{
BMinf(N=1000,T=10^5)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
