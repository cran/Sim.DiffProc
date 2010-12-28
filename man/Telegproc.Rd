\name{Telegproc}
\alias{Telegproc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Realization a Telegraphic Process 
}
\description{
Simulation a telegraphic process.
}
\usage{
Telegproc(t0, x0, T, lambda, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t0}{
initial time.
}
  \item{x0}{
state initial \code{(x0 = -1 or +1)}.
}
  \item{T}{
final time of the simulation.
}
  \item{lambda}{
exponential distribution with rate \code{lambda}.
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{Asys}} Evolution a Telegraphic Process. 
}
\examples{
## Simulation a telegraphic process
Telegproc(t0=0,x0=1,T=1,lambda=0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
