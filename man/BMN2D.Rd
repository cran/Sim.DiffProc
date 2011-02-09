\name{BMN2D}
\alias{BMN2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation Two-Dimensional Brownian Motion (by the Normal Distribution)
}
\description{
simulation 2-dimensional brownian otion in plane (O,X,Y).
}
\usage{
BMN2D(N, t0, T, x0, y0, Sigma, Step = FALSE, Output = FALSE)
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
  \item{x0}{
initial value of \code{BM1(t)} at time \code{t0}.
}
  \item{y0}{
initial value of \code{BM2(t)} at time \code{t0}.
}
  \item{Sigma}{
constant positive.
}
  \item{Step}{
if \code{Step = TRUE} ploting step by step.
}
  \item{Output}{
if \code{output = TRUE} write a \code{output} to an Excel 2007.
}
}
\details{
see , \code{\link{BMN}}
}
\value{
data.frame(time,W1(t),W2(t)) and plot of process 2-D.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\seealso{
\code{\link{BMN3D}} Simulation Three-Dimensional Brownian Motion.
}
\examples{
BMN2D(N=5000, t0=0, T=1, x0=0, y0=0,Sigma=0.2, 
     Step = FALSE, Output = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Diffusion Process Multidimensional}
