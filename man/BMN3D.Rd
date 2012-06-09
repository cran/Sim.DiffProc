\name{BMN3D}
\alias{BMN3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation Three-Dimensional Brownian Motion (by the Normal Distribution)
}
\description{
simulation 3-dimensional brownian otion in (O,X,Y,Z).
}
\usage{
BMN3D(N, t0, T, X0, Y0, Z0, Sigma, Output = FALSE)
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
  \item{X0}{
initial value of \code{BM1(t)} at time \code{t0}.
}
  \item{Y0}{
initial value of \code{BM2(t)} at time \code{t0}.
}
  \item{Z0}{
initial value of \code{BM3(t)} at time \code{t0}.
}
  \item{Sigma}{
constant positive.
}
  \item{Output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
see , \code{\link{BMN}}
}
\value{
data.frame(time,W1(t),W2(t),W3(t)) and plot of process 3-D.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{BMRW3D}} Simulation Three-Dimensional Brownian Motion.
}
\examples{
\donttest{
 BMN3D(N=500, t0=0, T=1, X0=0.5, Y0=0.5, Z0=0.5, 
      Sigma=0.3, Output = FALSE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Diffusion Process Multidimensional}
