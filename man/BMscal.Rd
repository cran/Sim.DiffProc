\name{BMscal}
\alias{BMscal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Brownian Motion Property (Invariance by scaling)
}
\description{
Brownian motion with different scales.
}
\usage{
BMscal(N, T, S1, S2, S3, output = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
size of process.
}
  \item{T}{
final time.
}
  \item{S1}{
constant (\code{scale 1}).
}
  \item{S2}{
constant (\code{scale 2}).
}
  \item{S3}{
constant (\code{scale 3}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
Brownian motion is invariance by change the scales,i.e \code{W(t) = (1/S) * W(S^2 * t)}, \code{S} is scale.
}
\value{
data.frame(w1,w2,w3) and plot of process. 
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\seealso{
\code{\link{BMinf}} brownian motion Property (time tends towards the infinite), \code{\link{BMcov}} empirical covariance for brownian motion, \code{\link{BMIrt}} brownian motion property(invariance by reversal of time).
}
\examples{
\donttest{
 BMscal(N=1000,T=10,S1=1,S2=1.1,S3=1.2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
