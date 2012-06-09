\name{MartExp}
\alias{MartExp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating The Exponential Martingales Process
}
\description{
Simulation the exponential martingales. 
}
\usage{
MartExp(N, t0, T, sigma, output = FALSE)
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
  \item{sigma}{
constant positive (\code{sigma is volatility}).
}
  \item{output}{
if \code{output = TRUE} write a \code{output} to an Excel (.csv).
}
}
\details{
That is to say \code{W(t)} a Brownian movement the following processes are continuous martingales :
\enumerate{\item \code{X(t) = W(t)^2 - t}. \item \code{Y(t) = exp( integral(f(s)dW(s),0,t)- 0.5 * integral(f(s)^2 ds,0,t) )}.}
}
\value{
data.frame(time,x,y) and plot of process.
}
\author{
Boukhetala Kamal, Guidoum Arsalane.
}
\examples{
\donttest{
## Exponential Martingales Process
 MartExp(N=1000,t0=0,T=1,sigma=2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Environment R }
\keyword{ Simulation }
