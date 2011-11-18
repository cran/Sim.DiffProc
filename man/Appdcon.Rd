\name{Appdcon}
\alias{Appdcon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximated Conditional Law a Diffusion Process
}
\description{
Approximated Conditional densities for \code{X(t)|X(t0)=X0} of a diffusion process.
}
\usage{
Appdcon(x, t, x0, t0, drift, diff, Output = FALSE,
         Methods = c("Euler", "Shoji-Ozaki", "Kessler"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of quantiles.
}
  \item{t}{
calcul at time \code{t}, or evolution in vector times.
}
  \item{x0}{
initial value of the process at time \code{t0}.
}
  \item{t0}{
initial time.
}
  \item{drift}{
drift coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{diff}{
diffusion coefficient: an expression of two variables \code{t} and \code{x}.
}
  \item{Output}{
if \code{Output = TRUE} write a \code{Output} to an Excel 2007.
}
  \item{Methods}{
Approximated methods,\code{see details}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
This function returns the value of the conditional density of \code{X(t)|X(t0)=X0} at point \code{x}.
}
\value{
data.frame(time,f(x(t)|x0)) at final time, and plot of evolution conditional Law.
}
\author{
boukhetala Kamal, guidoum Arsalane.
}
\examples{
## Euler methods
f <- expression(1*(1-x))
g <- expression(0.3)
Appdcon(x=seq(0,3,by=0.01), t = 2 , x0 = 1,t0=0, drift=f, diff=g)
## Kessler s'methods
f <- expression(1*(3-x))
g <- expression(x)
Appdcon(x=seq(0,5,by=0.01), t = seq(0,0.5,by=0.001) , 
       x0 = 1,t0=0, drift=f, diff=g,Methods="Kessler")

## Shoji-Ozaki methods
f <- expression(4*x)
g <- expression(0.3*x)
Appdcon(x=seq(0,3,by=0.01), t = seq(0,1.5,by=0.010) , 
        x0 = 1,t0=0, drift=f, diff=g,Methods="Shoji-Ozaki")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diffusion Process }
\keyword{ Environment R }
\keyword{ Simulation }
\keyword{ Conditional Law }
