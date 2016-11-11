\name{rsde1d}
\alias{rsde1d}
\alias{rsde1d.default}
\alias{dsde1d}
\alias{dsde1d.default}
\alias{print.dsde1d}
\alias{plot.dsde1d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Density and random generation for conditional law of 1-Dim SDE
}
\description{
Density and random generation for the conditional law of \code{X(t) | X(0)=x0} of the 1-dim SDE.
}
\usage{
rsde1d(object, ...)
dsde1d(object, ...)
\method{rsde1d}{default}(object, at, ...)
\method{dsde1d}{default}(object, at, ...)
\method{plot}{dsde1d}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object inheriting from class \code{\link{snssde1d}}.
}
  \item{at}{
time between \code{t0} and \code{T}. The default \code{at = T}.
}
  \item{x}{
an object inheriting from class \code{dsde1d}.
}
  \item{\dots}{
further arguments for (non-default) methods.
}
}
\details{
The function \code{rsde1d} returns a random variable \eqn{x_{t=at}}{x(t=at)} realize at time \eqn{t=at}{t=at} defined by : 
\deqn{x_{ t=at } = \{ t \geq 0 ; x = X_{ t=at } \}}{x(t=at) = { t>=0 ; x = X(t=at) }}

And \code{dsde1d} returns a kernel density approximation for conditional law \code{X(t) | X(0)=x0}.
with \eqn{ at }{at} is a fixed time between \code{t0} and \code{T}.
}
\value{
\code{dsde1d} {gives the density estimate of \code{X(t) | X(0)=x0}.}
\code{rsde1d} {generates random of \code{X(t) | X(t)=x0}.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}
\seealso{
\code{\link[yuima]{rng}} random number generators in "yuima" package. 

\code{\link[sde]{dcSim}} Pedersen's simulated transition density in "sde" package.

\code{\link[sde]{rcBS}}, \code{\link[sde]{rcCIR}}, \code{\link[sde]{rcOU}} and \code{\link[sde]{rsOU}}  in package "sde".

\code{\link[sde]{dcBS}}, \code{\link[sde]{dcCIR}}, \code{\link[sde]{dcOU}} and \code{\link[sde]{dsOU}}  in package "sde". 
}
\examples{

## Example 1:  
## dX(t) = (-2*(X(t)<=0)+2*(X(t)>=0)) *dt + 0.5 o dW(t)
set.seed(1234)

f <- expression(-2*(x<=0)+2*(x>=0))
g <- expression(0.5)
res1 <- snssde1d(drift=f,diffusion=g,M=100,type="str")
x <- rsde1d(res1, at = 1)
dens1 <-  dsde1d(res1, at = 1)
dens1
plot(dens1)  # kernel estimated
plot(dens1,hist=TRUE,h=0.7) # histogramme

## Example 2:
## random numbers of standard Brownian motion W(t) at time = 0.5
set.seed(1234)

f <- expression(0)
g <- expression(1)
res2 <- snssde1d(drift=f,diffusion=g,M=100)
x <- rsde1d(res2, at = 0.5)
plot(dsde1d(res2, at = 0.5),dens=function(x) dnorm(x,0,sqrt(0.5)))
plot(dsde1d(res2, at = 0.5),dens=function(x) dnorm(x,0,sqrt(0.5)),hist=TRUE)

## Example 3: Transition Density of Brownian motion W(t) in [0,1]

\dontrun{
f <- expression(0); g <- expression(1)
B <- snssde1d(drift=f,diffusion=g,M=1000)
for (i in seq(B$t0,B$T,by=B$Dt)){
plot(dsde1d(B, at = i),main=paste0('Transition Density \n t = ',i))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
\keyword{random generators}
