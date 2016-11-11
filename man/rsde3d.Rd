

\name{rsde3d}
\alias{rsde3d}
\alias{rsde3d.default}
\alias{dsde3d}
\alias{dsde3d.default}
\alias{print.dsde3d}
\alias{plot.dsde3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density and random generation for conditional law of 3-Dim SDE's}
\description{
Density and random generation for the marginal conditional law of \code{(X(t),Y(t),Z(t) | X(0)=x0,Y(0)=y0,Z(0)=z0)} of the SDE's 3-d.
}
\usage{
rsde3d(object, \dots)
dsde3d(object, ...)
\method{rsde3d}{default}(object, at, \dots)
\method{dsde3d}{default}(object, at, ...)
\method{plot}{dsde3d}(x,hist=FALSE,...)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde3d}}.}
  \item{at}{time between \code{t0} and \code{T}. Random number generated of the process at time \code{at}. The default \code{at = T}.}
  \item{x}{an object inheriting from class \code{dsde3d}.}
  \item{hist}{if \code{hist=TRUE} plot histogram.} 
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{rsde3d} returns a random variable \eqn{x_{t=at},y_{t=at},z_{t=at}}{x(t=at),y(t=at),z(t=at)} realize at time \eqn{t=at}{t=at}.

And \code{dsde3d} returns a kernel marginal density approximation for conditional law \code{(X(t),Y(t),Z(t) | X(0)=x0,Y(0)=y0,Z(0)=z0)}. with \eqn{ at }{at} is a fixed time between \code{t0} and \code{T}.
}
\value{
\code{dsde3d} {gives the kernel marginal density approximation for conditional law of \code{(X(t),Y(t),Z(t) | X(0)=x0,Y(0)=y0,Z(0)=z0)}.}
\code{rsde3d} {generates random of the couple \code{(X(t),Y(t),Z(t) | X(0)=x0,Y(0)=y0,Z(0)=z0)}.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link[yuima]{rng}} random number generators in "yuima" package. 

\code{\link[sde]{rcBS}}, \code{\link[sde]{rcCIR}}, \code{\link[sde]{rcOU}} and \code{\link[sde]{rsOU}}  in package "sde". 
}
\examples{

## 3-dim SDE's

## Example 3: Ito sde 3-dim
## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 * dW1(t) 
## dY(t) = 4*(1-Y(t)) *X(t) dt + 0.2 * dW2(t) 
## dZ(t) = 4*(1-Z(t)) *Y(t) dt + 0.2 * dW3(t)       
## W1(t), W2(t) and W3(t) three independent Brownian motion  
set.seed(1234)

# SDE's 3d
fx <- expression(4*(-1-x)*y, 4*(1-y)*x, 4*(1-z)*y)
gx <- rep(expression(0.2),3)
mod3d <- snssde3d(drift=fx,diffusion=gx,x0=c(2,-2,0),M=25)

# random
r3d <- rsde3d(mod3d)
summary(r3d)

# Marginal density

denM <- dsde3d(mod3d)
denM
plot(denM)

## for Joint density see package 'sm' or 'ks'

# library(sm)
# sm.density(Data3d,display="rgl")

# library(ks)
# fhat <- kde(x=Data3d )
# plot(fhat, drawpoints=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
\keyword{random generators}
