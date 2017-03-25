

\name{rsde2d}
\alias{rsde2d}
\alias{rsde2d.default}
\alias{dsde2d}
\alias{dsde2d.default}
\alias{print.dsde2d}
\alias{plot.dsde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Density and random generation for conditional law of 2-Dim SDE's}
\description{
Density and random generation for the joint and marginal law of \code{(X(t),Y(t) | X(0)=x0,Y(0)=y0)} of the SDE's 2-d.
}
\usage{
rsde2d(object, \dots)
dsde2d(object, ...)
\method{rsde2d}{default}(object, at, \dots)
\method{dsde2d}{default}(object, pdf=c("Joint","Marginal"), at, ...)
\method{plot}{dsde2d}(x,display=c("persp","rgl","image","contour"),hist=FALSE,...)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object inheriting from class \code{\link{snssde2d}}.}
  \item{at}{time between \code{t0} and \code{T}. Random number generated of the process at time \code{at}. The default \code{at = T}.}
  \item{pdf}{probability density function \code{Joint} or \code{Marginal}.}  
  \item{x}{an object inheriting from class \code{dsde2d}.}
  \item{display}{display plots.} 
  \item{hist}{if \code{hist=TRUE} plot histogram.} 
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{rsde2d} returns a random variable \eqn{x_{t=at},y_{t=at}}{x(t=at),y(t=at)} realize at time \eqn{t=at}{t=at}.

And \code{dsde2d} returns a bivariate density approximation for conditional law \code{(X(t),Y(t) | X(0)=x0,Y(0)=y0)}. with \eqn{ at }{at} is a fixed time between \code{t0} and \code{T}.
}
\value{
\code{dsde2d} {gives the bivariate density approximation for conditional law \code{(X(t),Y(t) | X(0)=x0,Y(0)=y0)}.}
\code{rsde2d} {generates random of the couple \code{(X(t),Y(t) | X(0)=x0,Y(0)=y0)}.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link[yuima]{rng}} random number generators in "yuima" package. 

\code{\link[DiffusionRgqd]{BiGQD.density}} Generate the Transition Density of a Bivariate Generalized Quadratic Diffusion Model (2D GQD).
}
\examples{

## 2-dim SDE's
## Example:
## random numbers of two standard Brownian motion W1(t) and W2(t) at time = 0.5
set.seed(1234)

# SDE's 2d
fx <- expression(0,0)
gx <- expression(1,1)
mod2d <- snssde2d(drift=fx,diffusion=gx,M=100)

# random 
r2d <- rsde2d(mod2d,at=0.5)
summary(r2d)

# Marginal density 

denM <- dsde2d(mod2d,pdf="M")
denM
plot(denM)

# Joint density
denJ <- dsde2d(mod2d,pdf="J")
denJ
plot(denJ)
plot(denJ,display="contour")

## Example 3: Bivariate Transition Density of 2 Brownian motion (W1(t),W2(t)) in [0,1]

\dontrun{
fx <- expression(0,0)
gx <- expression(1,1)
B <- snssde2d(drift=fx,diffusion=gx,M=1000)
B
for (i in seq(B$Dt,B$T,by=B$Dt)){
plot(dsde2d(B, at = i),display="contour",main=paste0('Transition Density \n t = ',i))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
\keyword{random generators}
