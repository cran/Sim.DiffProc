

\name{rsde2d}
\alias{rsde2d}
\alias{rsde2d.default}
\alias{summary.rsde2d}
\alias{mean.rsde2d}
\alias{median.rsde2d}
\alias{quantile.rsde2d}
\alias{kurtosis.rsde2d}
\alias{skewness.rsde2d}
\alias{moment.rsde2d}
\alias{bconfint.rsde2d}
\alias{plot.rsde2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random Number Generators for 2-Dim SDE}
\description{The (S3) generic function \code{rsde2d} for simulate random number generators to generate 2-dim sde.}
\usage{
rsde2d(N, \dots)
\method{rsde2d}{default}(N = 1000, M = 100, x0 = 0, y0 = 0, t0 = 0, T = 1, Dt, tau = 0.5, 
   driftx, diffx, drifty, diffy, alpha = 0.5, mu = 0.5, type = c("ito", "str"), 
   method = c("euler", "milstein", "predcorr", "smilstein", "taylor", 
   "heun", "rk1", "rk2", "rk3"), \dots)				 
   
   
\method{summary}{rsde2d}(object, \dots)
\method{mean}{rsde2d}(x, \dots)
\method{median}{rsde2d}(x, \dots)
\method{quantile}{rsde2d}(x, \dots)
\method{kurtosis}{rsde2d}(x, \dots)
\method{skewness}{rsde2d}(x, \dots)
\method{moment}{rsde2d}(x, order = 2, \dots)
\method{bconfint}{rsde2d}(x, level=0.95, \dots)
\method{plot}{rsde2d}(x, \dots)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{size of sde.}
  \item{M}{number of random numbers to be geneated.}
  \item{x0, y0}{initial value of the process \eqn{X_{t}}{X(t)} and \eqn{Y_{t}}{Y(t)} at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{tau}{time between \code{t0} and \code{T}. Random number generated at \code{time=tau}.}
  \item{driftx, drifty}{drift coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)}.}
  \item{diffx, diffy}{diffusion coefficient: an \code{\link[base]{expression}} of three variables \code{t}, \code{x} and \code{y} for process \eqn{X_t}{X(t)} and \eqn{Y_t}{Y(t)}.}
  \item{alpha, mu}{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{sde of  the type Ito or Stratonovich.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde2d}}.}
  \item{x, object}{an object inheriting from class \code{"rsde2d"}.}  
  \item{order}{order of moment.}
  \item{level}{the confidence level required.}    
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{rsde2d} returns a two random variables \eqn{(x_{\tau},y_{\tau})}{(x(tau),y(tau))} realize at time \eqn{t=\tau}{t=tau} defined by : 
\deqn{x_{\tau} = \{ t \geq 0 ; x = X_{\tau} \}}{x(tau)={t>=0; x = X(tau)}}
\deqn{y_{\tau} = \{ t \geq 0 ; y = Y_{\tau} \}}{y(tau)={t>=0; y = Y(tau)}}
with \eqn{\tau}{tau} is a fixed time between \code{t0} and \code{T}.
}
\value{
\code{rsde2d} returns an object inheriting from \code{\link[base]{class}} \code{"rsde2d"}.
\item{x, y }{a vector of random numbers of 2-dim sde realize at time time \eqn{t=\tau}{t=tau}, the couple \eqn{(x_{\tau},y_{\tau})}{(x(tau),y(tau))}.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link{rsde1d}} simulation RNs in sde 1-dim.

\code{\link[yuima]{rng}} random number generators in "yuima" package. 

\code{\link[sde]{rcBS}}, \code{\link[sde]{rcCIR}}, \code{\link[sde]{rcOU}} and \code{\link[sde]{rsOU}}  in package "sde". 
}
\examples{

## Example 1:
## random numbers of two standard Brownian motion W1(t) and W2(t) at time = 1
set.seed(1234)

fx <- expression(0)
gx <- expression(1)
fy <- expression(0)
gy <- expression(1)
res1 <- rsde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,tau=1)
res1
summary(res1)
X <- cbind(res1$x,res1$y)
## library(sm)
## sm.density(X,display="persp")

## Example 2: 
## dX(t) = 5*(-1-Y(t))*X(t) * dt + 0.5 * dW1(t)          
## dY(t) = 5*(-1-X(t))*Y(t) * dt + 0.5 * dW2(t)
## W1(t) and W2(t) two independent Brownian motion
set.seed(1234)

fx <- expression(5*(-1-y)*x)
gx <- expression(0.5)
fy <- expression(5*(-1-x)*y)
gy <- expression(0.5)
res2 <- rsde2d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,tau=0.4876
               ,x0=2,y0=-2,M=50)
res2
summary(res2)
plot(res2,union=TRUE)
dev.new()
plot(res2,union=FALSE)
X <- cbind(res2$x,res2$y)
## sm.density(X,display="persp")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
\keyword{random generators}
