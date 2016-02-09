\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{rsde3d}
\alias{rsde3d}
\alias{rsde3d.default}
\alias{summary.rsde3d}
\alias{mean.rsde3d}
\alias{median.rsde3d}
\alias{quantile.rsde3d}
\alias{kurtosis.rsde3d}
\alias{skewness.rsde3d}
\alias{moment.rsde3d}
\alias{bconfint.rsde3d}
\alias{plot.rsde3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random Number Generators for 3-Dim SDE}
\description{The (S3) generic function \code{rsde3d} for simulate random number generators to generate 3-dim sde.}
\usage{
rsde3d(N, \dots)
\method{rsde3d}{default}(N = 1000, M = 100, x0 = 0, y0 = 0, z0 = 0, t0 = 0, T = 1, Dt, tau = 0.5,
   driftx, diffx, drifty, diffy, driftz, diffz, alpha = 0.5, mu = 0.5, 
   type = c("ito", "str"), method = c("euler", "milstein", "predcorr",
   "smilstein", "taylor", "heun", "rk1", "rk2", "rk3"), \dots)		

   
\method{summary}{rsde3d}(object, \dots)
\method{mean}{rsde3d}(x, \dots)
\method{median}{rsde3d}(x, \dots)
\method{quantile}{rsde3d}(x, \dots)
\method{kurtosis}{rsde3d}(x, \dots)
\method{skewness}{rsde3d}(x, \dots)
\method{moment}{rsde3d}(x, order = 2, \dots)
\method{bconfint}{rsde3d}(x, level=0.95, \dots)
\method{plot}{rsde3d}(x, \dots)
}
%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{N}{size of sde.}
  \item{M}{number of random numbers to be geneated.}
  \item{x0, y0, z0}{initial value of the process \eqn{X_{t}}{X(t)}, \eqn{Y_{t}}{Y(t)} and \eqn{Z_{t}}{Z(t)} at time \code{t0}.}
  \item{t0}{initial time.}
  \item{T}{final time.} 
  \item{Dt}{time step of the simulation (discretization). If it is \code{\link[base]{missing}} a default \eqn{\Delta t = \frac{T-t_{0}}{N}}{Dt = (T-t0)/N}.}
  \item{tau}{moment (time) between \code{t0} and \code{T}. Random number generated at \code{time=tau}.}
  \item{driftx, drifty, driftz}{drift coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}, \eqn{Y_t}{Y(t)} and \eqn{Z_t}{Z(t)}.}
  \item{diffx, diffy, diffz}{diffusion coefficient: an \code{\link[base]{expression}} of four variables \code{t}, \code{x}, \code{y} and \code{z} for process \eqn{X_t}{X(t)}, \eqn{Y_t}{Y(t)} and \eqn{Z_t}{Z(t)}.}
  \item{alpha, mu }{weight of the predictor-corrector scheme; the default \code{alpha = 0.5} and \code{mu = 0.5}.}
  \item{type}{sde of  the type Ito or Stratonovich.}
  \item{method}{numerical methods of simulation, the default \code{method = "euler"}; see \code{\link{snssde3d}}.}
  \item{x, object}{an object inheriting from class \code{"rsde2d"}.}
  \item{order}{order of moment.}
  \item{level}{the confidence level required.}    
  \item{\dots}{further arguments for (non-default) methods.}
}
\details{
The function \code{rsde3d} returns a three random variables \eqn{(x_{\tau},y_{\tau},z_{\tau})}{(x(tau),y(tau),z(tau))} realize at time \eqn{t=\tau}{t=tau} defined by :  
\deqn{x_{\tau} = \{ t \geq 0 ; x = X_{\tau} \}}{x(tau)={t>=0; x = X(tau)}}
\deqn{y_{\tau} = \{ t \geq 0 ; y = Y_{\tau} \}}{y(tau)={t>=0; y = Y(tau)}}
\deqn{z_{\tau} = \{ t \geq 0 ; z = Z_{\tau} \}}{z(tau)={t>=0; z = Z(tau)}}
with \eqn{\tau}{tau} is a fixed time between \code{t0} and \code{T}.
}
\value{
\code{rsde3d} returns an object inheriting from \code{\link[base]{class}} \code{"rsde3d"}.
\item{x, y, z }{a vector of random numbers of 3-dim sde realize at time time \eqn{t=\tau}{t=tau}, the triplet \eqn{(x_{\tau},y_{\tau},z_{\tau})}{(x(tau),y(tau),z(tau))}.}
}
\author{
A.C. Guidoum, K. Boukhetala.
}

\seealso{
\code{\link{rsde1d}} simulation RNs in sde 1-dim.

\code{\link[yuima]{rng}} random number generators in \CRANpkg{yuima} package. 

\code{\link[sde]{rcBS}}, \code{\link[sde]{rcCIR}}, \code{\link[sde]{rcOU}} and \code{\link[sde]{rsOU}}  in package \CRANpkg{sde}. 
}
\examples{
## Example 1: Ito sde 3-dim
## dX(t) = 4*(-1-X(t))*Y(t) dt + 0.2 * dW1(t) 
## dY(t) = 4*(1-Y(t)) *X(t) dt + 0.2 * dW2(t) 
## dZ(t) = 4*(1-Z(t)) *Y(t) dt + 0.2 * dW3(t)       
## W1(t), W2(t) and W3(t) three independent Brownian motion  
set.seed(1234)

fx <- expression(4*(-1-x)*y)
gx <- expression(0.2)
fy <- expression(4*(1-y)*x)
gy <- expression(0.2)
fz <- expression(4*(1-z)*y)
gz <- expression(0.2)

res <- rsde3d(driftx=fx,diffx=gx,drifty=fy,diffy=gy,driftz=fz,diffz=gz,
               x0=2,y0=-2,z0=0,tau=0.3,M=50)
res
summary(res)
plot(res,union=TRUE)
dev.new()
plot(res,union=FALSE)
X <- cbind(res$x,res$y,res$z)
## library(sm)
## sm.density(X,display="rgl")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sde}
\keyword{ts}
\keyword{mts}
\keyword{random generators}
